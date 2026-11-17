#!/bin/bash

# Group statistics on per-subject ROI tables (ASE_groupROIStats.m), then optional painting of the
# results onto template(s) (paintStats.sh --table). Replaces batch_groupStats.csh (voxelwise
# 3dttest++ on painted maps): the model is fit once per ROI, which is what the painted maps
# amounted to, but with missing ROIs left out rather than counted as 0 and with FDR across ROIs.
#
# Subjects come from the covariate table: every uncommented row passing --select. The data list
# is not used. Selected subjects without a usable fit, and subject folders with a fit but no row
# in the table, are reported (group_report.txt).
#
# Output: <outPath>/<timestamp>[_<fitTag>]/ group_ROI_Stats.csv, group_ROI_long.csv,
#         group_subjects.csv, group_report.txt (+ painted_<column>.nii.gz and painted_paint_report.txt
#         with --paint).

need_arg() {
	if [[ $# -lt 2 ]]; then
		echo "Error: $1 requires an argument" >&2
		exit 1
	fi
}
make_cell() {
	local out="{" first=1 item
	for item in "$@"; do
		((first)) || out+=","
		out+="'${item//\'/\'\'}'"
		first=0
	done
	echo "${out}}"
}
mq() { printf '%s' "${1//\'/\'\'}"; }
make_str() {   # a b c -> ["a","b","c"]  (MATLAB string array; empty -> strings(0,1))
	local out="[" first=1 item
	(($# == 0)) && { echo "strings(0,1)"; return; }
	for item in "$@"; do
		((first)) || out+=","
		out+="\"${item//\"/\"\"}\""
		first=0
	done
	echo "${out}]"
}

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repoDir="$(cd "$scriptDir/.." && pwd)"

allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
covTable=""                    # default: <allSubjDir>/cov_table_ase
outPath=""                     # default: <allSubjDir>/group_stats
selects=()
covariates="age"
categorical=""
stats="OEF,R2P,vCBV"
fitTag=""
fitsRelDir="session_ase/fits"
tableName="ase_fit_results_ROI_Stats.csv"
minVox=0
minN=5
settingsCheck="error"
paint=0
paintStatsList=""              # default: <stat>_adjMean and <stat>_<cov>_t for every stat/covariate
paintArgs=()
matDir=$repoDir
matlabBin=matlab
selectDefault="group=ctrl,mdt"

usage() {
cat <<EOF
Usage: $0 [options]

Subjects and model:
  --allSubjDir <path>     Directory containing the subject folders [$allSubjDir]
  --covTable <file>       Covariate table: header line (may start with #), one row per subject,
                          first column = subject ID; '#' rows are excluded [<allSubjDir>/cov_table_ase]
  --select <filter>       Row filter, repeatable (all must hold): col=a,b  col!=a  col>x  col>=x
                          col<x  col<=x   [default if none given: $selectDefault]
  --selectAll             Use every uncommented row (no default filter)
  --covariates <a,b>      Covariate columns ("" = none) [$covariates]
  --categorical <a,b>     Numeric columns to treat as categorical (e.g. sex); text columns always are
  --stats <a,b>           ROI-table columns to analyse [$stats]
  --minVox <n>            Use an ROI value only if its nVox >= n [$minVox]
  --minN <n>              Minimum subjects for a model (also >= coefficients + 2) [$minN]

Fits:
  --fitTag <label>        Use each subject's latest run tagged <label> (aseFit.sh --tag) [latest run]
  --fitsRelDir <path>     Fits folder relative to each subject [$fitsRelDir]
  --tableName <name>      [$tableName]
  --settingsCheck <m>     error | warn | off: fit settings must match across subjects [$settingsCheck]

Output:
  --outPath <path>        Results go to <outPath>/<timestamp>[_<fitTag>] [<allSubjDir>/group_stats]

Painting (optional; runs paintStats.sh --table on the result):
  --paint                 Paint the group table onto the template(s)
  --paintStats <a,b>      Columns to paint [<stat>_adjMean and every <stat>_<term>_t]
  --templateDir, --templateNames, --LUTs, --omitCodes, --omitNames, --fillValue
                          Passed to paintStats.sh (see paintStats.sh --help)

  --matDir <path>         MATLAB code location [$matDir]
  --matlabBin <path>      [$matlabBin]
  -h, --help
EOF
}

selectAll=0
while [[ $# -gt 0 ]]; do
	case "$1" in
		--allSubjDir)    need_arg "$@"; allSubjDir=$2;    shift 2 ;;
		--covTable)      need_arg "$@"; covTable=$2;      shift 2 ;;
		--outPath)       need_arg "$@"; outPath=$2;       shift 2 ;;
		--select)        need_arg "$@"; selects+=("$2");  shift 2 ;;
		--selectAll)     selectAll=1;                     shift ;;
		--covariates)    need_arg "$@"; covariates=$2;    shift 2 ;;
		--categorical)   need_arg "$@"; categorical=$2;   shift 2 ;;
		--stats)         need_arg "$@"; stats=$2;         shift 2 ;;
		--minVox)        need_arg "$@"; minVox=$2;        shift 2 ;;
		--minN)          need_arg "$@"; minN=$2;          shift 2 ;;
		--fitTag)        need_arg "$@"; fitTag=$2;        shift 2 ;;
		--fitsRelDir)    need_arg "$@"; fitsRelDir=$2;    shift 2 ;;
		--tableName)     need_arg "$@"; tableName=$2;     shift 2 ;;
		--settingsCheck) need_arg "$@"; settingsCheck=$2; shift 2 ;;
		--paint)         paint=1;                         shift ;;
		--paintStats)    need_arg "$@"; paintStatsList=$2; shift 2 ;;
		--templateDir|--templateNames|--LUTs|--omitCodes|--omitNames|--fillValue)
		                 need_arg "$@"; paintArgs+=("$1" "$2"); shift 2 ;;
		--matDir)        need_arg "$@"; matDir=$2;        shift 2 ;;
		--matlabBin)     need_arg "$@"; matlabBin=$2;     shift 2 ;;
		-h|--help)       usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done

# ---- checks ----
[[ -d $allSubjDir ]] || { echo "Error: allSubjDir not found: $allSubjDir" >&2; exit 1; }
allSubjDir="$(cd "$allSubjDir" && pwd)"
[[ -n $covTable ]] || covTable="$allSubjDir/cov_table_ase"
[[ -e $covTable ]] || { echo "Error: covariate table not found: $covTable" >&2; exit 1; }
covTable="$(cd "$(dirname "$covTable")" && pwd)/$(basename "$covTable")"
[[ -n $outPath ]] || outPath="$allSubjDir/group_stats"
if [[ -n $fitTag && ! $fitTag =~ ^[A-Za-z0-9._-]+$ ]]; then
	echo "Error: --fitTag may contain only letters, digits, \".\", \"_\" and \"-\" (got: $fitTag)" >&2; exit 1
fi
for n in "$minVox" "$minN"; do
	[[ $n =~ ^[0-9]+(\.[0-9]+)?$ ]] || { echo "Error: --minVox/--minN must be non-negative numbers (got: $n)" >&2; exit 1; }
done
case $settingsCheck in error|warn|off) ;; *) echo "Error: --settingsCheck must be error, warn or off" >&2; exit 1 ;; esac
if ((${#selects[@]} == 0 && selectAll == 0)); then selects=("$selectDefault"); fi
if ((selectAll)) && ((${#selects[@]} > 0)); then
	echo "Error: --selectAll cannot be combined with --select" >&2; exit 1
fi

IFS=',' read -r -a statArr <<< "$stats"
covArr=(); catArr=()
[[ -n $covariates ]]  && IFS=',' read -r -a covArr <<< "$covariates"
[[ -n $categorical ]] && IFS=',' read -r -a catArr <<< "$categorical"

runName=$(date +%Y%m%d_%H%M%S)
[[ -n $fitTag ]] && runName="${runName}_${fitTag}"
outDir="$outPath/$runName"
mkdir -p "$outDir" || exit 1
echo "Output: $outDir"

tmpDir=$(mktemp -d)
cleanup() {
	[[ -n $tmpDir && -d $tmpDir ]] && rm -rf "$tmpDir"
	rmdir "$outDir" 2>/dev/null && echo "Removed empty output folder: $outDir"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

mfile="$tmpDir/runGroup.m"
cat > "$mfile" <<EOF
addpath(genpath('$(mq "$matDir")'));
try
    ASE_groupROIStats('$(mq "$allSubjDir")', '$(mq "$covTable")', '$(mq "$outDir")', ...
        'select', $(make_str "${selects[@]}"), 'covariates', $(make_str "${covArr[@]}"), ...
        'categorical', $(make_str "${catArr[@]}"), 'stats', $(make_str "${statArr[@]}"), ...
        'fitTag', "$(mq "$fitTag")", 'fitsSubdir', "$(mq "$fitsRelDir")", 'tableName', "$(mq "$tableName")", ...
        'minVox', $minVox, 'minN', $minN, 'settingsCheck', "$settingsCheck");
catch ME
    fprintf(2, 'ERROR: %s\n', ME.message);
    exit(1);
end
EOF
"$matlabBin" -batch "run('$(mq "$mfile")')" || { echo "Group statistics failed." >&2; exit 1; }

if ((paint)); then
	if [[ -z $paintStatsList ]]; then
		# <stat>_adjMean and every model term's t (<stat>_<term>_t), read from the result's header
		cols=()
		IFS=',' read -r -a header < <(head -n 1 "$outDir/group_ROI_Stats.csv" | tr -d '\r"')
		for s in "${statArr[@]}"; do
			for h in "${header[@]}"; do
				[[ $h == "${s}_adjMean" || ( $h == "${s}_"*"_t" ) ]] && cols+=("$h")
			done
		done
		paintStatsList=$(IFS=,; echo "${cols[*]}")
	fi
	echo "Painting: $paintStatsList"
	bash "$scriptDir/paintStats.sh" --table "$outDir/group_ROI_Stats.csv" --outDir "$outDir" \
		--targetStats "$paintStatsList" --outName "painted_%STAT%" --matDir "$matDir" --matlabBin "$matlabBin" \
		"${paintArgs[@]}" || { echo "Painting failed (statistics were written)." >&2; exit 1; }
fi
echo "Done: $outDir"
