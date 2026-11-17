#!/bin/bash

# Paints ROI values onto template segmentation(s) (ASE_paintROIStats.m).
#
# Two modes:
#   per-subject (default): for each subject in the data list (or --singleSubjID), paint the ROI table
#       of the fit run chosen by findFitRun.m: the latest run tagged --fitTag, or the latest run when
#       no tag is given (--allRuns paints every matching run). Outputs go into that run folder.
#   table (--table <csv>): paint one table, e.g. a group_ROI_Stats.csv from aseGroupStats.sh.
#       Outputs go to --outDir (default: the table's folder).
#
# Regions are placed in LUT order and a region that overlaps an earlier one is dropped whole
# (combineSegVols keepOverlaps=false). Use --LUTs to give a painting LUT listing only the regions
# to paint. Each painting writes one report (<outName with %STAT% -> paint_report>.txt).

need_arg() {
	if [[ $# -lt 2 ]]; then
		echo "Error: $1 requires an argument" >&2
		exit 1
	fi
}

# make_cell a b c  ->  {'a','b','c'}   (MATLAB cell array literal)
make_cell() {
	local out="{" first=1 item
	for item in "$@"; do
		((first)) || out+=","
		out+="'${item//\'/\'\'}'"
		first=0
	done
	echo "${out}}"
}
mq() { printf '%s' "${1//\'/\'\'}"; }     # escape single quotes for a MATLAB char literal

tmpDir=""
cleanup() {
	if [[ -n $tmpDir && -d $tmpDir ]]; then rm -rf "$tmpDir"; fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM
scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repoDir="$(cd "$scriptDir/.." && pwd)"

allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
singleSubjID=""
dataList="datalist_ase"
fitsRelDir="session_ase/fits"
tableName="ase_fit_results_ROI_Stats.csv"
fitTag=""
allRuns=0
table=""
outDir=""
templateDir="/autofs/space/guasha_004/users/crbc/yen/templates"
templateNames=""
targetStats="OEF,R2P,vCBV"
LUTs="/autofs/space/guasha_004/users/crbc/yen/ase/GM_WM_long.txt"
matDir=$repoDir
outName="template_painted_%STAT%"
omitCodes=""
omitNames=""
fillValue="NaN"
v=false
matlabBin=matlab

usage() {
    cat << EOF
Usage: $0 [options]
Options:
  --allSubjDir <path>      Directory containing all subjects [$allSubjDir]
  --singleSubjID <id>      Paint one subject only
  --dataList <name>        Subject list in allSubjDir (first column = subject ID; '#' lines skipped) [$dataList]
  --fitsRelDir <path>      Fits folder relative to each subject [$fitsRelDir]
  --tableName <name>       ROI table name inside a fit run [$tableName]
  --fitTag <label>         Use the latest run tagged <label> (aseFit.sh --tag); default: latest run
  --allRuns                Paint every matching run instead of only the latest
  --table <csv>            Table mode: paint this one table (e.g. group_ROI_Stats.csv)
  --outDir <path>          Table mode: output folder [folder of --table]
  --templateDir <path>     Folder with the template segmentation(s) [$templateDir]
  --templateNames <a,b>    Template segmentation file name(s) in templateDir (required)
  --LUTs <a,b>             LUT file path(s) for the template(s); regions are placed in this order [$LUTs]
  --targetStats <a,b>      Table columns to paint [$targetStats]
                           (group tables: e.g. OEF_adjMean,OEF_age_t,OEF_age_q)
  --outName <name>         Output name; %STAT% is replaced by the column name [$outName]
  --omitCodes <a,b>        LUT codes not to paint
  --omitNames <a,b>        LUT names not to paint
  --fillValue <x>          Value outside painted regions [$fillValue]
  --matDir <path>          MATLAB code location [$matDir]
  --matlabBin <path>       MATLAB binary [$matlabBin]
  --verbose -v             Verbose
  -h, --help               Show this help and exit
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --allSubjDir)     need_arg "$@"; allSubjDir="$2"; shift 2;;
        --singleSubjID)   need_arg "$@"; singleSubjID="$2"; shift 2;;
        --dataList)       need_arg "$@"; dataList="$2"; shift 2;;
        --fitsRelDir|--tableRelDir) need_arg "$@"; fitsRelDir="$2"; shift 2;;
        --tableName)      need_arg "$@"; tableName="$2"; shift 2;;
        --fitTag)         need_arg "$@"; fitTag="$2"; shift 2;;
        --allRuns)        allRuns=1; shift;;
        --table)          need_arg "$@"; table="$2"; shift 2;;
        --outDir)         need_arg "$@"; outDir="$2"; shift 2;;
        --templateDir)    need_arg "$@"; templateDir="$2"; shift 2;;
        --templateNames)  need_arg "$@"; templateNames="$2"; shift 2;;
        --LUTs)           need_arg "$@"; LUTs="$2"; shift 2;;
        --targetStats)    need_arg "$@"; targetStats="$2"; shift 2;;
        --outName)        need_arg "$@"; outName="$2"; shift 2;;
        --matDir)         need_arg "$@"; matDir="$2"; shift 2;;
        --omitCodes)      need_arg "$@"; omitCodes="$2"; shift 2;;
        --omitNames)      need_arg "$@"; omitNames="$2"; shift 2;;
        --fillValue)      need_arg "$@"; fillValue="$2"; shift 2;;
        --matlabBin)      need_arg "$@"; matlabBin="$2"; shift 2;;
        --verbose|-v)     v=true; shift ;;
        -h|--help)        usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done

# ---- checks ----
if [[ -z $templateNames ]]; then echo "Error: --templateNames is required" >&2; exit 1; fi
IFS=',' read -r -a templateNameArray <<< "$templateNames"
tmpFiles=()
for name in "${templateNameArray[@]}"; do
	f="$templateDir/$name"
	[[ -e $f ]] || { echo "Error: template not found: $f" >&2; exit 1; }
	tmpFiles+=("$f")
done
IFS=',' read -r -a lutArray <<< "$LUTs"
for f in "${lutArray[@]}"; do
	[[ -e $f ]] || { echo "Error: LUT not found: $f" >&2; exit 1; }
done
if [[ -n $fitTag && ! $fitTag =~ ^[A-Za-z0-9._-]+$ ]]; then
	echo "Error: --fitTag may contain only letters, digits, \".\", \"_\" and \"-\" (got: $fitTag)" >&2; exit 1
fi
if [[ ! $fillValue =~ ^(NaN|nan|-?[0-9]+(\.[0-9]+)?)$ ]]; then
	echo "Error: --fillValue must be a number or NaN (got: $fillValue)" >&2; exit 1
fi

IFS=',' read -r -a statArray <<< "$targetStats"
segArgs=""
if [[ -n $omitCodes ]]; then segArgs+=",'omitCodes',[${omitCodes//,/ }]"; fi
if [[ -n $omitNames ]]; then
	IFS=',' read -r -a omitArr <<< "$omitNames"
	segArgs+=",'omitNames',$(make_cell "${omitArr[@]}")"
fi
commonArgs="$(make_cell "${tmpFiles[@]}"), $(make_cell "${lutArray[@]}")"
optArgs="'stats',$(make_cell "${statArray[@]}"),'defaultVal',${fillValue},'verbose',${v}${segArgs}"

tmpDir=$(mktemp -d)
matFile="$tmpDir/runPaint.m"
{
	printf "addpath(genpath('%s'));\n" "$(mq "$matDir")"
	printf "nFail = 0;\n"
} > "$matFile"

nqueued=0
if [[ -n $table ]]; then
	# ---- table mode ----
	[[ -e $table ]] || { echo "Error: table not found: $table" >&2; exit 1; }
	table="$(cd "$(dirname "$table")" && pwd)/$(basename "$table")"
	[[ -n $outDir ]] || outDir=$(dirname "$table")
	mkdir -p "$outDir" || exit 1
	{
		printf "try\n"
		printf "    ASE_paintROIStats(%s, '%s', '%s', %s);\n" "$commonArgs" "$(mq "$table")" "$(mq "$outDir/$outName")" "$optArgs"
		printf "catch ME\n    fprintf(2, 'ERROR: %%s\\\\n', ME.message); nFail = nFail + 1;\nend\n"
	} >> "$matFile"
	nqueued=1
else
	# ---- per-subject mode ----
	if [[ -z $singleSubjID ]]; then
		[[ -e $allSubjDir/$dataList ]] || { echo "Error: data list not found: $allSubjDir/$dataList" >&2; exit 1; }
		subjIDs=()
		while IFS= read -r line; do subjIDs+=("$line"); done < <(grep -v '^[[:space:]]*#' "$allSubjDir/$dataList" | awk 'NF {print $1}')
	else
		subjIDs=("$singleSubjID")
	fi
	nSubj=${#subjIDs[@]}
	subjIdx=0
	for subjID in "${subjIDs[@]}"; do
		((subjIdx++))
		fitsDir="$allSubjDir/$subjID/$fitsRelDir"
		if [[ ! -d $fitsDir ]]; then
			echo "No fits folder for $subjID ($fitsDir); skipped."
			continue
		fi
		{
			printf "try\n"
			if ((allRuns)); then
				printf "    d = dir('%s'); d = d([d.isdir] & ~startsWith({d.name}, '.'));\n" "$(mq "$fitsDir")"
				printf "    runs = string(fullfile('%s', {d.name}));\n" "$(mq "$fitsDir")"
				if [[ -n $fitTag ]]; then printf "    runs = runs(endsWith(runs, '_%s'));\n" "$fitTag"; fi
			else
				printf "    [runs, info] = findFitRun('%s', '%s');\n" "$(mq "$fitsDir")" "$fitTag"
				printf "    if runs == \"\", fprintf('%%s: %%s; skipped.\\\\n', '%s', info.reason); runs = strings(0,1); end\n" "$(mq "$subjID")"
			fi
			printf "    for r = runs(:).'\n"
			printf "        tbl = fullfile(r, '%s');\n" "$(mq "$tableName")"
			printf "        if ~isfile(tbl), fprintf('%%s: no %%s in %%s; skipped.\\\\n', '%s', '%s', r); continue; end\n" "$(mq "$subjID")" "$(mq "$tableName")"
			printf "        fprintf('[%d/%d] %s: %%s\\\\n', r);\n" "$subjIdx" "$nSubj" "$(mq "$subjID")"
			printf "        ASE_paintROIStats(%s, tbl, fullfile(r, '%s'), %s);\n" "$commonArgs" "$(mq "$outName")" "$optArgs"
			printf "    end\n"
			printf "catch ME\n    fprintf(2, 'ERROR for %%s: %%s\\\\n', '%s', ME.message); nFail = nFail + 1;\nend\n" "$(mq "$subjID")"
		} >> "$matFile"
		((nqueued++))
	done
fi
printf "if nFail > 0, exit(1); end\n" >> "$matFile"

if ((nqueued > 0)); then
	[[ $v == true ]] && cat "$matFile"
	echo "Opening MATLAB..."
	"$matlabBin" -batch "run('$(mq "$matFile")')"
	exit $?
else
	echo "Nothing queued; MATLAB not started."
fi
