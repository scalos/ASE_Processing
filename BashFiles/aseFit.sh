#!/bin/bash

# Batch ASE fitting. Uses the outputs of the prep script (ase_e?_proc.nii.gz, ase_hdr.json,
# <seg>_ase.nii.gz). No transforms are performed here.
# Precedence: script defaults < config file ("pipeline" section) < command-line flags.
# Fit options live in the config's "fit" section and are read by runASEFromConfig.m.
# Written for basic bash (no mapfile, no associative arrays, no external dependencies).

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repoDir="$(cd "$scriptDir/.." && pwd)"

# ---- defaults ----
allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
singleSubjID=""
dataList="datalist_ase"
nEchos=2
aseSubjDir=session_ase
aseProcName="ase_e%ECHO%_proc.nii.gz"
hdrFileName=ase_hdr.json
segNames="anat_fs_aparc_aseg_al.nii.gz,anat_fs_wmparc_al.nii.gz,anat_fs_bstem_al.nii.gz"
segSuffix="_ase"
segLUTs="/autofs/space/guasha_004/users/crbc/yen/ase/GM_WM_long.txt"
qualityMaskNames="anat_fs_xvent_mask_al.nii.gz"   # quality masks (before segSuffix); used only if present
fitsSubdir=fits
outPref=ase_fit_results
matlabBin=matlab
matdir=$repoDir
B0=""                    # empty = read B0_T (actual field, from the imaging frequency) from each subject's header
configFile=""
fitTag=""                # optional run label: run folder becomes <timestamp>_<tag>
fitOverrides=""          # accumulated "key":value items from --fit

omitFailLog=0
retryMostRecent=0
retryFromFile=""
logPrefix="fitASE_failures"

echoTok='%ECHO%'
nl=$'\n'
runTimestamp=$(date +%Y%m%d_%H%M%S)

usage() {
cat <<EOF
Usage: $0 [options]

Config:
  --config <file>             JSON config (default: ase_config.json next to this script, if present)
                              "pipeline" keys are read here; "fit" keys are read by MATLAB.
  --fit key=value             Override a "fit" option (repeatable). Value is JSON; bare words are
                              quoted automatically (e.g. --fit physBdClipVal=NaN --fit Hct=0.4
                              --fit physBd_OEF=[0,1]). Nested blocks are merged key by key, e.g.
                              --fit 'quality={"enabled":true,"snr_min":5}' keeps the other quality keys.

Pipeline options (override the config; defaults in brackets):
  --allSubjDir <path>         [$allSubjDir]
  --singleSubjID <id>         Process one subject only
  --dataList <file>           [$dataList]
  --aseSubjDir <dir>          [$aseSubjDir]
  --nEchos <n>                [$nEchos]
  --aseProcName <pattern>     Processed echo file name; must contain $echoTok when nEchos > 1 [$aseProcName]
  --hdrFileName <file>        [$hdrFileName]
  --segNames <a,b,...>        Segmentation names before the suffix is added; empty = no ROI stats [$segNames]
                              Missing ones are left out with a warning; none found = subject skipped
  --segSuffix <str>           [$segSuffix]
  --segLUTs <a,b,...>         LUT file path(s) [$segLUTs]
  --qualityMaskNames <a,...>  Quality masks on the ASE grid (before the suffix); passed to MATLAB as
                              qualityMaskFiles when present (used by fit.quality) [$qualityMaskNames]
  --fitsSubdir <name>         [$fitsSubdir]
  --tag <label>               Label for this run: the run folder becomes <fitsSubdir>/<timestamp>_<label>.
                              Letters, digits, ".", "_" and "-" only. Select it later with
                              aseGroupStats.sh / paintStats.sh --fitTag <label>.
  --outPref <str>             [$outPref]
  --matlabBin <path>          [$matlabBin]
  --matdir <path>             MATLAB code location [$matdir]
  --B0 <value>                Override field strength in T (default: B0_T in each subject's header,
                              i.e. ImagingFrequency / 42.5775; written by asePreProcess.sh)

Logging / retry:
  --omitFailLog               Do not keep a log file
  --retryMostRecent           Skip subjects completed in the most recent log
  --retryFromFile <file>      Skip subjects completed in the given log
  --help|-h
EOF
}

need_arg() {
	if [[ $# -lt 2 ]]; then
		echo "Error: $1 requires an argument" >&2
		exit 1
	fi
}

# mq <str>: escape single quotes for a MATLAB string literal
mq() { local sq="'"; printf '%s' "${1//$sq/$sq$sq}"; }

# make_cell a b c -> {'a','b','c'}
make_cell() {
	local out="" item
	for item in "$@"; do out+="${out:+,}'$(mq "$item")'"; done
	printf '{%s}' "$out"
}

# add_fit_override key=value
add_fit_override() {
	local kv=$1 k v
	if [[ $kv != *=* ]]; then
		echo "Error: --fit expects key=value (got: $kv)" >&2; exit 1
	fi
	k=${kv%%=*}; v=${kv#*=}
	if [[ ! $k =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then
		echo "Error: invalid --fit key: $k" >&2; exit 1
	fi
	if [[ $v =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ || $v == true || $v == false || $v == null \
	      || $v == \[* || $v == \{* || $v == \"* ]]; then
		:
	else
		v="\"$v\""
	fi
	fitOverrides+="${fitOverrides:+,}\"$k\":$v"
}

# cfg_get <pipeline-section-text> <key>: value, or return 1 if missing/null
cfg_get() {
	local line val
	line=$(printf '%s\n' "$1" | grep "^[[:space:]]*\"$2\"[[:space:]]*:" | head -n 1)
	[[ -n $line ]] || return 1
	val=${line#*:}
	val=$(printf '%s' "$val" | sed 's/^[[:space:]]*//; s/[[:space:]]*,[[:space:]]*$//; s/[[:space:]]*$//')
	[[ $val == null ]] && return 1
	if [[ $val == \"*\" ]]; then val=${val#\"}; val=${val%\"}; fi
	printf '%s' "$val"
}

is_num() { [[ $1 =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; }

log_raw() { if [[ -n $failFile && -e $failFile ]]; then printf '%s' "$1" >> "$failFile"; fi; }
record_fail() { printf '%s' "$1"; failLog+=$1; ((nfail++)); }

# ---- locate and load config (before flags, so flags win) ----
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
	if [[ ${args[i]} == --config ]]; then
		if ((i + 1 >= ${#args[@]})); then echo "Error: --config requires an argument" >&2; exit 1; fi
		configFile=${args[i+1]}
	fi
done
if [[ -z $configFile && -e $scriptDir/ase_config.json ]]; then
	configFile="$scriptDir/ase_config.json"
fi
if [[ -n $configFile ]]; then
	if [[ ! -e $configFile ]]; then echo "Error: config file not found: $configFile" >&2; exit 1; fi
	configFile="$(cd "$(dirname "$configFile")" && pwd)/$(basename "$configFile")"
	pipeSec=$(tr -d '\r' < "$configFile" | awk '/"pipeline"[[:space:]]*:/ {f=1; next} f && /^[[:space:]]*}/ {exit} f')
	if [[ -z $pipeSec ]]; then
		echo "Warning: no \"pipeline\" section found in $configFile (script defaults used)"
	fi
	for key in allSubjDir dataList nEchos aseSubjDir aseProcName hdrFileName segNames segSuffix \
	           segLUTs qualityMaskNames fitsSubdir outPref matlabBin matdir B0; do
		if val=$(cfg_get "$pipeSec" "$key"); then printf -v "$key" '%s' "$val"; fi
	done
	echo "Using config: $configFile"
else
	echo "No config file found; using script defaults (MATLAB defaults for fit options)."
fi

# ---- parse flags ----
while [[ $# -gt 0 ]]; do
	case "$1" in
		--config)         need_arg "$@"; shift 2 ;;     # already handled
		--fit)            need_arg "$@"; add_fit_override "$2"; shift 2 ;;
		--allSubjDir)     need_arg "$@"; allSubjDir=$2;     shift 2 ;;
		--singleSubjID)   need_arg "$@"; singleSubjID=$2;   shift 2 ;;
		--dataList)       need_arg "$@"; dataList=$2;       shift 2 ;;
		--aseSubjDir)     need_arg "$@"; aseSubjDir=$2;     shift 2 ;;
		--nEchos)         need_arg "$@"; nEchos=$2;         shift 2 ;;
		--aseProcName)    need_arg "$@"; aseProcName=$2;    shift 2 ;;
		--hdrFileName)    need_arg "$@"; hdrFileName=$2;    shift 2 ;;
		--segNames)       need_arg "$@"; segNames=$2;       shift 2 ;;
		--segSuffix)      need_arg "$@"; segSuffix=$2;      shift 2 ;;
		--segLUTs)        need_arg "$@"; segLUTs=$2;        shift 2 ;;
		--qualityMaskNames) need_arg "$@"; qualityMaskNames=$2; shift 2 ;;
		--fitsSubdir)     need_arg "$@"; fitsSubdir=$2;     shift 2 ;;
		--tag)            need_arg "$@"; fitTag=$2;         shift 2 ;;
		--outPref)        need_arg "$@"; outPref=$2;        shift 2 ;;
		--matlabBin)      need_arg "$@"; matlabBin=$2;      shift 2 ;;
		--matdir)         need_arg "$@"; matdir=$2;         shift 2 ;;
		--B0)             need_arg "$@"; B0=$2;             shift 2 ;;
		--retryFromFile)  need_arg "$@"; retryFromFile=$2;  shift 2 ;;
		--retryMostRecent) retryMostRecent=1; shift ;;
		--omitFailLog)    omitFailLog=1;      shift ;;
		-h|--help)        usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done

if [[ ! $nEchos =~ ^[0-9]+$ ]] || ((nEchos < 1)); then
	echo "Error: nEchos must be a positive integer (got: $nEchos)" >&2; exit 1
fi
if [[ -n $fitTag && ! $fitTag =~ ^[A-Za-z0-9._-]+$ ]]; then
	echo "Error: --tag may contain only letters, digits, \".\", \"_\" and \"-\" (got: $fitTag)" >&2; exit 1
fi
runName=$runTimestamp
[[ -n $fitTag ]] && runName="${runTimestamp}_${fitTag}"
if [[ -n $B0 ]] && ! is_num "$B0"; then
	echo "Error: B0 must be numeric (got: $B0)" >&2; exit 1
fi
if ((nEchos > 1)) && [[ $aseProcName != *"$echoTok"* ]]; then
	echo "Error: aseProcName must contain $echoTok when nEchos > 1 (got: $aseProcName)" >&2; exit 1
fi

overridesJson=""
[[ -n $fitOverrides ]] && overridesJson="{${fitOverrides}}"

segList=(); lutList=(); qmList=()
[[ -n $segNames ]] && IFS=',' read -ra segList <<< "$segNames"
[[ -n $qualityMaskNames ]] && IFS=',' read -ra qmList <<< "$qualityMaskNames"
[[ -n $segLUTs ]]  && IFS=',' read -ra lutList <<< "$segLUTs"
nSeg=${#segList[@]}
if ((nSeg > 0 && ${#lutList[@]} == 0)); then
	echo "Error: segNames given but segLUTs is empty" >&2; exit 1
fi

# ---- cleanup (any exit) ----
tmpDir=""; queuedDirs=(); queuedSubjs=(); queuedSubjBases=(); failFile=""
cleanup() {
	local d
	if [[ -n $tmpDir && -d $tmpDir ]]; then rm -rf "$tmpDir"; fi
	for d in "${queuedDirs[@]}"; do rmdir "$d" 2>/dev/null && echo "Removed empty run folder: $d"; done
	for d in "${queuedSubjBases[@]}"; do rmdir "$d/$fitsSubdir" 2>/dev/null; done
	if [[ -n $failFile && -e $failFile && $omitFailLog -eq 1 ]]; then rm -f "$failFile"; fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

# ---- subject list ----
# resolve to absolute: the loop cd's into each subject dir, which would break a relative path
if [[ ! -d $allSubjDir ]]; then echo "Error: allSubjDir not found: $allSubjDir" >&2; exit 1; fi
allSubjDir="$(cd "$allSubjDir" && pwd)"
if [[ -n $retryFromFile && -e $retryFromFile ]]; then
	retryFromFile="$(cd "$(dirname "$retryFromFile")" && pwd)/$(basename "$retryFromFile")"
fi
for ((li = 0; li < ${#lutList[@]}; li++)); do
	lut=${lutList[li]}
	if [[ $lut != /* && -d $(dirname "$lut") ]]; then
		lutList[li]="$(cd "$(dirname "$lut")" && pwd)/$(basename "$lut")"
	fi
done
if [[ -n $matdir && $matdir != /* && -d $matdir ]]; then matdir="$(cd "$matdir" && pwd)"; fi
# locate matlab: shell aliases (e.g. zsh) are not visible to bash, so fall back to a standard install
if ! command -v "$matlabBin" > /dev/null 2>&1; then
	found=""
	if [[ $matlabBin == matlab ]]; then
		found=$(ls -1d /Applications/MATLAB_R*.app/bin/matlab 2>/dev/null | sort | tail -n 1)
	fi
	if [[ -n $found && -x $found ]]; then
		echo "Note: '$matlabBin' not on PATH; using $found"
		matlabBin=$found
	else
		echo "Error: MATLAB executable not found: $matlabBin (set --matlabBin or matlabBin in the config)" >&2
		exit 1
	fi
fi
# accept a subject folder path as well as a bare ID
if [[ -n $singleSubjID ]]; then
	singleSubjID=${singleSubjID%/}
	singleSubjID=${singleSubjID##*/}
fi
cd "$allSubjDir" || { echo "Error: cannot cd to $allSubjDir" >&2; exit 1; }
if [ -z "$singleSubjID" ]; then
	subjlist=()
	while IFS= read -r line; do
		subjlist+=("$line")
	done < <(grep -v '#' "$allSubjDir/$dataList" | awk 'NF {print $1}')
else
	subjlist=("$singleSubjID")
fi
nsubj=${#subjlist[@]}

# ---- retry handling ----
retryFile=""
if ((retryMostRecent == 1)); then
	retryFile=$(ls -1t "$allSubjDir"/${logPrefix}_*.txt 2>/dev/null | head -n 1)
fi
if [[ -n $retryFromFile ]]; then
	if [[ -n $retryFile && $retryFile != "$retryFromFile" ]]; then
		echo "Warning: --retryMostRecent and --retryFromFile refer to different files; using --retryFromFile ($retryFromFile)"
	fi
	retryFile=$retryFromFile
fi
if [[ -n $retryFile && ! -e $retryFile ]]; then
	echo "Error: retry requested but specified file not found: $retryFile" >&2
	exit 1
fi
completedSubjs=()
if [[ -n $retryFile ]]; then
	while IFS= read -r line; do
		completedSubjs+=("$line")
	done < <(awk '/^\+/ {print substr($1, 2)}' "$retryFile")
	printf 'Skipping completed subjects from file: %s\nCompleted subjects:\n' "$retryFile"
	for subj in "${completedSubjs[@]}"; do printf '  - %s\n' "$subj"; done
fi

failFile="$allSubjDir/${logPrefix}_${runTimestamp}.txt"
echo "Run started: $(date)" > "$failFile"
echo "Completed subjects:" >> "$failFile"

tmpDir=$(mktemp -d)
mfile="$tmpDir/run_ase_fit_batch.m"
echo "addpath(genpath('$(mq "$matdir")'));" > "$mfile"

failLog=""; nfail=0; exitCode=0

echo "=================================================="
echo "Preparing ASE fitting for $nsubj subject(s)..."
subj_i=0
for x in "${subjlist[@]}"; do
	((subj_i++))
	if ((${#completedSubjs[@]} > 0)); then
		if printf '%s\n' "${completedSubjs[@]}" | grep -qxF "$x"; then
			echo "Skipping $x (already completed)"
			log_raw "+$x${nl}"
			continue
		fi
	fi

	subjBase="$allSubjDir/$x/$aseSubjDir"
	okay=1; missingReasons=""; warnMsg=""

	if [[ ! -d $subjBase ]]; then
		record_fail "Skipping $x:${nl}  - subject directory not found: $subjBase${nl}"
		continue
	fi
	cd "$subjBase" || continue
	prepHint="(run the prep script first)"

	# ---- header ----
	hdrPth="$subjBase/$hdrFileName"
	if [[ ! -e $hdrPth ]]; then
		okay=0; missingReasons+="  - header file not found: $hdrPth $prepHint${nl}"
	fi

	# ---- ase files: one processed file per echo -> nEchos x 1 cell (rows = echos) ----
	aseList=""
	for ((ei = 1; ei <= nEchos; ei++)); do
		fn=${aseProcName//"$echoTok"/$ei}
		if [[ ! -e $subjBase/$fn ]]; then
			okay=0; missingReasons+="  - echo $ei file not found: $subjBase/$fn $prepHint${nl}"
			continue
		fi
		aseList+="${aseList:+;}'$(mq "$subjBase/$fn")'"
	done
	asePths="{${aseList}}"

	# ---- B0 ----
	subjB0=$B0
	if [[ -z $subjB0 && -e $hdrPth ]]; then
		subjB0=$(tr -d '\r' < "$hdrPth" | sed -n 's/^[[:space:]]*"B0_T"[[:space:]]*:[[:space:]]*\([-0-9.eE+]*\).*/\1/p' | head -n 1)
		if ! is_num "$subjB0"; then
			okay=0; missingReasons+="  - no B0_T in $hdrPth (header predates it): rerun asePreProcess.sh --regenHdr, or pass --B0${nl}"
		fi
	fi

	# ---- segmentation (already in ase space) and LUTs ----
	segArgs=""
	if ((nSeg > 0)); then
		# a missing segmentation is left out with a warning (its ROIs are absent from the table);
		# the subject is skipped only when none of them is found
		segPthList=()
		segMissing=""
		for segName in "${segList[@]}"; do
			base=${segName%.gz}; base=${base%.nii}
			segOut="${base}${segSuffix}.nii.gz"
			if [[ -e "$subjBase/$segOut" ]]; then
				segPthList+=("$subjBase/$segOut")
			else
				segMissing+="  - transformed segmentation not found: $subjBase/$segOut $prepHint${nl}"
			fi
		done
		if ((${#segPthList[@]} == 0)); then
			okay=0; missingReasons+="$segMissing"
		elif [[ -n $segMissing ]]; then
			warnMsg+="${segMissing//  - /  - warning: left out of the ROI table, }"
		fi
		for lut in "${lutList[@]}"; do
			if [[ ! -e $lut ]]; then
				okay=0; missingReasons+="  - LUT file not found: $lut${nl}"
			fi
		done
		segArgs=",'segFiles',$(make_cell "${segPthList[@]}"),'segLUTs',$(make_cell "${lutList[@]}")"
	fi

	# ---- quality masks (already on the ASE grid); optional: missing ones are only warned about ----
	qmArgs=""
	if ((${#qmList[@]} > 0)); then
		qmPthList=()
		for qmName in "${qmList[@]}"; do
			base=${qmName%.gz}; base=${base%.nii}
			qmOut="${base}${segSuffix}.nii.gz"
			if [[ -e "$subjBase/$qmOut" ]]; then
				qmPthList+=("$subjBase/$qmOut")
			else
				warnMsg+="  - quality mask not found (quality criteria using it are skipped): $subjBase/$qmOut${nl}"
			fi
		done
		((${#qmPthList[@]} > 0)) && qmArgs=",'qualityMaskFiles',$(make_cell "${qmPthList[@]}")"
	fi

	if ((okay == 0)); then
		record_fail "Skipping $x:${nl}${missingReasons}"
		continue
	fi
	if [[ -n $warnMsg ]]; then
		printf 'Warnings for %s:\n%s' "$x" "$warnMsg"
		log_raw "Warnings for $x:${nl}${warnMsg}"
	fi

	# ---- run folder and matlab call ----
	runOutDir="$subjBase/$fitsSubdir/$runName"
	mkdir -p "$runOutDir"
	queuedDirs+=("$runOutDir"); queuedSubjs+=("$x"); queuedSubjBases+=("$subjBase")
	outPrefFull="$runOutDir/$outPref"

	cfgArg=""; ovArg=""
	[[ -n $configFile ]]    && cfgArg=",'cfgFile','$(mq "$configFile")'"
	[[ -n $overridesJson ]] && ovArg=",'overrides','$(mq "$overridesJson")'"

	{
		printf '%s\n' "try"
		printf '%s\n' "    runWithProgress(@() runASEFromConfig(${subjB0},'$(mq "$hdrPth")',${asePths},'$(mq "$runOutDir")'${cfgArg}${ovArg}${segArgs}${qmArgs},'outPref','$(mq "$outPrefFull")'), '$(mq "$x")', ${subj_i}, ${nsubj});"
		printf '%s\n' "    fID = fopen('$(mq "$failFile")', 'at');"
		printf '%s\n' "    if fID ~= -1; fprintf(fID, '+%s\n', '$(mq "$x")'); fclose(fID); end"
		printf '%s\n' "catch ME"
		printf '%s\n' "    fprintf(2, 'ERROR for %s: %s\n', '$(mq "$x")', ME.message);"
		printf '%s\n' "    fID = fopen('$(mq "$failFile")', 'at');"
		printf '%s\n' "    if fID ~= -1; fprintf(fID, 'Failed %s (MATLAB error):\n  - %s\n', '$(mq "$x")', strrep(ME.message, newline, ' ')); fclose(fID); end"
		printf '%s\n' "end"
	} >> "$mfile"
done
cd "$allSubjDir" || exit 1

# ---- run matlab once; failures are isolated by try/catch per subject ----
nqueued=${#queuedSubjs[@]}
matStatus=0
if ((nqueued > 0)); then
	echo "Opening MATLAB for $nqueued subject(s)..."
	"$matlabBin" -batch "run('$(mq "$mfile")')"
	matStatus=$?
else
	echo "No subjects queued; skipping matlab."
fi

# ---- tally results from the log ----
nMatFail=0; matFailList=""
for x in "${queuedSubjs[@]}"; do
	if ! grep -qxF "+$x" "$failFile"; then
		((nMatFail++)); matFailList+="  - $x${nl}"
	fi
done
if ((matStatus != 0 && nMatFail > 0)); then
	if ((matStatus == 126 || matStatus == 127)); then
		failLog+="MATLAB could not be launched (status $matStatus): $matlabBin${nl}"
		echo "Error: MATLAB could not be launched (status $matStatus): $matlabBin" >&2
	else
		failLog+="MATLAB exited with status $matStatus${nl}"
	fi
fi

echo ""
echo "=================================================="
echo "Run complete: $nqueued subject(s) queued, $((nqueued - nMatFail)) fitted, $nMatFail failed in MATLAB, $nfail skipped before MATLAB"
if ((nfail == 0 && nMatFail == 0)); then
	echo "No failures."
	rm -f "$failFile"
else
	exitCode=1
	if ((nfail > 0)); then echo "Skipped before MATLAB:"; printf '%s' "$failLog"; fi
	if ((nMatFail > 0)); then echo "Failed in MATLAB (details in log):"; printf '%s' "$matFailList"; fi
	if ((omitFailLog == 0)); then
		{
			echo ""
			echo "$nfail skipped, $nMatFail failed in MATLAB"
			echo ""
			printf '%s' "$failLog"
		} >> "$failFile"
		echo "Failure details written to: $failFile"
	else
		echo "Failure details omitted (omitFailLog=1)"
	fi
fi
echo "=================================================="
exit $exitCode