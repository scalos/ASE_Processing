#!/bin/bash

# exits the script if the flag has no value after it
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
		out+="'${item}'"
		first=0
	done
	echo "${out}}"
}

# ---- cleanup runs on any exit (normal, error, Ctrl-C) ----
cleanup() {
	if [[ -n $tmpDir && -d $tmpDir ]]; then
		rm -rf "$tmpDir"
		echo "Removed temporary directory: $tmpDir"
	fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM
scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repoDir="$(cd "$scriptDir/.." && pwd)"

allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
singleSubjID=""
tableRelDir="session_ase/fits"
tableName="ase_fit_results_ROI_Stats.csv"
dataList="datalist_ase"
templateDir="/autofs/space/guasha_004/users/crbc/yen/templates"
templateNames=""
targetStats="OEF,R2P,vCBV"
LUTs="/autofs/space/guasha_004/users/crbc/yen/ase/GM_WM_long.txt"
matDir=$repoDir
outName="template_painted_%STAT%"
omitCodes=()
omitNames=""
v=false
matlabBin=matlab

usage() {
    cat << EOF
Usage: $0 [options]
Options:
  --allSubjDir <path>      Path to the directory containing all subjects (default: $allSubjDir)
  --singleSubjID <id>      ID of the single subject to process (if provided only this subject will be processed) (default: $singleSubjID)
  --tableRelDir <path>     Path to stats table relative to subject directories (default: $tableRelDir)
  --tableName <name>       Name of the stats table (default: $tableName)
  --dataList <name>        Name of the data list file. First col is subj id, subsequent are scan names. (default: $dataList)
  --templateDir <path>     Path to the directory containing template files (default: $templateDir)
  --templateNames <names>  Names of the template files (default: $templateNames)
  --LUTs <names>           Names of the Look-Up Table files (default: $LUTs)
  --matDir <path>          Path to directory containing matlab scripts (default: $matDir)
  --targetStats <names>    Names of stats to be painted (default: $targetStats)
  --outName <name>         Name of the output file (default: $outName)
  --omitCodes <codes>      Comma-separated list of codes to omit (default: None)
  --omitNames <names>      Comma-separated list of names to omit (default: None)
  --matlabBin <path>       Path to the MATLAB binary (default: $matlabBin)
  --verbose -v             verbosity
  -h, --help               Show this help message and exit
EOF
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --allSubjDir|--allSubjDir=*)       need_arg "$@"; allSubjDir="$2"; shift 2;;
        --singleSubjID|--singleSubjID=*)   need_arg "$@"; singleSubjID="$2"; shift 2;;
        --tableRelDir|--tableRelDir=*)     need_arg "$@"; tableRelDir="$2"; shift 2;;
        --tableName|--tableName=*)         need_arg "$@"; tableName="$2"; shift 2;;
        --dataList|--dataList=*)           need_arg "$@"; dataList="$2"; shift 2;;
        --templateDir|--templateDir=*)     need_arg "$@"; templateDir="$2"; shift 2;;
        --templateNames|--templateNames=*) need_arg "$@"; templateNames="$2"; shift 2;;
        --LUTs|--LUTs=*)                   need_arg "$@"; LUTs="$2"; shift 2;;
        --targetStats|--targetStats=*)     need_arg "$@"; targetStats="$2"; shift 2;;
        --outName|--outName=*)             need_arg "$@"; outName="$2"; shift 2;;
        --matDir|--matDir=*)               need_arg "$@"; matDir="$2"; shift 2;;
        --omitCodes|--omitCodes=*)         need_arg "$@"; omitCodes="$2"; shift 2;;
        --omitNames|--omitNames=*)         need_arg "$@"; omitNames="$2"; shift 2;;
        --matlabBin|--matlabBin=*)         need_arg "$@"; matlabBin="$2"; shift 2;;
        --verbose|-v)       v=true; shift ;;
        -h|--help)           usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done

if [ -z "$singleSubjID" ]; then
    subjIDs=()
    while IFS= read -r line; do
    subjIDs+=("$line")
    done < <(grep -v '#' "$allSubjDir/$dataList" | awk '{print $1}')
else
    subjIDs=("$singleSubjID")
fi

tmpDir=$(mktemp -d)
matFile="$tmpDir/runPaint.m"
printf "addpath(genpath('$matDir'));\n" > "$matFile"
printf "dbstop if error;" >> "$matFile"

IFS=',' read -r -a templateNameArray <<< "$templateNames"
tmpFiles=()
for templateName in "${templateNameArray[@]}"; do
    tmpFiles+=("$templateDir/$templateName")
done
IFS=',' read -r -a LUTs <<< "$LUTs"
LUTs_cell=$(make_cell "${LUTs[@]}")
tmpFiles_cell=$(make_cell "${tmpFiles[@]}")
#begin loop over subjects:
subjIdx=0
nSubj=${#subjIDs[@]}
nqueued=0
for subjID in "${subjIDs[@]}"; do
    ((subjIdx++))
    subjDir="$allSubjDir/$subjID"
    #discover all stats tables
    if [[ -d "$subjDir" ]]; then
        statTbls=$(find "$subjDir/$tableRelDir" -name "$tableName")
        if [[ -n "$statTbls" ]]; then
            #Found Stats tables. Now loop over them:
            for tbl in "$statTbls"; do
                dirpath=$(dirname "$tbl")
                outPath="$dirpath/$outName"
                printf '%s\n' "runWithProgress(@() ASE_paintROIStats($tmpFiles_cell, $LUTs_cell, \"$tbl\", \"$outPath\",\"verbose\",$v), '$subjID', $subjIdx, $nSubj)" >> "$matFile"
                ((nqueued++))
            done
        else
            if [ "$v" == true ]; then echo "No stats tables found for subject: $subjID"; fi
        fi
    else
        if [ "$v" == true ]; then echo "No directory found for subject: $subjID"; fi
    fi
done

echo "exit;" >> "$matFile"
if (( nqueued > 0 )); then
    cat $matFile
    mFileName=$(basename "$matFile" .m);
	echo "Opening MATLAB..."
	"$matlabBin" -sd "$tmpDir" -nodesktop -nosplash -r "$mFileName"
else
    echo "No subjects queued; skipping matlab"
fi