#!/bin/bash

# exits the script if the flag has no value after it
need_arg() {
	if [[ $# -lt 2 ]]; then
		echo "Error: $1 requires an argument" >&2
		exit 1
	fi
}

allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
singleSubjDir=""
dataList="datalist_ase"
fitDir="session_ase/fits"
clearAllFits=0
keepLatestN=0
rmLatestN=0

usage() {
    cat << EOF
Usage: $0 [options]
Options:
  --allSubjDir <path>       
  --singleSubjID <path>   Path to the directory of the single subject to process (if provided only this subject will be processed)
  --dataList <name>        Name of the data list file (required). First col is subj id, subsequent are scan names.
  --fitDir <path>          Path to the directory containing fit files relative to subject directories (required)
  --clearAllFits           Clear all fit output for all subjects
  --keepLatestN <num>      Keep the latest N fit runs per subject
  --rmLatestN <num>        Remove the latest N fit runs per subject
  --verbose -v             verbosity
  -h, --help               Show this help message and exit
EOF
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --allSubjDir|--allSubjDir=*)       need_arg "$@"; allSubjDir="$2"; shift 2;;
        --singleSubjID|--singleSubjID=*)   need_arg "$@"; singleSubjID="$2"; shift 2;;
        --dataList|--dataList=*)           need_arg "$@"; dataList="$2"; shift 2;;
        --fitDir|--fitDir=*)               need_arg "$@"; fitDir="$2"; shift 2;;
        --clearAllFits)                    clearAllFits=1; shift ;;
        --keepLatestN|--keepLatestN=*)     need_arg "$@"; keepLatestN="$2"; shift 2;;
        --rmLatestN|--rmLatestN=*)         need_arg "$@"; rmLatestN="$2"; shift 2;;
        --verbose|-v)       v=true; shift ;;
        -h|--help)           usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done

cd "$allSubjDir" || { echo "Error: cannot cd to $allSubjDir" >&2; exit 1; }
if [ -z "$singleSubjID" ]; then
    subjlist=()
	while IFS= read -r line; do
		subjlist+=("$line")
	done < <(grep -v '#' "$allSubjDir/$dataList" | awk '{print $1}')
else
    subjlist=("$singleSubjID")
fi


# ---- maintenance actions: run and exit before touching matlab/transforms ----
if ((clearAllFits == 1)); then
	echo "Clearing all fit output for ${#subjlist[@]} subject(s)..."
	for subj in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$subj/$fitDir"
		if [[ -d $fitsPth ]]; then
			rm -rf "$fitsPth"
			echo "  Removed: $fitsPth"
		fi
	done
	echo "Done."
	exit 0
fi

if ((keepLatestN > 0)); then
	echo "Pruning to the $keepLatestN most recent run(s) per subject..."
	for subj in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$subj/$fitDir"
		[[ -d $fitsPth ]] || continue
		shopt -s nullglob
        runDirs=( "$fitsPth"/* )
        shopt -u nullglob
		nRuns=${#runDirs[@]}
		((nRuns <= keepLatestN)) && continue
		nToDelete=$((nRuns - keepLatestN))
		for ((i = 0; i < nToDelete; i++)); do   # bash arrays are 0-based
			rm -rf "${runDirs[i]}"
			echo "  Removed: ${runDirs[i]}"
		done
	done
	echo "Done."
	exit 0
fi

if ((rmLatestN > 0)); then
	echo "Pruning the most recent $rmLatestN run(s) per subject..."
	for subj in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$subj/$fitDir"
		[[ -d $fitsPth ]] || continue
		shopt -s nullglob
        runDirs=( "$fitsPth"/* )
        shopt -u nullglob
		nRuns=${#runDirs[@]}
		((nRuns < rmLatestN)) && continue
		nToDelete=$((nRuns - rmLatestN))
		for ((i = nRuns - 1; i >= nToDelete; i--)); do
			rm -rf "${runDirs[i]}"
			echo "  Removed: ${runDirs[i]}"
		done
	done
	echo "Done."
	exit 0
fi