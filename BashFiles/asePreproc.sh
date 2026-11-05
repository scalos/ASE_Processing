#!/bin/bash

# exits the script if the flag has no value after it
need_arg() {
	if [[ $# -lt 2 ]]; then
		echo "Error: $1 requires an argument" >&2
		exit 1
	fi
}

allSubjDir=""
singleSubjDir=""
sessionName="session_ase"

usage() {
    cat << EOF
Usage: $0 [options]
Options:
  --allSubjDir <path>      Path to the directory containing all subjects (required if --singleSubjDir is not provided)
  --singleSubjDir <path>   Path to the directory containing a single subject (required if --allSubjDir is not provided)
  --sessionName <name>     Name of the session (default: $sessionName)
  -h, --help                Show this help message and exit
EOF
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --allSubjDir|--allSubjDir=*)       need_arg "$@"; allSubjDir="$2"; shift 2;;
        --singleSubjDir|--singleSubjDir=*) need_arg "$@"; singleSubjDir="$2"; shift 2;;
        --sessionName|--sessionName=*)     need_arg "$@"; sessionName="$2"; shift 2;;
        -h|--help)           usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1
			;;
	esac
done