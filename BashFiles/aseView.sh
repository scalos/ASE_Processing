#!/bin/bash
# aseView.sh: open ASE_Viewer on a fit result from the command line.
#
# Usage: aseView.sh [options] <manifest.json | directory>
#   A directory is searched (recursively) for the newest *_manifest.json, so you can pass
#   a fit run folder, a subject's fits/ folder, or the subject's session folder.
#   MATLAB runs without the desktop and exits when the viewer window is closed.
#
# Options:
#   --matlabBin <path>   MATLAB executable (default: matlab on PATH, else newest /Applications/MATLAB_R*.app)
#   --matdir <path>      MATLAB code location (default: repo containing this script)
#   -h|--help

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
matdir="$(cd "$scriptDir/.." && pwd)"
matlabBin=matlab
target=""

usage() { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
	case "$1" in
		--matlabBin) [[ $# -ge 2 ]] || { echo "Error: --matlabBin requires an argument" >&2; exit 1; }; matlabBin=$2; shift 2 ;;
		--matdir)    [[ $# -ge 2 ]] || { echo "Error: --matdir requires an argument" >&2; exit 1; };    matdir=$2;    shift 2 ;;
		-h|--help)   usage; exit 0 ;;
		-*)          echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
		*)           if [[ -n $target ]]; then echo "Error: only one manifest/directory allowed" >&2; exit 1; fi
		             target=$1; shift ;;
	esac
done
if [[ -z $target ]]; then usage >&2; exit 1; fi

# ---- resolve the manifest ----
if [[ -d $target ]]; then
	manifest=$(find "$target" -name '*_manifest.json' -type f -print0 2>/dev/null \
		| xargs -0 ls -1t 2>/dev/null | head -n 1)
	if [[ -z $manifest ]]; then echo "Error: no *_manifest.json found under $target" >&2; exit 1; fi
elif [[ -f $target ]]; then
	manifest=$target
else
	echo "Error: not found: $target" >&2; exit 1
fi
manifest="$(cd "$(dirname "$manifest")" && pwd)/$(basename "$manifest")"
if [[ ! -d $matdir ]]; then echo "Error: matdir not found: $matdir" >&2; exit 1; fi
matdir="$(cd "$matdir" && pwd)"

# ---- locate matlab (shell aliases aren't visible to bash) ----
if ! command -v "$matlabBin" > /dev/null 2>&1; then
	found=""
	[[ $matlabBin == matlab ]] && found=$(ls -1d /Applications/MATLAB_R*.app/bin/matlab 2>/dev/null | sort | tail -n 1)
	if [[ -n $found && -x $found ]]; then
		matlabBin=$found
	else
		echo "Error: MATLAB executable not found: $matlabBin (use --matlabBin)" >&2; exit 1
	fi
fi

# mq <str>: escape single quotes for a MATLAB string literal
mq() { local sq="'"; printf '%s' "${1//$sq/$sq$sq}"; }

echo "Opening viewer for: $manifest"
cmd="try"
cmd+=" addpath(genpath('$(mq "$matdir")'));"
cmd+=" V = ASE_Viewer('$(mq "$manifest")'); V.initWindow;"
cmd+=" fprintf('Viewer open; close the window to exit.\n');"
cmd+=" uiwait(V.window);"
cmd+=" catch ME; fprintf(2,'%s\n',getReport(ME)); exit(1); end;"
cmd+=" exit(0);"

exec "$matlabBin" -nodesktop -nosplash -r "$cmd"
