#!/bin/bash
#
# Build an ASE pipeline tree (the layout asePreProcess.sh expects) from the o2co2_ase project,
# using symlinks so no raw data is copied. Run on the cluster.
#
#   source (o2co2_ase)                                   destination (pipeline layout)
#   ------------------                                   -----------------------------
#   datalist_ase, cov_table_ase                       -> <dst>/datalist_ase, cov_table_ase        (links)
#   raw_data/<subj>/dcm/<series>/                     -> <dst>/<subj>/raw_data/dcm/<series>       (links)
#   analyzed/<subj>/dcm.index.dat                     -> <dst>/<subj>/dcm.index.dat               (link)
#   analyzed/<subj>/session_afni_nortc/<file>         -> <dst>/<subj>/session_afni_nortc/<file>   (links)
#
# Everything comes from o2co2_ase; nothing is copied. Real (not linked) directories: <dst>,
# <dst>/<subj>, raw_data, raw_data/dcm, session_afni_nortc, session_ase. The pipeline writes into
# session_ase and, when Dimon converts a series, into raw_data/dcm (GERT_Reco_*, dimon.files.run.*),
# so those must not point back into o2co2_ase. session_afni_nortc is a real folder of per-file
# links so nothing can be created in the original.
#
# Files added to o2co2_ase later (e.g. anat_fs_bstem_al+orig, which no subject has yet) are
# picked up by rerunning this script; asePreProcess.sh then transfers them on its next run.
#
# Safe to rerun: existing links are left alone; links pointing elsewhere are reported.
#
# Usage: linkO2co2Ase.sh --src <o2co2_ase> [--dst <dir>] [--subj <id>]... [--dryRun]

set -u

src=""
dst=""
topFiles="datalist_ase cov_table_ase"
dryRun=0
subjList=()

usage() {
	sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
	cat <<EOF

Options:
  --src <dir>       o2co2_ase folder (required)
  --dst <dir>       output folder (default: <parent of src>/o2co2_ase_proc)
  --subj <id>       only this subject (repeatable; default: every subject in src/analyzed with raw data)
  --dryRun          print what would be done, change nothing
EOF
}

while (($# > 0)); do
	case $1 in
		--src)      src=$2; shift 2 ;;
		--dst)      dst=$2; shift 2 ;;
		--subj)     subjList+=("$2"); shift 2 ;;
		--dryRun)   dryRun=1; shift ;;
		-h|--help)  usage; exit 0 ;;
		*) echo "Error: unknown option $1" >&2; usage >&2; exit 1 ;;
	esac
done

[[ -n $src ]] || { echo "Error: --src is required" >&2; exit 1; }
[[ -d $src ]] || { echo "Error: source not found: $src" >&2; exit 1; }
src="$(cd "$src" && pwd)"
for d in analyzed raw_data; do
	[[ -d $src/$d ]] || { echo "Error: $src/$d not found; is --src the o2co2_ase folder?" >&2; exit 1; }
done
[[ -n $dst ]] || dst="$(dirname "$src")/o2co2_ase_proc"
case "$dst" in /*) ;; *) dst="$(pwd)/$dst" ;; esac
case "$dst/" in "$src"/*) echo "Error: --dst must not be inside --src" >&2; exit 1 ;; esac

nl=$'\n'
problems=""
nLink=0; nDir=0; nKept=0

run() {   # run a command, or print it with --dryRun
	if ((dryRun == 1)); then echo "  [dry] $*"; else "$@"; fi
}

make_dir() {   # make_dir <dir>: real directory (refuse if a link or file is in the way)
	if [[ -L $1 || (-e $1 && ! -d $1) ]]; then
		problems+="  - $1 exists but is not a real directory${nl}"; return 1
	fi
	[[ -d $1 ]] && return 0
	run mkdir -p "$1" && ((++nDir))
}

make_link() {   # make_link <target> <linkPath>
	local target=$1 link=$2
	if [[ -L $link ]]; then
		if [[ $(readlink "$link") == "$target" ]]; then ((++nKept)); return 0; fi
		problems+="  - $link already links to $(readlink "$link") (wanted $target); left alone${nl}"; return 1
	fi
	if [[ -e $link ]]; then
		problems+="  - $link exists and is not a link; left alone${nl}"; return 1
	fi
	run ln -s "$target" "$link" && ((++nLink))
}

if ((${#subjList[@]} == 0)); then
	for d in "$src"/analyzed/*/; do
		s=$(basename "$d")
		[[ -d $src/raw_data/$s/dcm ]] && subjList+=("$s")
	done
fi
((${#subjList[@]} > 0)) || { echo "Error: no subjects found (analyzed/<subj> with raw_data/<subj>/dcm)" >&2; exit 1; }

echo "Source:      $src"
echo "Destination: $dst"
echo "Subjects:    ${#subjList[@]}"
((dryRun == 1)) && echo "(dry run: nothing is changed)"

make_dir "$dst" || { printf '%s' "$problems" >&2; exit 1; }

for f in $topFiles; do
	if [[ -e $src/$f ]]; then make_link "$src/$f" "$dst/$f"
	else problems+="  - $src/$f not found (needed at the top of the pipeline folder)${nl}"; fi
done

for s in "${subjList[@]}"; do
	a="$src/analyzed/$s"
	r="$src/raw_data/$s/dcm"
	o="$dst/$s"
	echo "[$s]"
	if [[ ! -d $a || ! -d $r ]]; then
		problems+="  - $s: missing $a or $r; skipped${nl}"; continue
	fi
	make_dir "$o" || continue
	make_dir "$o/raw_data/dcm" || continue
	make_dir "$o/session_afni_nortc" || continue
	make_dir "$o/session_ase" || continue

	# dicom series: numbered folders only (skips e.g. dcm/deleteme)
	nSer=0
	for d in "$r"/*/; do
		ser=$(basename "$d")
		[[ $ser =~ ^[0-9]+$ ]] || continue
		make_link "$r/$ser" "$o/raw_data/dcm/$ser" && ((++nSer))
	done
	((nSer > 0)) || problems+="  - $s: no numbered series folders in $r${nl}"

	if [[ -e $a/dcm.index.dat ]]; then make_link "$a/dcm.index.dat" "$o/dcm.index.dat"
	else problems+="  - $s: no dcm.index.dat (the comment script will fall back to all dicoms per series)${nl}"; fi

	if [[ -d $a/session_afni_nortc ]]; then
		for f in "$a"/session_afni_nortc/* "$a"/session_afni_nortc/.[!.]*; do
			[[ -e $f || -L $f ]] || continue
			make_link "$f" "$o/session_afni_nortc/$(basename "$f")"
		done
	else
		problems+="  - $s: no session_afni_nortc in $a${nl}"
	fi
done

echo
((dryRun == 1)) && echo -n "(dry run, would be) "
echo "Done: $nDir folders made, $nLink links made, $nKept already in place."
if [[ -n $problems ]]; then
	echo "Problems:"; printf '%s' "$problems"
	exit 2
fi
