#!/bin/bash
# Siemens_ase_XA30.sh: builds the per-frame image-comment header for an ASE series.
# Bash port of the original tcsh script (same DICOM tags, same output format).
#
# Output: one line per DICOM, sorted (stably) by acquisition number:
#     <dicom file> <acq number> <acq time in ms> <image comment>
# Image comment (0020,4000) has spaces and '*' removed and ',' replaced by '_'.
#
# Written for Siemens XA30 (Bay 4, after the upgrade). If a scanner upgrade changes the
# DICOM layout or tags, edit this script (or add a sibling) and point commentScript at it.
#
# Files: every DICOM in the series folder (one folder per series, as dcmunpack writes them).
# The folder must hold a single series: all files must share the DICOM series number (0020,0011)
# and study (0020,000D), else the script stops. The dicom index is no longer used to pick files:
# it describes only one exam, and a series from a second exam unpacked into the same dcm folder
# (renumbered, e.g. folder 116 for series 16, bch045_250823) is missing from it.

usage() {
cat <<EOF
Usage: $0 direcdcm direcdcmindx seqno scanner [outFile]

  direcdcm      dicom path (folder holding this series' dicoms)
  direcdcmindx  accepted for compatibility, not used (files are taken from direcdcm)
  seqno         sequence (series) number, e.g. 22 or 022 (the folder's number; if the DICOMs carry
                another series number, e.g. a renumbered second exam, a note is printed)
  scanner       1:Trio 2:Skyra 3:Prisma (accepted for compatibility; not used yet)
  outFile       output file (default: ./dcm_hdr_<seqno, 3 digits> in the current directory)

Notes:
  (1) only tested on dicom data acquired at Bay 4 after the upgrade to XA30
  (2) assumes one brain volume (per time point) per dicom file
EOF
}

if [[ $# -lt 4 || $# -gt 5 ]]; then usage; exit 1; fi

direcdcm=$1
direcdcmindx=$2
seqno=$3
scanner=$4

if [[ ! $seqno =~ ^[0-9]+$ ]]; then
	echo "Error: seqno must be a number (got: $seqno)" >&2; exit 1
fi
if [[ ! $scanner =~ ^[0-9]+$ ]]; then
	echo "Error: scanner must be a number (got: $scanner)" >&2; exit 1
fi
if ! command -v dicom_hinfo > /dev/null 2>&1; then
	echo "Error: dicom_hinfo (AFNI) not found on PATH" >&2; exit 1
fi
if [[ ! -d $direcdcm ]]; then
	echo "Error: dicom directory not found: $direcdcm" >&2; exit 1
fi

seqno_tmp=$(printf '%03d' $((10#$seqno)))
outFile=${5:-"$PWD/dcm_hdr_${seqno_tmp}"}

rowsFile=$(mktemp)
sortedFile=$(mktemp)
trap 'rm -f "$rowsFile" "$sortedFile"' EXIT

# ---- list of dicom file names (relative to direcdcm): every dicom in the series folder ----
nameList=""
shopt -s nullglob
for p in "$direcdcm"/*MR.* "$direcdcm"/Trio* "$direcdcm"/*IM*; do
	[[ -f $p ]] && nameList+="${p##*/}"$'\n'
done
shopt -u nullglob
nameList=$(printf '%s' "$nameList" | awk 'NF && !seen[$0]++')   # drop duplicates from overlapping patterns
if [[ -z $nameList ]]; then
	echo "Error: no dicom files (MR.*, Trio*, *IM*) found in $direcdcm" >&2; exit 1
fi

# ---- one row per dicom ----
# (-no_name: dicom_hinfo otherwise prefixes the file path, which breaks field parsing if it has spaces)
# numeric tags in one call: series number (0020,0012 acquisition number, 0018,9074 acquisition time,
# 0020,000D study UID)
n=0; firstSer=""; firstStudy=""; firstFile=""
while IFS= read -r x; do
	[[ -n $x ]] || continue
	path="$direcdcm/$x"
	if [[ ! -r $path ]]; then
		echo "Error: cannot read dicom file: $path" >&2; exit 1
	fi

	read -r ser acq tt study < <(dicom_hinfo -no_name -tag 0020,0011 0020,0012 0018,9074 0020,000d "$path")

	# one series per folder
	if [[ -z $firstFile ]]; then
		firstSer=$ser; firstStudy=$study; firstFile=$x
	elif [[ $ser != "$firstSer" || $study != "$firstStudy" ]]; then
		echo "Error: $direcdcm holds more than one series: $firstFile is series $firstSer (study $firstStudy)," >&2
		echo "       $x is series $ser (study $study)" >&2
		exit 1
	fi

	# acquisition number (0020,0012)
	if [[ ! $acq =~ ^[0-9]+$ ]]; then
		echo "Error: could not read acquisition number (0020,0012) from $path" >&2; exit 1
	fi
	acq=$((10#$acq))

	# acquisition time (0018,9074: YYYYMMDDHHMMSS.FFFFFF) -> ms since midnight
	tt=$(printf '%s' "$tt" | cut -c9-30)
	if [[ ! $tt =~ ^[0-9]{6}[.][0-9]{4} ]]; then
		echo "Error: could not parse frame acquisition time (0018,9074) from $path (got '$tt')" >&2; exit 1
	fi
	hr=${tt:0:2}; mn=${tt:2:2}; sc=${tt:4:2}; tick=${tt:7:4}
	tms=$(( 10#$hr*3600000 + 10#$mn*60000 + 10#$sc*1000 + 10#$tick/10 ))

	# image comment (0020,4000)
	cmt=$(dicom_hinfo -no_name -tag 0020,4000 -full_entry "$path" | sed -e 's/ //g' -e 's/,/_/g' -e 's/\*//g')

	((n++))
	printf '%s %s %s %s %s%s\n' "$acq" "$n" "$x" "$acq" "$tms" "${cmt:+ $cmt}" >> "$rowsFile"
done <<< "$nameList"

if [[ ! $firstSer =~ ^[0-9]+$ ]]; then
	echo "Error: could not read the series number (0020,0011) from $direcdcm/$firstFile" >&2; exit 1
fi
if ((10#$firstSer != 10#$seqno)); then
	echo "Note: folder $seqno holds DICOM series $firstSer (renumbered when unpacked, e.g. a second exam); $n files used"
fi

# ---- stable sort by acquisition number (ties keep original order), drop sort keys ----
sort -n -k1,1 -k2,2 "$rowsFile" | cut -d' ' -f3- > "$sortedFile"

mkdir -p "$(dirname "$outFile")"
if ! cp "$sortedFile" "$outFile.tmp" || ! mv "$outFile.tmp" "$outFile"; then
	echo "Error: could not write $outFile" >&2; exit 1
fi
echo "Wrote $n rows to $outFile"