#!/bin/bash

# asePreProcess.sh: prepares ASE subject directories for fitting, starting from raw DICOMs.
# Per subject, per echo series:
#   1. convert dicoms -> 4D nifti (all frames, acquisition order) + dcm2niix json sidecar; de-oblique
#   2. build the dTE / z-moment / image-comment header (separate script, see commentScript)
#   3. build the combined json header (all echos)
#   4. motion correct (3dvolreg) -> <aseProcName>
#   5. transfer anatomical/segmentation volumes from transferSrcDir (+orig -> nifti)   [if missing]
#   6. register ase_base to the anatomical (3dAllineate) and make the inverse matrix   [if missing]
#   7. transform segmentations into ASE space (nearest neighbor)
# Existing outputs are skipped unless regenerated. Regenerating an upstream stage forces its
# downstream stages (and overwrites). Settings: script defaults < config "pipeline" section < flags.
# Stages 5-7 are skipped when segNames is empty.

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------- defaults (all can be set in the config "pipeline" section) ----------------
allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
singleSubjID=""
dataList="datalist_ase"
echo1Col=28
echoSpacing=2
nEchos=2
aseSubjDir=session_ase
rawDataDir=raw_data/dcm
dcmIndexFile=dcm.index.dat          # subject-relative (sits next to session_ase)
scannerCode=3
commentScript=Siemens_ase_XA30.sh
dcm2niixMinFiles=500
dimonOpts=""
rawNii="ase_%SERIES%.nii"
pattern_dicomHdr="dcm_hdr_%SERIES%.json"
pattern_dteHdr="dcm_hdr_%SERIES%_dtes.txt"
aseProcName="ase_e%ECHO%_proc.nii.gz"
motionParamName="ase_e%ECHO%_mc.1D"
hdrFileName=ase_hdr.json
volregOpts="-Fourier -twopass -twodup"
baseName=ase_base.nii.gz
baseFrame=0
baseAlName=ase_base_al.nii.gz
regBaseName=anat_fs_brain_al.nii.gz
allineateOpts="-twopass -cost lpa+ -autoweight -fineblur 1 -cmass"
tformName=ase_base_al_mat.aff12.1D
tformApplyName=ase_base_al_inv_mat.aff12.1D
transferSrcDir=session_afni_nortc
transferNames="anat_fs_brain_al.nii.gz,anat_fs_brain_al_mask.nii.gz,anat_fs_xvent_mask_al.nii.gz,anat_fs_bstem_al.nii.gz"
segNames="anat_fs_aparc_aseg_al.nii.gz,anat_fs_wmparc_al.nii.gz"
segSuffix="_ase"

pipelineKeys="allSubjDir dataList echo1Col echoSpacing nEchos aseSubjDir rawDataDir dcmIndexFile scannerCode commentScript dcm2niixMinFiles dimonOpts rawNii pattern_dicomHdr pattern_dteHdr aseProcName motionParamName hdrFileName volregOpts baseName baseFrame baseAlName regBaseName allineateOpts tformName tformApplyName transferSrcDir transferNames segNames segSuffix"

# regeneration / logging / retry
regenNifti=0; regenComments=0; regenHdr=0; regenMotion=0
regenTransfer=0; regenRegister=0; regenTransform=0
omitFailLog=0; retryMostRecent=0; retryFromFile=""
logPrefix="prepASE_failures"
configFile=""

seriesTok='%SERIES%'
echoTok='%ECHO%'
nl=$'\n'
runTimestamp=$(date +%Y%m%d_%H%M%S)

usage() {
cat <<EOF
Usage: $0 [options]

  --config <file>        JSON config; its "pipeline" section sets defaults
                         (default: ase_config.json next to this script, if present)
  --set key=value        Override any pipeline key (repeatable), e.g. --set segSuffix=_ase
                         Keys: $pipelineKeys
  --allSubjDir <path>    Top-level directory containing subject folders
  --singleSubjID <id>    Process one subject only (default: all in dataList)
  --dataList <file>      Subject list in allSubjDir
  --aseSubjDir <dir>     Subject-relative ase directory
  --nEchos <n>           Number of echos

Regeneration (overwrites; downstream stages are forced too):
  --regenNifti           dicom -> nifti + sidecar (forces header and motion correction)
  --regenComments        dTE / image-comment header (forces combined header)
  --regenHdr             combined json header only
  --regenMotion          motion correction only
  --regenTransfer        re-copy anatomical/segmentation files (forces segmentation transforms)
  --regenRegister        re-register ase_base to anatomical (forces segmentation transforms)
  --regenTransform       segmentation transforms only

Logging / retry:
  --omitFailLog | --retryMostRecent | --retryFromFile <file> | --help
EOF
}

need_arg() {
	if [[ $# -lt 2 ]]; then echo "Error: $1 requires an argument" >&2; exit 1; fi
}

set_key() {
	local kv=$1 k v
	if [[ $kv != *=* ]]; then echo "Error: --set expects key=value (got: $kv)" >&2; exit 1; fi
	k=${kv%%=*}; v=${kv#*=}
	case " $pipelineKeys " in
		*" $k "*) ;;
		*) echo "Error: unknown pipeline key: $k" >&2; exit 1 ;;
	esac
	printf -v "$k" '%s' "$v"
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
	pipeSec=$(tr -d '\r' < "$configFile" | awk '/"pipeline"[[:space:]]*:/ {f=1; next} f && /^[[:space:]]*}/ {exit} f')
	if [[ -z $pipeSec ]]; then
		echo "Warning: no \"pipeline\" section found in $configFile (script defaults used)"
	fi
	for key in $pipelineKeys; do
		if val=$(cfg_get "$pipeSec" "$key"); then printf -v "$key" '%s' "$val"; fi
	done
	echo "Using config: $configFile"
fi

# ---- parse flags ----
while [[ $# -gt 0 ]]; do
	case "$1" in
		--config)          need_arg "$@"; shift 2 ;;
		--set)             need_arg "$@"; set_key "$2"; shift 2 ;;
		--allSubjDir)      need_arg "$@"; allSubjDir=$2;     shift 2 ;;
		--singleSubjID)    need_arg "$@"; singleSubjID=$2;   shift 2 ;;
		--dataList)        need_arg "$@"; dataList=$2;       shift 2 ;;
		--aseSubjDir)      need_arg "$@"; aseSubjDir=$2;     shift 2 ;;
		--nEchos)          need_arg "$@"; nEchos=$2;         shift 2 ;;
		--retryFromFile)   need_arg "$@"; retryFromFile=$2;  shift 2 ;;
		--regenNifti)      regenNifti=1;      shift ;;
		--regenComments)   regenComments=1;   shift ;;
		--regenHdr)        regenHdr=1;        shift ;;
		--regenMotion)     regenMotion=1;     shift ;;
		--regenTransfer)   regenTransfer=1;   shift ;;
		--regenRegister)   regenRegister=1;   shift ;;
		--regenTransform)  regenTransform=1;  shift ;;
		--omitFailLog)     omitFailLog=1;     shift ;;
		--retryMostRecent) retryMostRecent=1; shift ;;
		-h|--help)         usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use --help for details)" >&2
			exit 1 ;;
	esac
done

# ---- resolve paths to absolute (script cd's into subject dirs later) ----
if [[ ! -d $allSubjDir ]]; then echo "Error: allSubjDir not found: $allSubjDir" >&2; exit 1; fi
allSubjDir="$(cd "$allSubjDir" && pwd)"
if [[ -n $retryFromFile && -e $retryFromFile ]]; then
	retryFromFile="$(cd "$(dirname "$retryFromFile")" && pwd)/$(basename "$retryFromFile")"
fi

# ---- validate settings ----
is_int() { [[ $1 =~ ^[0-9]+$ ]]; }
for k in nEchos echo1Col echoSpacing scannerCode dcm2niixMinFiles baseFrame; do
	if ! is_int "${!k}"; then echo "Error: $k must be a non-negative integer (got: ${!k})" >&2; exit 1; fi
done
if ((nEchos < 1)); then echo "Error: nEchos must be >= 1" >&2; exit 1; fi
for k in rawNii pattern_dicomHdr pattern_dteHdr; do
	if [[ ${!k} != *"$seriesTok"* ]]; then echo "Error: $k must contain $seriesTok (got: ${!k})" >&2; exit 1; fi
done
for k in aseProcName motionParamName; do
	if ((nEchos > 1)) && [[ ${!k} != *"$echoTok"* ]]; then
		echo "Error: $k must contain $echoTok when nEchos > 1 (got: ${!k})" >&2; exit 1
	fi
done

# ---------------- helpers ----------------
is_num() { [[ $1 =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; }
subst_series() { local s=${1//"$seriesTok"/$2}; printf '%s' "$s"; }
subst_echo()   { local s=${1//"$echoTok"/$2};   printf '%s' "$s"; }
strip_nii()    { local b=${1%.gz}; b=${b%.nii}; printf '%s' "$b"; }

# json_num <file> <key>: value of a top-level "key": number line (dcm2niix layout)
json_num() {
	tr -d '\r' < "$1" | sed -n "s/^[[:space:]]*\"$2\"[[:space:]]*:[[:space:]]*\([-+0-9.eE]*\)[[:space:]]*,\{0,1\}[[:space:]]*$/\1/p" | head -n 1
}
dte_field() { tr -d '\r' < "$1" | sed -n "$2"; }
count_lines() { printf '%s\n' "$1" | grep -c .; }

dteExpr='s/.*_DTE=\([-0-9.]*\)us_.*/\1/p'
zmExpr='s/.*ZMoment=\([^T]*\)Tms\/m.*/\1/p'
maxExpr='s/.*_MaxDTE=\([-0-9.]*\)us.*/\1/p'
stepExpr='s/.*_dTE=\([-0-9.]*\)us.*/\1/p'

log_raw() { if ((omitFailLog == 0)); then printf '%s' "$1" >> "$failFile"; fi; }
record_fail() { printf '%s' "$1"; failLog+=$1; log_raw "$1"; ((nfail++)); }

# find the comment script
commentScriptPath=""
if [[ $commentScript == /* ]]; then
	commentScriptPath=$commentScript
elif [[ -e $scriptDir/$commentScript ]]; then
	commentScriptPath=$scriptDir/$commentScript
elif command -v "$commentScript" > /dev/null 2>&1; then
	commentScriptPath=$(command -v "$commentScript")
fi

# segmentation / transfer lists
segArr=(); transferArr=()
[[ -n $segNames ]] && IFS=',' read -ra segArr <<< "$segNames"
[[ -n $transferNames ]] && IFS=',' read -ra transferArr <<< "$transferNames"
nSeg=${#segArr[@]}
xferList=()
if ((nSeg > 0)); then
	for name in "${transferArr[@]}" "${segArr[@]}"; do
		dup=0
		for e in "${xferList[@]}"; do [[ $e == "$name" ]] && dup=1; done
		((dup)) || xferList+=("$name")
	done
fi
nXfer=${#xferList[@]}

xfer_index() { local j; for ((j = 0; j < nXfer; j++)); do [[ ${xferList[j]} == "$1" ]] && { printf '%d' "$j"; return 0; }; done; return 1; }

# avail <name>: file exists in subjBase or is scheduled to be transferred
avail() {
	local j
	[[ -e $subjBase/$1 ]] && return 0
	for ((j = 0; j < nXfer; j++)); do
		[[ ${xferList[j]} == "$1" && ${xferDo[j]} == 1 ]] && return 0
	done
	return 1
}

# tools
needTools="dcm2niix 3dinfo 3drefit 3dvolreg 3dcalc"
((nSeg > 0)) && needTools+=" 3dAllineate cat_matvec 3dAFNItoNIFTI"
for t in $needTools; do
	if ! command -v "$t" > /dev/null 2>&1; then echo "Error: $t not found on PATH" >&2; exit 1; fi
done

# ---------------------------------------------------------------------------
# planning: which stages run for a subject, and are their inputs available?
# ---------------------------------------------------------------------------
check_subject() {
	local ei series nii side dte proc mc nraw j k name base idx
	missingReasons=""; ok=1
	echoNums=(); doNiiE=(); doSideE=(); doCommE=(); doMotionE=(); xferDo=(); segDo=()
	anyNii=0; anySide=0; anyComm=0; anyMotion=0
	doRegister=0; doInverse=0; doHdr=0

	if [[ ! -d $subjDir ]]; then
		missingReasons+="  - subject directory not found: $subjDir${nl}"
		return 1
	fi

	local datarow
	datarow=$(awk -v s="$subjID" '$1==s {print; exit}' "$allSubjDir/$dataList")
	datatmp=( $datarow )
	e1Num=${datatmp[echo1Col-1]//ase_/}
	if [[ ! $e1Num =~ ^[0-9]+$ ]]; then
		missingReasons+="  - could not read echo-1 series number from column $echo1Col of $dataList (got '${datatmp[echo1Col-1]}')${nl}"
		return 1
	fi
	local width=${#e1Num}

	for ((ei = 0; ei < nEchos; ei++)); do
		series=$(printf "%0${width}d" $((10#$e1Num + ei * echoSpacing)))
		echoNums+=("$series")
		nii=$(subst_series "$rawNii" "$series")
		side=$(subst_series "$pattern_dicomHdr" "$series")
		dte=$(subst_series "$pattern_dteHdr" "$series")
		proc=$(subst_echo "$aseProcName" "$((ei + 1))")
		mc=$(subst_echo "$motionParamName" "$((ei + 1))")

		doNiiE[ei]=0; doSideE[ei]=0; doCommE[ei]=0; doMotionE[ei]=0
		if ((regenNifti == 1)) || [[ ! -e $subjBase/$nii ]];  then doNiiE[ei]=1;  anyNii=1;  fi
		if ((regenNifti == 1)) || [[ ! -e $subjBase/$side ]]; then doSideE[ei]=1; anySide=1; fi
		if ((regenComments == 1)) || [[ ! -e $subjBase/$dte ]]; then doCommE[ei]=1; anyComm=1; fi
		if ((regenMotion == 1 || doNiiE[ei] == 1)) || [[ ! -e $subjBase/$proc || ! -e $subjBase/$mc ]]; then
			doMotionE[ei]=1; anyMotion=1
		fi

		if ((doNiiE[ei] == 1 || doSideE[ei] == 1 || doCommE[ei] == 1)); then
			if [[ ! -d $rawBase/$series ]]; then
				ok=0; missingReasons+="  - raw dicom directory not found: $rawBase/$series${nl}"
			elif ((doNiiE[ei] == 1)); then
				nraw=$(ls -1 "$rawBase/$series"/MR.* 2>/dev/null | wc -l | tr -d ' ')
				if ((nraw == 0)); then
					ok=0; missingReasons+="  - no MR.* dicom files in $rawBase/$series${nl}"
				fi
			fi
		fi
	done

	if ((anyComm == 1)); then
		if [[ -z $commentScriptPath || ! -x $commentScriptPath ]]; then
			ok=0; missingReasons+="  - comment script not found or not executable: $commentScript${nl}"
		fi
		if [[ ! -e $indexFile ]]; then
			warnMsg+="  - warning: dicom index file not found ($indexFile); the comment script will use all dicoms in each series folder${nl}"
		fi
	fi

	if ((regenHdr == 1 || anyNii == 1 || anySide == 1 || anyComm == 1)) || [[ ! -e $subjBase/$hdrFileName ]]; then
		doHdr=1
	fi

	# ---- anatomical transfer, registration, segmentation transforms ----
	if ((nSeg > 0)); then
		for ((j = 0; j < nXfer; j++)); do
			name=${xferList[j]}; xferDo[j]=0
			if ((regenTransfer == 1)) || [[ ! -e $subjBase/$name ]]; then
				xferDo[j]=1
				base=$(strip_nii "$name")
				if [[ ! -e $transferSrc/${base}+orig.HEAD ]]; then
					ok=0; missingReasons+="  - source dataset for $name not found: $transferSrc/${base}+orig.HEAD${nl}"
				elif ! compgen -G "$transferSrc/${base}+orig.BRIK*" > /dev/null; then
					ok=0; missingReasons+="  - source dataset for $name has .HEAD but no .BRIK: $transferSrc/${base}+orig${nl}"
				fi
			fi
		done

		if ((regenRegister == 1)) || [[ ! -e $subjBase/$tformName ]]; then doRegister=1; fi
		if [[ $tformApplyName != "$tformName" ]]; then
			if ((doRegister == 1)) || [[ ! -e $subjBase/$tformApplyName ]]; then doInverse=1; fi
		fi
		if ((doRegister == 1)) && ! avail "$regBaseName"; then
			ok=0; missingReasons+="  - registration target $regBaseName not found in $subjBase (and not listed in transferNames)${nl}"
		fi

		for ((k = 0; k < nSeg; k++)); do
			name=${segArr[k]}; base=$(strip_nii "$name")
			idx=$(xfer_index "$name")
			segDo[k]=0
			if ((regenTransform == 1 || doRegister == 1)) || [[ ! -e $subjBase/${base}${segSuffix}.nii.gz ]] || ((xferDo[idx] == 1)); then
				segDo[k]=1
			fi
		done
	fi

	((ok == 1))
}

# ---------------------------------------------------------------------------
# stages (run with cwd = subjBase; set procErr and return 1 on failure)
# ---------------------------------------------------------------------------
# take_single_output <dir> <glob> <dest>: move the one file matching glob in dir to dest, remove dir
take_single_output() {
	local d=$1 g=$2 dest=$3 matches=()
	takeErr=""
	shopt -s nullglob
	matches=("$d"/$g)
	shopt -u nullglob
	if ((${#matches[@]} != 1)); then
		takeErr="expected 1 output matching $g, got ${#matches[@]}"
		((${#matches[@]} > 0)) && takeErr+=" ($(basename -a "${matches[@]}" | tr '\n' ' '))"
		rm -rf "$d"; return 1
	fi
	mv "${matches[0]}" "$dest" || { takeErr="could not move output to $dest"; rm -rf "$d"; return 1; }
	rm -rf "$d"
}

stage_nifti() {
	local ei series nii side niiBase sideBase nraw nframes zopt tmpOut
	for ((ei = 0; ei < nEchos; ei++)); do
		series=${echoNums[ei]}
		nii=$(subst_series "$rawNii" "$series")
		side=$(subst_series "$pattern_dicomHdr" "$series")

		if ((doNiiE[ei] == 1)); then
			echo "[$subjID] Converting series $series -> $nii"
			niiBase=$(strip_nii "$nii")
			nraw=$(ls -1 "$rawBase/$series"/MR.* 2>/dev/null | wc -l | tr -d ' ')
			rm -f "$subjBase/$nii" "$rawBase/$nii"
			# enhanced (multi-frame) dicoms hold many slices per file; Dimon can't assemble them
			nframes=1
			if command -v dicom_hinfo > /dev/null 2>&1; then
				nframes=$(dicom_hinfo -no_name -tag 0028,0008 "$(ls -1 "$rawBase/$series"/MR.* 2>/dev/null | head -n 1)" 2>/dev/null | tr -d ' ')
				is_int "$nframes" || nframes=1
			fi
			if ((nraw > dcm2niixMinFiles || nframes > 1)); then
				if ! command -v dcm2niix_afni > /dev/null 2>&1; then
					procErr+="  - dcm2niix_afni not found (needed for $nraw dicoms, $nframes frames/file)${nl}"; return 1
				fi
				zopt=""; [[ $nii == *.gz ]] && zopt="-z y"
				# dcm2niix may append suffixes (e.g. _e2 for multi-echo); write to a temp dir and rename
				tmpOut=$(mktemp -d "$subjBase/.dcm2niix.XXXXXX")
				if ! dcm2niix_afni -b n $zopt -f "$niiBase" -o "$tmpOut" "$rawBase/$series"; then
					rm -rf "$tmpOut"
					procErr+="  - dcm2niix_afni failed for $rawBase/$series${nl}"; return 1
				fi
				if ! take_single_output "$tmpOut" "*.nii*" "$subjBase/$nii"; then
					procErr+="  - dcm2niix_afni for $rawBase/$series: $takeErr${nl}"; return 1
				fi
			else
				if ! command -v Dimon > /dev/null 2>&1; then
					procErr+="  - Dimon not found (needed for $nraw dicoms)${nl}"; return 1
				fi
				if ! ( cd "$rawBase" && Dimon -infile_pattern "./$series/MR.*" -gert_create_dataset $dimonOpts -gert_to3d_prefix "$nii" ); then
					rm -f "$rawBase/$nii"
					procErr+="  - Dimon failed for $rawBase/$series${nl}"; return 1
				fi
				if [[ ! -e $rawBase/$nii ]] || ! mv "$rawBase/$nii" "$subjBase/$nii"; then
					procErr+="  - Dimon did not produce $rawBase/$nii${nl}"; return 1
				fi
			fi
			if [[ ! -e $subjBase/$nii ]]; then
				procErr+="  - expected nifti not found after conversion: $subjBase/$nii${nl}"; return 1
			fi
		fi

		if ((doSideE[ei] == 1)); then
			echo "[$subjID] Writing dcm2niix sidecar for series $series -> $side"
			sideBase=${side%.json}
			rm -f "$subjBase/$side"
			tmpOut=$(mktemp -d "$subjBase/.dcm2niix.XXXXXX")
			if ! dcm2niix -b o -f "$sideBase" -o "$tmpOut" "$rawBase/$series" > /dev/null; then
				rm -rf "$tmpOut"
				procErr+="  - dcm2niix (sidecar) failed for $rawBase/$series${nl}"; return 1
			fi
			if ! take_single_output "$tmpOut" "*.json" "$subjBase/$side"; then
				procErr+="  - dcm2niix (sidecar) for $rawBase/$series: $takeErr${nl}"; return 1
			fi
			if [[ ! -e $subjBase/$side ]]; then
				procErr+="  - sidecar not found after dcm2niix: $subjBase/$side (check pattern_dicomHdr)${nl}"; return 1
			fi
		fi

		# de-oblique (idempotent; same test as the old script)
		if 3dinfo "$nii" 2> /dev/null | grep -q 'Data Axes Tilt:.*Oblique'; then
			echo "[$subjID] De-obliquing $nii"
			if ! 3drefit -deoblique "$nii"; then
				procErr+="  - 3drefit -deoblique failed for $nii${nl}"; return 1
			fi
		fi
	done
	return 0
}

stage_comments() {
	local ei series dte
	for ((ei = 0; ei < nEchos; ei++)); do
		((doCommE[ei] == 1)) || continue
		series=${echoNums[ei]}
		dte=$(subst_series "$pattern_dteHdr" "$series")
		echo "[$subjID] Reading image comments for series $series -> $dte"
		if ! "$commentScriptPath" "$rawBase/$series" "$indexFile" "$series" "$scannerCode" "$subjBase/$dte"; then
			rm -f "$subjBase/$dte"
			procErr+="  - comment script failed for series $series${nl}"; return 1
		fi
	done
	return 0
}

# build_header <outfile>: one json for all echos (frames in acquisition/header order)
build_header() {
	local out=$1
	local ei series echoLabel dte side nii nLines ntNii
	local dteRaw zmRaw maxRaw stepRaw v nu
	local b0e fe te te_ms dteMs maxMs stepMs zmCsv proc sth sthLine
	local B0="" freq="" blocks="" block
	local firstDte="" firstZm=""
	hdrErr=""

	for ((ei = 0; ei < nEchos; ei++)); do
		series=${echoNums[ei]}; echoLabel=$((ei + 1))
		nii=$(subst_series "$rawNii" "$series")
		side="$subjBase/$(subst_series "$pattern_dicomHdr" "$series")"
		dte="$subjBase/$(subst_series "$pattern_dteHdr" "$series")"
		proc=$(subst_echo "$aseProcName" "$echoLabel")

		b0e=$(json_num "$side" MagneticFieldStrength)
		fe=$(json_num "$side" ImagingFrequency)
		te=$(json_num "$side" EchoTime)
		for v in "$b0e" "$fe" "$te"; do
			if ! is_num "$v"; then
				hdrErr+="  - missing/invalid MagneticFieldStrength, ImagingFrequency or EchoTime in $side${nl}"; return 1
			fi
		done
		if [[ -z $B0 ]]; then
			B0=$b0e; freq=$fe
		else
			[[ $b0e == "$B0" ]]  || warnMsg+="  - warning: MagneticFieldStrength for echo $echoLabel ($b0e) differs from echo 1 ($B0)${nl}"
			[[ $fe == "$freq" ]] || warnMsg+="  - warning: ImagingFrequency for echo $echoLabel ($fe) differs from echo 1 ($freq)${nl}"
		fi
		te_ms=$(awk -v t="$te" 'BEGIN{printf "%.10g", t*1000}')
		# slice thickness (mm; needed for z-shim correction). Optional: omitted from the header if absent.
		sth=$(json_num "$side" SliceThickness)
		sthLine=""
		if is_num "$sth"; then
			sthLine=$(printf '\n      "SliceThickness_mm": %s,' "$sth")
		else
			warnMsg+="  - warning: SliceThickness not found in $side (z-shim correction will be unavailable)${nl}"
		fi

		nLines=$(tr -d '\r' < "$dte" | grep -c '_DTE=')
		if ((nLines == 0)); then hdrErr+="  - no '_DTE=' entries found in $dte${nl}"; return 1; fi
		ntNii=$(3dinfo -nt "$nii" 2> /dev/null)   # relative: AFNI splits paths on spaces
		if [[ ! $ntNii =~ ^[0-9]+$ ]] || ((ntNii != nLines)); then
			hdrErr+="  - frame count mismatch for series $series: $nii has '$ntNii' frames but $dte has $nLines rows${nl}"; return 1
		fi
		dteRaw=$(dte_field "$dte" "$dteExpr")
		zmRaw=$(dte_field "$dte" "$zmExpr")
		maxRaw=$(dte_field "$dte" "$maxExpr")
		stepRaw=$(dte_field "$dte" "$stepExpr")
		for v in "$dteRaw" "$zmRaw" "$maxRaw" "$stepRaw"; do
			if [[ $(count_lines "$v") -ne $nLines ]]; then
				hdrErr+="  - could not parse DTE/ZMoment/MaxDTE/dTE from every line of $dte${nl}"; return 1
			fi
		done
		for v in $dteRaw $zmRaw $maxRaw $stepRaw; do
			if ! is_num "$v"; then hdrErr+="  - non-numeric value '$v' parsed from $dte${nl}"; return 1; fi
		done
		nu=$(printf '%s\n' "$maxRaw"  | sort -u | grep -c .)
		((nu > 1)) && warnMsg+="  - warning: MaxDTE varies between lines of $dte (using first)${nl}"
		nu=$(printf '%s\n' "$stepRaw" | sort -u | grep -c .)
		((nu > 1)) && warnMsg+="  - warning: dTE step varies between lines of $dte (using first)${nl}"

		dteMs=$(printf '%s\n' "$dteRaw" | awk '{printf "%.10g\n", $1/1000}' | paste -sd, -)
		maxMs=$(printf '%s\n' "$maxRaw"  | head -n 1 | awk '{printf "%.10g", $1/1000}')
		stepMs=$(printf '%s\n' "$stepRaw" | head -n 1 | awk '{printf "%.10g", $1/1000}')
		zmCsv=$(printf '%s\n' "$zmRaw" | paste -sd, -)

		if [[ -z $firstDte ]]; then
			firstDte=$dteMs; firstZm=$zmCsv
		else
			[[ $dteMs == "$firstDte" ]] || warnMsg+="  - warning: dTE values for echo $echoLabel differ from echo 1${nl}"
			[[ $zmCsv == "$firstZm" ]]  || warnMsg+="  - warning: ZMoment values for echo $echoLabel differ from echo 1${nl}"
		fi

		block=$(printf '    {\n      "echo": %d,\n      "series": "%s",\n      "EchoTime_ms": %s,%s\n      "nFrames": %d,\n      "MaxDTE_ms": %s,\n      "StepDTE_ms": %s,\n      "dTE_ms": [%s],\n      "ZMoment": [%s],\n      "file": "%s"\n    }' \
			"$echoLabel" "$series" "$te_ms" "$sthLine" "$nLines" "$maxMs" "$stepMs" "$dteMs" "$zmCsv" "$proc")
		blocks+="${blocks:+,$nl}$block"
	done

	{
		printf '{\n'
		printf '  "subject": "%s",\n' "$subjID"
		printf '  "created": "%s",\n' "$runTimestamp"
		printf '  "units": {"MagneticFieldStrength": "T", "ImagingFrequency": "MHz", "EchoTime": "ms", "dTE": "ms", "MaxDTE": "ms", "StepDTE": "ms", "ZMoment": "Tms/m"},\n'
		printf '  "frame_order": "acquisition (header row) order; frame i of each echo file matches dTE_ms[i] and ZMoment[i]; not sorted",\n'
		printf '  "MagneticFieldStrength": %s,\n' "$B0"
		printf '  "ImagingFrequency": %s,\n' "$freq"
		printf '  "nEchos": %d,\n' "$nEchos"
		printf '  "echos": [\n%s\n  ]\n' "$blocks"
		printf '}\n'
	} > "$out.tmp" && mv "$out.tmp" "$out"
}

stage_header() {
	((doHdr == 1)) || return 0
	echo "[$subjID] Building header: $hdrFileName"
	if ! build_header "$subjBase/$hdrFileName"; then
		rm -f "$subjBase/$hdrFileName.tmp"
		procErr+=$hdrErr; return 1
	fi
	return 0
}

# base volume for motion correction / registration: frame <baseFrame> of echo 1 (created only if missing)
ensure_base() {
	local nii1
	[[ -e $baseName ]] && return 0
	nii1=$(subst_series "$rawNii" "${echoNums[0]}")
	echo "[$subjID] Creating $baseName from $nii1[$baseFrame]"
	if ! 3dcalc -a "${nii1}[${baseFrame}]" -expr a -prefix "$baseName"; then
		rm -f "$baseName"
		procErr+="  - could not create $baseName from $nii1${nl}"; return 1
	fi
	return 0
}

stage_motion() {
	local ei series nii proc mc
	((anyMotion == 1)) || return 0
	ensure_base || return 1
	for ((ei = 0; ei < nEchos; ei++)); do
		((doMotionE[ei] == 1)) || continue
		series=${echoNums[ei]}
		nii=$(subst_series "$rawNii" "$series")
		proc=$(subst_echo "$aseProcName" "$((ei + 1))")
		mc=$(subst_echo "$motionParamName" "$((ei + 1))")
		if 3dinfo -same_grid "$baseName" "$nii" 2> /dev/null | grep -qw 0; then
			warnMsg+="  - warning: $nii is not on the same grid as the existing $baseName; existing registration may not match${nl}"
		fi
		echo "[$subjID] Motion correcting echo $((ei + 1)) -> $proc"
		rm -f "$mc"
		if ! 3dvolreg $volregOpts -base "$baseName" -dfile "$mc" -prefix "$proc" -overwrite "$nii"; then
			rm -f "$proc" "$mc"
			procErr+="  - 3dvolreg failed for $nii${nl}"; return 1
		fi
		if [[ ! -e $proc ]]; then
			procErr+="  - expected output not found after 3dvolreg: $proc${nl}"; return 1
		fi
	done
	return 0
}

stage_transfer() {
	local j name base
	((nSeg > 0)) || return 0
	for ((j = 0; j < nXfer; j++)); do
		((xferDo[j] == 1)) || continue
		name=${xferList[j]}; base=$(strip_nii "$name")
		echo "[$subjID] Transferring $name from $transferSrcDir"
		if ! ( cd "$subjDir" && 3dAFNItoNIFTI -overwrite -prefix "$aseSubjDir/$name" "$transferSrcDir/${base}+orig" ); then
			rm -f "$subjBase/$name"
			procErr+="  - 3dAFNItoNIFTI failed for $transferSrcDir/${base}+orig${nl}"; return 1
		fi
		if [[ ! -e $subjBase/$name ]]; then
			procErr+="  - expected file not found after transfer: $subjBase/$name${nl}"; return 1
		fi
	done
	return 0
}

stage_register() {
	((nSeg > 0)) || return 0
	if ((doRegister == 1)); then
		ensure_base || return 1
		echo "[$subjID] Registering $baseName to $regBaseName"
		rm -f "$tformName"
		if ! 3dAllineate -prefix "$baseAlName" -base "$regBaseName" -master "$regBaseName" -source "$baseName" \
				$allineateOpts -1Dmatrix_save "$tformName" -overwrite; then
			rm -f "$tformName" "$baseAlName"
			procErr+="  - 3dAllineate registration failed${nl}"; return 1
		fi
		if [[ ! -s $tformName ]]; then
			procErr+="  - transform file not produced: $tformName${nl}"; return 1
		fi
	fi
	if ((doInverse == 1)); then
		echo "[$subjID] Creating inverse transform: $tformApplyName"
		rm -f "$tformApplyName"
		if ! cat_matvec "$tformName" -I > "$tformApplyName" || [[ ! -s $tformApplyName ]]; then
			rm -f "$tformApplyName"
			procErr+="  - cat_matvec failed for $tformName${nl}"; return 1
		fi
	fi
	return 0
}

stage_segments() {
	local k seg base out t0 t1
	((nSeg > 0)) || return 0
	for ((k = 0; k < nSeg; k++)); do
		((segDo[k] == 1)) || continue
		seg=${segArr[k]}; base=$(strip_nii "$seg"); out=${base}${segSuffix}.nii.gz
		echo "[$subjID] Transforming $seg -> $out (matrix: $tformApplyName)"
		t0=$(date +%s)
		if ! 3dAllineate -1Dmatrix_apply "$tformApplyName" -source "$seg" -master "$seg" -prefix "$out" -final NN -overwrite; then
			rm -f "$out"
			procErr+="  - 3dAllineate failed for $seg${nl}"; return 1
		fi
		if [[ ! -e $out ]]; then
			procErr+="  - expected output not found after 3dAllineate: $out${nl}"; return 1
		fi
		t1=$(date +%s)
		echo "[$subjID]   done in $((t1 - t0))s"
	done
	return 0
}

process_subject() {
	procErr=""
	mkdir -p "$subjBase" && cd "$subjBase" || { procErr+="  - cannot create/enter $subjBase${nl}"; return 1; }
	stage_nifti    || return 1
	stage_comments || return 1
	stage_header   || return 1
	stage_motion   || return 1
	stage_transfer || return 1
	stage_register || return 1
	stage_segments || return 1
	return 0
}

# ---------------------------------------------------------------------------
# subject list, retry, logging
# ---------------------------------------------------------------------------
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
if ((omitFailLog == 0)); then
	echo "Run started: $(date)" > "$failFile"
	echo "Completed subjects:" >> "$failFile"
fi

failLog=""; nfail=0; nok=0

echo "=================================================="
echo "Preprocessing ASE data for $nsubj subject(s)..."
subj_i=0
for subjID in "${subjlist[@]}"; do
	((subj_i++))
	if ((${#completedSubjs[@]} > 0)); then
		if printf '%s\n' "${completedSubjs[@]}" | grep -qxF "$subjID"; then
			echo "Skipping $subjID (already completed)"
			log_raw "+$subjID${nl}"
			continue
		fi
	fi

	subjDir="$allSubjDir/$subjID"
	subjBase="$subjDir/$aseSubjDir"
	rawBase="$subjDir/$rawDataDir"
	indexFile="$subjDir/$dcmIndexFile"
	transferSrc="$subjDir/$transferSrcDir"
	warnMsg=""

	if ! check_subject; then
		record_fail "Skipping $subjID:${nl}${missingReasons}"
		continue
	fi
	if ! process_subject; then
		record_fail "Failed $subjID:${nl}${procErr}"
		continue
	fi

	if [[ -n $warnMsg ]]; then
		printf 'Warnings for %s:\n%s' "$subjID" "$warnMsg"
		log_raw "Warnings for $subjID:${nl}${warnMsg}"
	fi
	log_raw "+$subjID${nl}"
	((nok++))
	echo "Finished $subjID ($subj_i/$nsubj)"
done
cd "$allSubjDir" || exit 1

echo "=================================================="
echo "Done: $nok succeeded, $nfail failed/skipped."
if ((nfail > 0)); then printf '\nFailure summary:\n%s' "$failLog"; fi
if ((omitFailLog == 0)); then echo "Log: $failFile"; fi
((nfail == 0))