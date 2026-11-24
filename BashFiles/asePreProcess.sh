#!/bin/bash

# asePreProcess.sh: prepares ASE subject directories for fitting, starting from raw DICOMs.
# Per subject, per echo series:
#   1. convert dicoms -> 4D nifti (all frames, acquisition order) + dcm2niix json sidecar; de-oblique
#   2. build the dTE / z-moment / image-comment header (separate script, see commentScript)
#   3. build the combined json header (all echos)
#   4. motion correct -> <aseProcName>. motionMode=echo1 (default): 3dvolreg on echo 1 only, and
#      the same per-frame transforms applied to the other echoes (all echoes of a frame come from
#      the same excitation, so they share its head position; avoids cross-TE registration).
#      motionMode=perEcho: 3dvolreg on each echo against the echo-1 base (legacy).
#   5. anatomy. anatSource=fsrecon (default): convert the FreeSurfer volumes (fsVolumes) from
#      <fsSubjDir>/mri with mri_convert, and build the registration target (T1 x step(aseg),
#      3dUnifize) and the xvent tissue mask (aparc+aseg > 0 minus the ventricles) on FreeSurfer's
#      own grid. anatSource=afni (legacy): copy the BOLD-aligned "_al" volumes from transferSrcDir
#      (+orig -> nifti). See docs/anat-route-report.md.                                [if missing]
#   6. register ase_base to the anatomical (3dAllineate, registerDof parameters) and make the
#      inverse matrix. The settings are recorded next to the matrix; a change re-registers.
#   7. transform segmentations into ASE space (nearest neighbor)
#   7b. optional (anatSource=afni only): build a tissue fraction (tissueFracName) from SynthSeg
#       posteriors (tissuePosteriorName): 1 - sum(background, CSF, ventricle posteriors) [if missing]
#   8. resample quality masks (qualityMaskNames) onto the ASE grid (linear -> fractions 0..1)
# Existing outputs are skipped unless regenerated. Regenerating an upstream stage forces its
# downstream stages (and overwrites). Settings: script defaults < config "pipeline" section < flags.
# Stages 5-8 are skipped when segNames is empty.

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
rawNii="ase_%SERIES%.nii"
pattern_dicomHdr="dcm_hdr_%SERIES%.json"
pattern_dteHdr="dcm_hdr_%SERIES%_dtes.txt"
aseProcName="ase_e%ECHO%_proc.nii.gz"
motionParamName="ase_e%ECHO%_mc.1D"
hdrFileName=ase_hdr.json
volregOpts="-Fourier -twopass -twodup"
motionMode="echo1"                  # echo1 | perEcho (see stage 4)
motionMatName="ase_e1_mc_mat.aff12.1D"   # echo-1 per-frame transforms (motionMode=echo1)
motionApplyInterp="wsinc5"          # 3dAllineate -final for echoes > 1 (motionMode=echo1)
baseName=ase_base.nii.gz
baseFrame=0
baseAlName=ase_base_al.nii.gz
# ---- anatomy (stage 5) ----
anatSource=fsrecon                  # fsrecon | afni (legacy: BOLD-aligned volumes from transferSrcDir)
# fsrecon: FreeSurfer subject folder, subject-relative or absolute; %SUBJ% is replaced by the subject ID
fsSubjDir=fsrecon
# fsrecon: <file in fsSubjDir/mri>=<name in aseSubjDir>, comma-separated. T1.mgz and aseg.mgz make the
# registration target, aparc+aseg.mgz the xvent mask; the others are only segmentations.
fsVolumes="T1.mgz=anat_fs_T1.nii.gz,aseg.mgz=anat_fs_aseg.nii.gz,aparc+aseg.mgz=anat_fs_aparc_aseg.nii.gz,wmparc.mgz=anat_fs_wmparc.nii.gz,brainstemSsLabels.v13.FSvoxelSpace.mgz=anat_fs_bstem.nii.gz"
xventName=anat_fs_xvent_mask.nii.gz          # fsrecon: aparc+aseg > 0 minus xventExcludeLabels
xventExcludeLabels="4,5,14,15,43,44"         # lateral, inferior-lateral, 3rd and 4th ventricles (as the old xvent)
# "auto" = the names for anatSource (fsrecon: anat_fs_brain, anat_fs_aparc_aseg, anat_fs_wmparc,
# anat_fs_bstem, anat_fs_xvent_mask; afni: the "_al" names). An empty segNames means no ROI stats.
regBaseName=auto
segNames=auto
qualityMaskNames=auto                        # tissue fractions resampled to the ASE grid; empty = none
registerDof=12                      # ase_base -> anatomy: 12 (affine, default), 9 (+ scale), 6 (rigid)
allineateOpts="-twopass -cost lpa+ -autoweight -fineblur 1 -cmass"   # no -warp here: use registerDof
tformName=ase_base_al_mat.aff12.1D
tformApplyName=ase_base_al_inv_mat.aff12.1D
# afni (legacy) only:
transferSrcDir=session_afni_nortc
transferNames="anat_fs_brain_al.nii.gz,anat_fs_brain_al_mask.nii.gz,anat_fs_xvent_mask_al.nii.gz,anat_fs_bstem_al.nii.gz"
segSuffix="_ase"
# Optional SynthSeg tissue fraction (off by default). To use it instead of xvent, set
#   tissuePosteriorName=anat_fs_al_posteriors.nii.gz  and  qualityMaskNames=anat_fs_tissueFrac_al.nii.gz
# (and qualityMaskNames for aseFit.sh). Note SynthSeg's cortex extends into FreeSurfer's sulcal CSF.
tissuePosteriorName=""                             # SynthSeg posteriors (NIfTI) in transferSrcDir; empty = off (afni only)
tissueFracName="anat_fs_tissueFrac_al.nii.gz"      # built from the posteriors (anatomical grid)
posteriorNonTissueVols="0,3,4,11,12,16,21,22"      # 0-based sub-bricks: background, L lat/inf-lat vent,
                                                   #   3rd, 4th vent, CSF (label 24), R lat/inf-lat vent
posteriorNVols=33                                  # expected sub-bricks (SynthSeg standard label set)

pipelineKeys="allSubjDir dataList echo1Col echoSpacing nEchos aseSubjDir rawDataDir dcmIndexFile scannerCode commentScript rawNii pattern_dicomHdr pattern_dteHdr aseProcName motionParamName hdrFileName volregOpts motionMode motionMatName motionApplyInterp baseName baseFrame baseAlName anatSource fsSubjDir fsVolumes xventName xventExcludeLabels registerDof regBaseName allineateOpts tformName tformApplyName transferSrcDir transferNames segNames segSuffix qualityMaskNames tissuePosteriorName tissueFracName posteriorNonTissueVols posteriorNVols"

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
  --regenAnat           redo stage 5: FreeSurfer conversion (fsrecon) or transfer (afni); forces
                         registration and segmentation transforms. --regenTransfer is the old name.
  --regenRegister        re-register ase_base to anatomical (forces segmentation transforms)
  --regenTransform       segmentation transforms and quality-mask resampling only

Logging / retry:
  --omitFailLog | --retryMostRecent | --retryFromFile <file> | --help
EOF
}

need_arg() {
	if [[ $# -lt 2 ]]; then echo "Error: $1 requires an argument" >&2; exit 1; fi
}

# keys that no longer do anything (2026-11-20: Dimon removed, dcm2niix converts every series)
removedKeys="dcm2niixMinFiles dimonOpts"

set_key() {
	local kv=$1 k v
	if [[ $kv != *=* ]]; then echo "Error: --set expects key=value (got: $kv)" >&2; exit 1; fi
	k=${kv%%=*}; v=${kv#*=}
	case " $removedKeys " in
		*" $k "*) echo "Warning: pipeline key $k was removed (dcm2niix now converts every series); ignored" >&2; return 0 ;;
	esac
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
	for key in $removedKeys; do
		if cfg_get "$pipeSec" "$key" > /dev/null; then
			echo "Warning: config key $key was removed (dcm2niix now converts every series); ignored" >&2
		fi
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
		--regenTransfer|--regenAnat) regenTransfer=1; shift ;;
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
for k in nEchos echo1Col echoSpacing scannerCode baseFrame; do
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

# anatomy source and the names that depend on it
case $anatSource in
	fsrecon)
		[[ $regBaseName == auto ]] && regBaseName=anat_fs_brain.nii.gz
		[[ $segNames == auto ]] && segNames="anat_fs_aparc_aseg.nii.gz,anat_fs_wmparc.nii.gz,anat_fs_bstem.nii.gz"
		[[ $qualityMaskNames == auto ]] && qualityMaskNames=$xventName
		if [[ -n $tissuePosteriorName ]]; then
			echo "Error: tissuePosteriorName (SynthSeg tissue fraction) is only available with anatSource=afni" >&2; exit 1
		fi
		if [[ ! $xventExcludeLabels =~ ^[0-9]+(,[0-9]+)*$ ]]; then
			echo "Error: xventExcludeLabels must be comma-separated label numbers (got: $xventExcludeLabels)" >&2; exit 1
		fi ;;
	afni)
		[[ $regBaseName == auto ]] && regBaseName=anat_fs_brain_al.nii.gz
		[[ $segNames == auto ]] && segNames="anat_fs_aparc_aseg_al.nii.gz,anat_fs_wmparc_al.nii.gz"
		[[ $qualityMaskNames == auto ]] && qualityMaskNames=anat_fs_xvent_mask_al.nii.gz ;;
	*) echo "Error: anatSource must be fsrecon or afni (got: $anatSource)" >&2; exit 1 ;;
esac
case $registerDof in
	12) regWarp=affine_general ;;
	9)  regWarp=shift_rotate_scale ;;
	6)  regWarp=shift_rotate ;;
	*)  echo "Error: registerDof must be 12, 9 or 6 (got: $registerDof)" >&2; exit 1 ;;
esac
if [[ " $allineateOpts " == *" -warp "* ]]; then
	echo "Error: allineateOpts contains -warp; set registerDof instead" >&2; exit 1
fi
# what stage 6 was run with; recorded next to the matrix so a change forces a new registration
regSettings="anatSource=$anatSource regBaseName=$regBaseName registerDof=$registerDof allineateOpts=$allineateOpts"
regStampName="${tformName%.aff12.1D}_settings.txt"

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

gammaBarMHz=42.577478518      # proton gyromagnetic ratio / 2pi (MHz/T), CODATA 2018
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
segArr=(); transferArr=(); qmArr=()
[[ -n $segNames ]] && IFS=',' read -ra segArr <<< "$segNames"
[[ -n $transferNames ]] && IFS=',' read -ra transferArr <<< "$transferNames"
[[ -n $qualityMaskNames ]] && IFS=',' read -ra qmArr <<< "$qualityMaskNames"
nSeg=${#segArr[@]}
nQm=${#qmArr[@]}
# stage-5 items: xferList = names in aseSubjDir; xferSrcList = their source (afni: +orig dataset
# name; fsrecon: file in <fsSubjDir>/mri)
xferList=(); xferSrcList=()
fsT1="" fsAseg="" fsAparc=""
if ((nSeg > 0)) && [[ $anatSource == afni ]]; then
	for name in "${transferArr[@]}" "${segArr[@]}" "${qmArr[@]}"; do
		[[ -n $tissuePosteriorName && $name == "$tissueFracName" ]] && continue   # built, not transferred
		dup=0
		for e in "${xferList[@]}"; do [[ $e == "$name" ]] && dup=1; done
		((dup)) || { xferList+=("$name"); xferSrcList+=("$(strip_nii "$name")+orig"); }
	done
elif ((nSeg > 0)); then
	IFS=',' read -ra fsItems <<< "$fsVolumes"
	for item in "${fsItems[@]}"; do
		if [[ $item != *=* || -z ${item%%=*} || -z ${item#*=} ]]; then
			echo "Error: fsVolumes entries must be <file>=<name> (got: $item)" >&2; exit 1
		fi
		xferSrcList+=("${item%%=*}"); xferList+=("${item#*=}")
		case ${item%%=*} in
			T1.mgz)         fsT1=${item#*=} ;;
			aseg.mgz)       fsAseg=${item#*=} ;;
			aparc+aseg.mgz) fsAparc=${item#*=} ;;
		esac
	done
	if [[ -z $fsT1 || -z $fsAseg ]]; then
		echo "Error: fsVolumes must include T1.mgz and aseg.mgz (they make the registration target)" >&2; exit 1
	fi
fi
nXfer=${#xferList[@]}
# fsrecon: the xvent mask is built when a later stage uses it
needXvent=0
if [[ $anatSource == fsrecon ]] && ((nSeg > 0)); then
	for e in "${qmArr[@]}" "${segArr[@]}"; do [[ $e == "$xventName" ]] && needXvent=1; done
	if ((needXvent == 1)) && [[ -z $fsAparc ]]; then
		echo "Error: $xventName is used but fsVolumes has no aparc+aseg.mgz entry" >&2; exit 1
	fi
fi

xfer_index() { local j; for ((j = 0; j < nXfer; j++)); do [[ ${xferList[j]} == "$1" ]] && { printf '%d' "$j"; return 0; }; done; return 1; }

# xfer_required <name>: a later stage reads it (registration target, segmentation, quality mask).
# Other transferNames entries are optional: a missing source is a warning, not a failure.
# fsrecon: the target's inputs (T1, aseg) and, when the xvent mask is used, aparc+aseg are required;
# a missing segmentation is left out with a warning (as aseFit.sh does).
xfer_required() {
	local e
	if [[ $anatSource == fsrecon ]]; then
		[[ $1 == "$fsT1" || $1 == "$fsAseg" ]] && return 0
		((needXvent == 1)) && [[ $1 == "$fsAparc" ]] && return 0
		for e in "${qmArr[@]}"; do [[ $e == "$1" ]] && return 0; done
		return 1
	fi
	[[ $1 == "$regBaseName" ]] && return 0
	for e in "${segArr[@]}" "${qmArr[@]}"; do [[ $e == "$1" ]] && return 0; done
	return 1
}

# avail <name>: file exists in subjBase or is scheduled to be transferred
avail() {
	local j
	[[ -e $subjBase/$1 ]] && return 0
	[[ $anatSource == fsrecon && $1 == "$regBaseName" && $brainDo == 1 ]] && return 0
	[[ $anatSource == fsrecon && $1 == "$xventName" && $xventDo == 1 ]] && return 0
	[[ -n $tissuePosteriorName && $1 == "$tissueFracName" && $tfDo == 1 ]] && return 0
	for ((j = 0; j < nXfer; j++)); do
		[[ ${xferList[j]} == "$1" && ${xferDo[j]} == 1 ]] && return 0
	done
	return 1
}

# tools
needTools="dcm2niix dcm2niix_afni 3dinfo 3drefit 3dvolreg 3dcalc"
[[ $motionMode == echo1 ]] && ((nEchos > 1)) && needTools+=" 3dAllineate"
((nSeg > 0)) && needTools+=" 3dAllineate cat_matvec"
((nSeg > 0)) && [[ $anatSource == afni ]] && needTools+=" 3dAFNItoNIFTI"
((nSeg > 0)) && [[ $anatSource == fsrecon ]] && needTools+=" mri_convert 3dUnifize"
for t in $needTools; do
	if ! command -v "$t" > /dev/null 2>&1; then
		echo "Error: $t not found on PATH" >&2
		[[ $t == mri_convert ]] && echo "       anatSource=fsrecon needs FreeSurfer (source \$FREESURFER_HOME/SetUpFreeSurfer.sh), or set anatSource=afni" >&2
		exit 1
	fi
done
echo "Anatomy: anatSource=$anatSource, registration $registerDof parameters ($regWarp)"

# AFNI build in use, and whether its 3dAllineate knows the registration cost (old builds,
# e.g. 18.0.22, lack lpa+ and would otherwise fail on every subject)
afniVer=$(3dvolreg -ver 2>&1 | head -n 1 | sed "s/^++ 3dvolreg: //")
echo "AFNI: $(dirname "$(command -v 3dvolreg)") (${afniVer:-version unknown})"
if ((nSeg > 0)) && [[ $allineateOpts =~ -cost[[:space:]]+([^[:space:]]+) ]]; then
	allineateCost=${BASH_REMATCH[1]}
	if ! 3dAllineate -HELP 2>/dev/null | grep -qF -- "$allineateCost"; then
		echo "Error: this AFNI's 3dAllineate does not know the cost '$allineateCost' (from allineateOpts)." >&2
		echo "       Put a newer AFNI first on PATH, e.g. /usr/pubsw/packages/AFNI/26.0.10 on the cluster." >&2
		exit 1
	fi
fi

# ---------------------------------------------------------------------------
# planning: which stages run for a subject, and are their inputs available?
# ---------------------------------------------------------------------------
check_subject() {
	local ei series nii side dte proc mc nraw j k name base idx srcErr
	missingReasons=""; ok=1
	echoNums=(); doNiiE=(); doSideE=(); doCommE=(); doMotionE=(); xferDo=(); segDo=(); qmDo=()
	anyNii=0; anySide=0; anyComm=0; anyMotion=0
	doRegister=0; doInverse=0; doHdr=0; tfDo=0; brainDo=0; xventDo=0

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

	if [[ $motionMode == echo1 ]] && ((anyMotion == 1)); then
		# echoes > 1 reuse echo 1's transforms, so all echoes are redone together. Existing outputs
		# (e.g. from motionMode=perEcho) are left alone unless regenerated with --regenMotion.
		for ((ei = 0; ei < nEchos; ei++)); do doMotionE[ei]=1; done
	fi

	if ((anyComm == 1)); then
		if [[ -z $commentScriptPath || ! -x $commentScriptPath ]]; then
			ok=0; missingReasons+="  - comment script not found or not executable: $commentScript${nl}"
		fi
	fi

	if ((regenHdr == 1 || anyNii == 1 || anySide == 1 || anyComm == 1)) || [[ ! -e $subjBase/$hdrFileName ]]; then
		doHdr=1
	fi

	# ---- anatomical transfer, registration, segmentation transforms ----
	if ((nSeg > 0)); then
		local fsMissing=0
		if [[ $anatSource == fsrecon && ! -d $fsBase/mri ]]; then fsMissing=1; fi
		for ((j = 0; j < nXfer; j++)); do
			name=${xferList[j]}; xferDo[j]=0
			if ((regenTransfer == 1)) || [[ ! -e $subjBase/$name ]]; then
				xferDo[j]=1
				srcErr=""
				if [[ $anatSource == fsrecon ]]; then
					if ((fsMissing == 1)); then
						srcErr="FreeSurfer subject folder not found: $fsBase/mri (fsSubjDir=$fsSubjDir)"
						[[ $missingReasons == *"$srcErr"* ]] && srcErr="(same)"
					elif [[ ! -e $fsBase/mri/${xferSrcList[j]} ]]; then
						srcErr="FreeSurfer volume for $name not found: $fsBase/mri/${xferSrcList[j]}"
					fi
				elif [[ ! -e $transferSrc/${xferSrcList[j]}.HEAD ]]; then
					srcErr="source dataset for $name not found: $transferSrc/${xferSrcList[j]}.HEAD"
				elif ! compgen -G "$transferSrc/${xferSrcList[j]}.BRIK*" > /dev/null; then
					srcErr="source dataset for $name has .HEAD but no .BRIK: $transferSrc/${xferSrcList[j]}"
				fi
				if [[ -n $srcErr ]]; then
					if xfer_required "$name"; then
						ok=0; [[ $srcErr == "(same)" ]] || missingReasons+="  - $srcErr${nl}"
					elif [[ $srcErr == "(same)" || $fsMissing == 1 ]]; then
						xferDo[j]=0
					else
						# nothing later needs it (or it is a segmentation, left out of the ROIs): skip it;
						# rerun once the source exists
						xferDo[j]=0
						warnMsg+="  - warning: optional anatomy skipped, $srcErr${nl}"
					fi
				fi
			fi
		done

		# fsrecon: registration target and xvent mask, built on FreeSurfer's grid from the converted volumes
		if [[ $anatSource == fsrecon ]]; then
			if ((regenTransfer == 1 || xferDo[$(xfer_index "$fsT1")] == 1 || xferDo[$(xfer_index "$fsAseg")] == 1)) || \
					[[ ! -e $subjBase/$regBaseName ]]; then
				brainDo=1
			fi
			if ((needXvent == 1)); then
				if ((regenTransfer == 1 || xferDo[$(xfer_index "$fsAparc")] == 1)) || [[ ! -e $subjBase/$xventName ]]; then
					xventDo=1
				fi
			fi
		fi

		if ((regenRegister == 1 || brainDo == 1)) || [[ ! -e $subjBase/$tformName ]]; then doRegister=1; fi
		# a matrix made with other settings (e.g. the old BOLD-aligned anatomy) is not reused
		if ((doRegister == 0)); then
			if [[ ! -e $subjBase/$regStampName ]]; then
				if [[ $anatSource == fsrecon ]]; then
					doRegister=1
					warnMsg+="  - note: existing $tformName has no record of its settings; re-registering to the FreeSurfer anatomy${nl}"
				fi
			elif [[ $(cat "$subjBase/$regStampName") != "$regSettings" ]]; then
				doRegister=1
				warnMsg+="  - note: registration settings changed since $tformName was made; re-registering${nl}"
			fi
		fi
		if [[ $tformApplyName != "$tformName" ]]; then
			if ((doRegister == 1)) || [[ ! -e $subjBase/$tformApplyName ]]; then doInverse=1; fi
		fi
		if ((doRegister == 1)) && ! avail "$regBaseName"; then
			ok=0; missingReasons+="  - registration target $regBaseName not found in $subjBase (and not listed in transferNames)${nl}"
		fi

		for ((k = 0; k < nSeg; k++)); do
			name=${segArr[k]}; base=$(strip_nii "$name")
			segDo[k]=0
			if ! avail "$name"; then
				# fsrecon: left out of the ROI table (aseFit.sh warns too); afni: caught as a missing transfer
				[[ $anatSource == fsrecon ]] && warnMsg+="  - warning: segmentation $name not available; left out${nl}"
				continue
			fi
			srcNew=0
			if [[ $anatSource == fsrecon && $name == "$xventName" ]]; then srcNew=$xventDo
			elif idx=$(xfer_index "$name"); then srcNew=${xferDo[idx]}; fi
			if ((regenTransform == 1 || doRegister == 1 || srcNew == 1)) || [[ ! -e $subjBase/${base}${segSuffix}.nii.gz ]]; then
				segDo[k]=1
			fi
		done

		# tissue fraction from the SynthSeg posteriors (rebuilt with the transfers)
		if [[ -n $tissuePosteriorName ]]; then
			if ((regenTransfer == 1)) || [[ ! -e $subjBase/$tissueFracName ]]; then
				tfDo=1
				if [[ ! -e $transferSrc/$tissuePosteriorName ]]; then
					ok=0; missingReasons+="  - tissue posteriors not found: $transferSrc/$tissuePosteriorName${nl}"
				fi
			fi
		fi

		# quality masks: same matrix, but onto the ASE grid (master = baseName)
		for ((k = 0; k < nQm; k++)); do
			name=${qmArr[k]}; base=$(strip_nii "$name")
			srcNew=0
			if [[ -n $tissuePosteriorName && $name == "$tissueFracName" ]]; then
				srcNew=$tfDo
			elif [[ $anatSource == fsrecon && $name == "$xventName" ]]; then
				srcNew=$xventDo
			elif idx=$(xfer_index "$name"); then
				srcNew=${xferDo[idx]}
			fi
			if ! avail "$name"; then
				ok=0; missingReasons+="  - quality mask $name not available (not in fsVolumes / transferNames, not built)${nl}"
				continue
			fi
			qmDo[k]=0
			if ((regenTransform == 1 || doRegister == 1)) || [[ ! -e $subjBase/${base}${segSuffix}.nii.gz ]] || ((srcNew == 1)); then
				qmDo[k]=1
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
	local ei series nii side niiBase sideBase zopt tmpOut
	for ((ei = 0; ei < nEchos; ei++)); do
		series=${echoNums[ei]}
		nii=$(subst_series "$rawNii" "$series")
		side=$(subst_series "$pattern_dicomHdr" "$series")

		if ((doNiiE[ei] == 1)); then
			echo "[$subjID] Converting series $series -> $nii"
			niiBase=$(strip_nii "$nii")
			rm -f "$subjBase/$nii" "$rawBase/$nii"
			# dcm2niix for every series: it handles classic, mosaic and enhanced (multi-frame)
			# DICOMs and is needed anyway for the json sidecar. (The old cluster script used Dimon
			# for small series; on XA30 data both give identical pixels and only the header origin
			# differs: docs/figures/preproc_registration/reg1_conversion.png.)
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
		# the script's own messages go to the log: its errors on failure, its notes as warnings
		cmtLog=$(mktemp)
		if ! "$commentScriptPath" "$rawBase/$series" "$indexFile" "$series" "$scannerCode" "$subjBase/$dte" > "$cmtLog" 2>&1; then
			cat "$cmtLog"
			rm -f "$subjBase/$dte"
			procErr+="  - comment script failed for series $series:${nl}$(sed 's/^/      /' "$cmtLog")${nl}"
			rm -f "$cmtLog"; return 1
		fi
		cat "$cmtLog"
		while IFS= read -r line; do
			warnMsg+="  - note: ${line#Note: }${nl}"
		done < <(grep '^Note:' "$cmtLog")
		rm -f "$cmtLog"
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
			# echo k's z-moments are echo 1's scaled by the ratio of their dephasing times (e.g. 76/44 ms),
			# zero at the same frames; warn only if they are not a constant multiple
			zmCheck=$(awk -v a="$firstZm" -v b="$zmCsv" 'BEGIN {
				n = split(a, x, ","); m = split(b, y, ",")
				if (n != m) { print "a different number of frames"; exit }
				r = ""
				for (i = 1; i <= n; i++) {
					if ((x[i] + 0 == 0) != (y[i] + 0 == 0)) { print "zero at different frames"; exit }
					if (x[i] + 0 != 0) {
						q = y[i] / x[i]
						if (r == "") r = q
						else if (q / r > 1.01 || q / r < 0.99) { print "not a constant multiple"; exit }
					}
				}
				print "ok" }')
			[[ $zmCheck == ok ]] || warnMsg+="  - warning: ZMoment values for echo $echoLabel vs echo 1: $zmCheck (expected echo 1's scaled by the dephasing-time ratio)${nl}"
		fi

		block=$(printf '    {\n      "echo": %d,\n      "series": "%s",\n      "EchoTime_ms": %s,%s\n      "nFrames": %d,\n      "MaxDTE_ms": %s,\n      "StepDTE_ms": %s,\n      "dTE_ms": [%s],\n      "ZMoment": [%s],\n      "file": "%s"\n    }' \
			"$echoLabel" "$series" "$te_ms" "$sthLine" "$nLines" "$maxMs" "$stepMs" "$dteMs" "$zmCsv" "$proc")
		blocks+="${blocks:+,$nl}$block"
	done

	# actual field strength from the proton imaging frequency (MagneticFieldStrength is the nominal
	# value, e.g. 3 for a 2.89 T Prisma); used by aseFit.sh for t_crit and the OEF scale
	b0T=$(awk -v f="$freq" -v g="$gammaBarMHz" 'BEGIN { printf "%.6f", f / g }')
	{
		printf '{\n'
		printf '  "subject": "%s",\n' "$subjID"
		printf '  "created": "%s",\n' "$runTimestamp"
		printf '  "units": {"MagneticFieldStrength": "T", "ImagingFrequency": "MHz", "B0_T": "T", "EchoTime": "ms", "dTE": "ms", "MaxDTE": "ms", "StepDTE": "ms", "ZMoment": "Tms/m"},\n'
		printf '  "frame_order": "acquisition (header row) order; frame i of each echo file matches dTE_ms[i] and ZMoment[i]; not sorted",\n'
		printf '  "MagneticFieldStrength": %s,\n' "$B0"
		printf '  "ImagingFrequency": %s,\n' "$freq"
		printf '  "B0_T": %s,\n' "$b0T"
		printf '  "B0_source": "ImagingFrequency / %s MHz/T (proton gamma/2pi)",\n' "$gammaBarMHz"
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
# motionMode=echo1: estimate motion on echo 1, apply the same per-frame transforms to all echoes
stage_motion_echo1() {
	local ei series nii proc mc mc1
	series=${echoNums[0]}
	nii=$(subst_series "$rawNii" "$series")
	proc=$(subst_echo "$aseProcName" 1)
	mc1=$(subst_echo "$motionParamName" 1)
	echo "[$subjID] Motion correcting echo 1 -> $proc (transforms -> $motionMatName)"
	rm -f "$mc1" "$motionMatName"
	if ! 3dvolreg $volregOpts -base "$baseName" -dfile "$mc1" -1Dmatrix_save "$motionMatName" \
			-prefix "$proc" -overwrite "$nii"; then
		rm -f "$proc" "$mc1" "$motionMatName"
		procErr+="  - 3dvolreg failed for $nii${nl}"; return 1
	fi
	if [[ ! -e $proc || ! -s $motionMatName ]]; then
		procErr+="  - expected output not found after 3dvolreg: $proc / $motionMatName${nl}"; return 1
	fi
	for ((ei = 1; ei < nEchos; ei++)); do
		series=${echoNums[ei]}
		nii=$(subst_series "$rawNii" "$series")
		proc=$(subst_echo "$aseProcName" "$((ei + 1))")
		mc=$(subst_echo "$motionParamName" "$((ei + 1))")
		if [[ $(3dinfo -nv "$nii" 2>/dev/null) != $(3dinfo -nv "$(subst_series "$rawNii" "${echoNums[0]}")" 2>/dev/null) ]]; then
			procErr+="  - echo $((ei + 1)) has a different number of frames than echo 1; cannot reuse its transforms${nl}"; return 1
		fi
		echo "[$subjID] Applying echo-1 transforms to echo $((ei + 1)) -> $proc ($motionApplyInterp)"
		if ! 3dAllineate -1Dmatrix_apply "$motionMatName" -source "$nii" -master "$baseName" \
				-final "$motionApplyInterp" -prefix "$proc" -overwrite; then
			rm -f "$proc"
			procErr+="  - 3dAllineate (apply echo-1 motion) failed for $nii${nl}"; return 1
		fi
		cp -f "$mc1" "$mc"                     # same motion as echo 1, by construction
	done
	return 0
}

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
	case $motionMode in
		echo1)   stage_motion_echo1; return ;;
		perEcho) ;;
		*) procErr+="  - unknown motionMode '$motionMode' (echo1 | perEcho)${nl}"; return 1 ;;
	esac
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
	if [[ $anatSource == fsrecon ]]; then stage_fsrecon; return; fi
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

# anatSource=fsrecon: FreeSurfer volumes -> NIfTI (no resampling: they share FreeSurfer's 256^3 grid),
# then the registration target and the xvent mask on that grid (docs/anat-route-report.md)
stage_fsrecon() {
	local j name src expr="" lab
	for ((j = 0; j < nXfer; j++)); do
		((xferDo[j] == 1)) || continue
		name=${xferList[j]}; src="$fsBase/mri/${xferSrcList[j]}"
		echo "[$subjID] Converting $src -> $name"
		rm -f "$name"
		if ! mri_convert "$src" "$subjBase/$name" > /dev/null; then
			rm -f "$name"
			procErr+="  - mri_convert failed for $src${nl}"; return 1
		fi
		if [[ ! -e $name ]]; then
			procErr+="  - expected file not found after mri_convert: $subjBase/$name${nl}"; return 1
		fi
	done
	if ((brainDo == 1)); then
		echo "[$subjID] Building registration target $regBaseName: $fsT1 x step($fsAseg), 3dUnifize"
		rm -f "$regBaseName" .anat_brain_tmp.nii.gz
		if ! 3dcalc -a "$fsT1" -b "$fsAseg" -expr 'a*step(b)' -datum float -prefix .anat_brain_tmp.nii.gz -overwrite || \
				! 3dUnifize -prefix "$regBaseName" -overwrite .anat_brain_tmp.nii.gz; then
			rm -f "$regBaseName" .anat_brain_tmp.nii.gz
			procErr+="  - could not build $regBaseName${nl}"; return 1
		fi
		rm -f .anat_brain_tmp.nii.gz
	fi
	if ((xventDo == 1)); then
		for lab in ${xventExcludeLabels//,/ }; do expr+="${expr:++}equals(a,$lab)"; done
		echo "[$subjID] Building $xventName: $fsAparc > 0, minus labels $xventExcludeLabels"
		if ! 3dcalc -a "$fsAparc" -expr "step(a)-step($expr)" -datum byte -prefix "$xventName" -overwrite; then
			rm -f "$xventName"
			procErr+="  - could not build $xventName${nl}"; return 1
		fi
	fi
	return 0
}

stage_register() {
	((nSeg > 0)) || return 0
	if ((doRegister == 1)); then
		ensure_base || return 1
		echo "[$subjID] Registering $baseName to $regBaseName ($registerDof parameters)"
		rm -f "$tformName" "$regStampName"
		if ! 3dAllineate -prefix "$baseAlName" -base "$regBaseName" -master "$regBaseName" -source "$baseName" \
				$allineateOpts -warp "$regWarp" -1Dmatrix_save "$tformName" -overwrite; then
			rm -f "$tformName" "$baseAlName"
			procErr+="  - 3dAllineate registration failed${nl}"; return 1
		fi
		if [[ ! -s $tformName ]]; then
			procErr+="  - transform file not produced: $tformName${nl}"; return 1
		fi
		printf '%s' "$regSettings" > "$regStampName"
	elif [[ ! -e $regStampName && $anatSource == afni ]]; then
		printf '%s' "$regSettings" > "$regStampName"      # existing legacy registration: record its settings
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

# tissue fraction = 1 - (background + CSF + ventricle posteriors), on the anatomical grid
stage_tissueFrac() {
	local src nv letters="abcdefghijklmnopqrstuvwxyz" args=() expr="" i v vols
	((nSeg > 0 && tfDo == 1)) || return 0
	src="$transferSrcDir/$tissuePosteriorName"          # relative to subjDir (AFNI dislikes spaces in paths)
	nv=$(cd "$subjDir" && 3dinfo -nv "$src" 2>/dev/null)
	if [[ $nv != "$posteriorNVols" ]]; then
		procErr+="  - $tissuePosteriorName has ${nv:-?} sub-bricks, expected $posteriorNVols (SynthSeg label set)${nl}"; return 1
	fi
	IFS=',' read -ra vols <<< "$posteriorNonTissueVols"
	for ((i = 0; i < ${#vols[@]}; i++)); do
		v=${vols[i]}
		if [[ ! $v =~ ^[0-9]+$ ]] || ((v >= nv)) || ((i >= 26)); then
			procErr+="  - invalid posteriorNonTissueVols entry '$v'${nl}"; return 1
		fi
		args+=("-${letters:i:1}" "${src}[$v]")
		expr+="${expr:++}${letters:i:1}"
	done
	echo "[$subjID] Building $tissueFracName from $tissuePosteriorName (1 - sub-bricks $posteriorNonTissueVols)"
	if ! ( cd "$subjDir" && 3dcalc "${args[@]}" -expr "max(0,min(1,1-($expr)))" -datum float \
			-prefix "$aseSubjDir/$tissueFracName" -overwrite ); then
		rm -f "$subjBase/$tissueFracName"
		procErr+="  - 3dcalc failed building $tissueFracName${nl}"; return 1
	fi
	if [[ ! -e $subjBase/$tissueFracName ]]; then
		procErr+="  - expected output not found after 3dcalc: $subjBase/$tissueFracName${nl}"; return 1
	fi
	return 0
}

# quality masks -> ASE grid with linear interpolation, so the result is a fraction (0-1)
stage_qualityMasks() {
	local k qm base out
	((nSeg > 0 && nQm > 0)) || return 0
	for ((k = 0; k < nQm; k++)); do
		((qmDo[k] == 1)) || continue
		ensure_base || return 1
		qm=${qmArr[k]}; base=$(strip_nii "$qm"); out=${base}${segSuffix}.nii.gz
		echo "[$subjID] Resampling $qm -> $out on the ASE grid (linear; matrix: $tformApplyName)"
		if ! 3dAllineate -1Dmatrix_apply "$tformApplyName" -source "$qm" -master "$baseName" -prefix "$out" -final linear -overwrite; then
			rm -f "$out"
			procErr+="  - 3dAllineate (quality mask) failed for $qm${nl}"; return 1
		fi
		if [[ ! -e $out ]]; then
			procErr+="  - expected output not found after 3dAllineate: $out${nl}"; return 1
		fi
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
	stage_tissueFrac || return 1
	stage_qualityMasks || return 1
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
	fsBase=${fsSubjDir//%SUBJ%/$subjID}
	[[ $fsBase == /* ]] || fsBase="$subjDir/$fsBase"
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