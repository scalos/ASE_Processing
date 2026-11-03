#!/bin/bash
# ---- defaults (accessible from command line) ----
allSubjDir="/autofs/space/guasha_004/users/crbc/yen/ase"
dlist=datalist_ase
echo1Col=28
matdir="/autofs/space/guasha_004/users/crbc/yen/ase/sc_analysis_code"
B0=3
regenTransform=0
interpMethod=linear
tformName=ase_base_al_mat.aff12.1D
# segmentation-mode inputs (used when -roiMode seg)
segNamesStr="anat_fs_aparc_aseg_al.nii.gz,anat_fs_wmparc_al.nii.gz,anat_fs_bstem_al.nii.gz"
segLUTsStr="/autofs/space/guasha_004/users/crbc/yen/ase/GM_WM_long.txt"
masterSeg=""   # grid used as -master for 3dAllineate; defaults to first entry of segNames if blank

# ROI mode: files | seg | none
roiMode=seg
roiDirName=allROIs   # used only when -roiMode files

# repetition discovery. %ECHO% is substituted with the echo number (x3) per subject.
repPattern="ase_%ECHO%_reg_r?.nii.gz"  # NOTE: only supports r1-r9 due to the single-character wildcard

# output layout
fitsSubdir=fits
writeOptsStr="maps,coeffs,stats"
overwriteOptsStr="maps,coeffs,stats"
tauRangeStr="0,43"   # leave empty ("") to omit tauRange entirely
tCritBuffer=3        # ms on either side of critical time for ASE_Fitting_V3

runTimestamp=$(date +%Y%m%d_%H%M%S)

# maintenance actions (run instead of the main pipeline, then exit)
clearAllFits=0
keepLatestN=0        # 0 = disabled
rmLatestN=0          # 0 = disabled
# todo: add purgeEmpty

# ---- defaults (rarely changed, still CLI-accessible) ----
aseSubjDir=session_ase   # subject-relative path to ase data (no leading/trailing slash)
hdrPattern=dcm_hdr_%ECHO%
matlabBin=matlab
omitFailLog=0
retryFromFile=""     # failure file listing completed subjects; empty = disabled
retryMostRecent=0

tok='%ECHO%'
nl=$'\n'

usage() {
cat <<EOF
Usage: $0 [options]

Environment/Pipeline (Dependent on user / pipeline structure):
  -allSubjDir <path>         Top-level directory containing subject folders (default: $allSubjDir)
  -dlist <file>              Subject list filename, located in allSubjDir (default: $dlist)
  -matdir <path>             Path to MATLAB analysis code (default: $matdir)
  -aseSubjDir <dir>          Subject-relative path to ase data (default: $aseSubjDir)
  -hdrPattern <str>          Pattern (with %ECHO% token) for header files (default: $hdrPattern)
  -matlabBin <path>          Path to MATLAB binary (default: $matlabBin)
  -echo1Col <value>          Column in header file corresponding to echo 1 (default: $echo1Col)

ASE Transforms:
  -repPattern <glob>         Glob (with %ECHO% token) to discover ASE repetition files (default: $repPattern)
  -interp <method>           3dAllineate interpolation: NN, linear, cubic, quintic, wsinc5 (default: $interpMethod)
  -tformName <file>          File name for transform from subj space -> standard space (default: $tformName)
  -regenTransform            Force regeneration of anat-space ASE transforms (default: off)

ROI source:
  -roiMode <files|seg|none>  How to derive ROIs (default: $roiMode)
  -roiDirName <dir>          [files mode] per-subject dir of individual ROI masks (default: $roiDirName)
  -segNames <a,b,...>        [seg mode] segmentation volume filenames, relative to aseSubjDir (default: $segNamesStr)
  -segLUTs <a,b,...>         [seg mode] LUT file path(s) matching segNames (default: $segLUTsStr)
  -masterSeg <file>          Grid used as -master for the ASE transform step (default: first entry of -segNames)

Fitting options:
  -tauRange <lo,hi>          Tau range in ms passed to ASE_Fitting_V3 (default: [$tauRangeStr]; empty = omit)
  -B0 <value>                Field strength in Tesla passed to ASE_Fitting_V3 (default: $B0)
  -tCritBuffer <value>       Time buffer in ms on either side of critical time (default: $tCritBuffer)

Fitting outputs:
  -fitsSubdir <name>         Subdirectory (per-subject) that holds timestamped run folders (default: $fitsSubdir)
  -writeOpts <a,b,...>       maps,coeffs,stats subset to write (default: $writeOptsStr)
  -overwriteOpts <a,b,...>   maps,coeffs,stats subset allowed to overwrite (default: $overwriteOptsStr)

Maintenance (run standalone, then exit):
  -clearAllFits              Delete every subject's entire fitsSubdir, then exit
  -keepLatest <N>            Keep only the N most recent run folders per subject, delete the rest, then exit
  -rmLatest <N>              Remove the N most recent run folders per subject, then exit

Other:
  -help                      Show this message
  -omitFailLog               Do not write a failure log file (default: write log)
  -retryFromFile <file>      Failure file listing subjects already completed (one per line); empty = disabled
  -retryMostRecent           Retry using the most recent failure file in allSubjDir (default: off)
EOF
}

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

# ---- parse command-line flags ----
while [[ $# -gt 0 ]]; do
	case "$1" in
		-allSubjDir)     need_arg "$@"; allSubjDir=$2;      shift 2 ;;
		-dlist)          need_arg "$@"; dlist=$2;           shift 2 ;;
		-matdir)         need_arg "$@"; matdir=$2;          shift 2 ;;
		-B0)             need_arg "$@"; B0=$2;              shift 2 ;;
		-interp)         need_arg "$@"; interpMethod=$2;    shift 2 ;;
		-tformName)      need_arg "$@"; tformName=$2;       shift 2 ;;
		-segNames)       need_arg "$@"; segNamesStr=$2;     shift 2 ;;
		-segLUTs)        need_arg "$@"; segLUTsStr=$2;      shift 2 ;;
		-masterSeg)      need_arg "$@"; masterSeg=$2;       shift 2 ;;
		-roiMode)        need_arg "$@"; roiMode=$2;         shift 2 ;;
		-roiDirName)     need_arg "$@"; roiDirName=$2;      shift 2 ;;
		-repPattern)     need_arg "$@"; repPattern=$2;      shift 2 ;;
		-fitsSubdir)     need_arg "$@"; fitsSubdir=$2;      shift 2 ;;
		-writeOpts)      need_arg "$@"; writeOptsStr=$2;    shift 2 ;;
		-overwriteOpts)  need_arg "$@"; overwriteOptsStr=$2; shift 2 ;;
		-tauRange)       need_arg "$@"; tauRangeStr=$2;     shift 2 ;;
		-aseSubjDir)     need_arg "$@"; aseSubjDir=$2;      shift 2 ;;
		-hdrPattern)     need_arg "$@"; hdrPattern=$2;      shift 2 ;;
		-keepLatest)     need_arg "$@"; keepLatestN=$2;     shift 2 ;;
		-rmLatest)       need_arg "$@"; rmLatestN=$2;       shift 2 ;;
		-matlabBin)      need_arg "$@"; matlabBin=$2;       shift 2 ;;
		-echo1Col)       need_arg "$@"; echo1Col=$2;        shift 2 ;;
		-tCritBuffer)    need_arg "$@"; tCritBuffer=$2;     shift 2 ;;
		-retryFromFile)  need_arg "$@"; retryFromFile=$2;   shift 2 ;;
		-regenTransform) regenTransform=1; shift ;;
		-clearAllFits)   clearAllFits=1;   shift ;;
		-omitFailLog)    omitFailLog=1;    shift ;;
		-retryMostRecent) retryMostRecent=1; shift ;;
		-h|-help)        usage; exit 0 ;;
		*)
			echo "Unknown option: $1" >&2
			echo "Usage: $0 [options]  (use -help for details)" >&2
			exit 1
			;;
	esac
done

# ---- cleanup runs on any exit (normal, error, Ctrl-C) ----
cleanup() {
	if [[ -n $mfile && -e $mfile ]]; then
		rm -f "$mfile"
		echo "Removed temporary matlab script: $mfile"
	fi
	if [[ -n $runOutDir && -d $runOutDir ]]; then
		# rmdir only succeeds on an empty directory
		rmdir "$runOutDir" 2>/dev/null && echo "Removed empty runOutDir: $runOutDir"
	fi
	if [[ -n $failFile && -e $failFile && $omitFailLog -eq 1 ]]; then
		rm -f "$failFile"
		echo "Removed failure file: $failFile"
	fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

####################### MAINTENANCE ###############################
cd "$allSubjDir" || { echo "Error: cannot cd to $allSubjDir" >&2; exit 1; }
mapfile -t subjlist < <(grep -v '#' "$allSubjDir/$dlist" | awk '{print $1}')

# ---- maintenance actions: run and exit before touching matlab/transforms ----
if ((clearAllFits == 1)); then
	echo "Clearing all fit output for ${#subjlist[@]} subject(s)..."
	for x in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$x/$aseSubjDir/$fitsSubdir"
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
	for x in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$x/$aseSubjDir/$fitsSubdir"
		[[ -d $fitsPth ]] || continue
		mapfile -t runDirs < <(ls -1 "$fitsPth")
		nRuns=${#runDirs[@]}
		((nRuns <= keepLatestN)) && continue
		nToDelete=$((nRuns - keepLatestN))
		for ((i = 0; i < nToDelete; i++)); do   # bash arrays are 0-based
			rm -rf "$fitsPth/${runDirs[i]}"
			echo "  Removed: $fitsPth/${runDirs[i]}"
		done
	done
	echo "Done."
	exit 0
fi

if ((rmLatestN > 0)); then
	echo "Pruning the most recent $rmLatestN run(s) per subject..."
	for x in "${subjlist[@]}"; do
		fitsPth="$allSubjDir/$x/$aseSubjDir/$fitsSubdir"
		[[ -d $fitsPth ]] || continue
		mapfile -t runDirs < <(ls -1 "$fitsPth")
		nRuns=${#runDirs[@]}
		((nRuns <= rmLatestN)) && continue
		nToDelete=$((nRuns - rmLatestN))
		for ((i = nRuns - 1; i >= nToDelete; i--)); do
			rm -rf "$fitsPth/${runDirs[i]}"
			echo "  Removed: $fitsPth/${runDirs[i]}"
		done
	done
	echo "Done."
	exit 0
fi

##########################################################################
# ---- validate roiMode up front ----
if [[ $roiMode != files && $roiMode != seg && $roiMode != none ]]; then
	echo "Error: -roiMode must be one of: files, seg, none (got: $roiMode)" >&2
	exit 1
fi

# ---- split comma-separated lists into arrays ----
IFS=',' read -ra segNames <<< "$segNamesStr"
IFS=',' read -ra segLUTs  <<< "$segLUTsStr"
nSegNames=${#segNames[@]}
nSegLUTs=${#segLUTs[@]}

if [[ $roiMode == seg ]]; then
	if ((nSegNames == 0)); then
		echo "Error: -roiMode seg requires at least one entry in -segNames" >&2
		exit 1
	fi
	if ((nSegLUTs == 0)); then
		echo "Error: -roiMode seg requires at least one entry in -segLUTs" >&2
		exit 1
	fi
fi

if [[ -z $masterSeg ]]; then
	if ((nSegNames > 0)); then
		masterSeg=${segNames[0]}
	else
		masterSeg=anat_fs_aparc_aseg_al.nii.gz
	fi
fi

segLUTsCell=""
[[ $roiMode == seg ]] && segLUTsCell=$(make_cell "${segLUTs[@]}")

mfile="$allSubjDir/run_ase_fit_batch.m"
rm -f "$mfile"
echo "addpath(genpath('${matdir}'));" >> "$mfile"
echo "dbstop if error;" >> "$mfile"

# writeOpts/overwriteOpts are plain keyword lists (no spaces expected)
writeOptsCell="{'$(echo "$writeOptsStr" | sed "s/,/','/g")'}"
overwriteOptsCell="{'$(echo "$overwriteOptsStr" | sed "s/,/','/g")'}"

tauRangeArg=""
[[ -n $tauRangeStr ]] && tauRangeArg=",'tauRange',[${tauRangeStr}]"

tCritBufferArg=""
[[ -n $tCritBuffer ]] && tCritBufferArg=",'tCrit_buffer',[${tCritBuffer}]"

nsubj=${#subjlist[@]}

# running totals for estimating avg time per transform
totTime=0
compCount=0
avgTime=0

failLog=""
nfail=0
nqueued=0

# ---- retry handling ----
retryFile=""
if ((retryMostRecent == 1)); then
	retryFile=$(ls -1t "$allSubjDir"/fitASE_failures_*.txt 2>/dev/null | head -n 1)
fi
if [[ -n $retryFromFile ]]; then
	if [[ -n $retryFile && $retryFile != "$retryFromFile" ]]; then
		echo "Warning: -retryMostRecent and -retryFromFile refer to different files; using -retryFromFile ($retryFromFile)"
	fi
	retryFile=$retryFromFile
fi
if [[ -n $retryFile && ! -e $retryFile ]]; then
	echo "Error: retry requested but specified file not found: $retryFile" >&2
	exit 1
fi
completedSubjs=()
if [[ -n $retryFile ]]; then
	mapfile -t completedSubjs < <(awk '/^\+/ {print substr($1, 2)}' "$retryFile")
	printf 'Skipping completed subjects from file: %s\nCompleted subjects:\n' "$retryFile"
	for subj in "${completedSubjs[@]}"; do
		printf '  - %s\n' "$subj"
	done
fi
failFile="$allSubjDir/fitASE_failures_${runTimestamp}.txt"
echo "Run started: $(date)" > "$failFile"
echo "Completed subjects:" >> "$failFile"

echo "=================================================="
echo "Starting batch ASE fitting for $nsubj subject(s)..."
subj_i=0
for x in "${subjlist[@]}"; do
	((subj_i++))
	if ((${#completedSubjs[@]} > 0)); then
		if printf '%s\n' "${completedSubjs[@]}" | grep -qxF "$x"; then
			echo "Skipping $x (already completed)"
			echo "+$x" >> "$failFile"
			continue
		fi
	fi

	subjBase="$allSubjDir/$x/$aseSubjDir"

	okay=1
	missingReasons=""

	if [[ ! -d $subjBase ]]; then
		msg="Skipping $x:${nl}  - subject directory not found: $subjBase${nl}"
		printf '%s' "$msg"
		failLog+=$msg
		((nfail++))
		continue
	fi

	cd "$subjBase" || continue

	# csh arrays are 1-based, bash arrays are 0-based -> echo1Col-1
	datatmp=( $(grep "$x" "$allSubjDir/$dlist") )
	x3=${datatmp[echo1Col-1]//ase_/}   # first echo number

	hdrPth="$subjBase/${hdrPattern//"$tok"/$x3}"
	matPth="$subjBase/$tformName"

	if [[ ! -e $hdrPth ]]; then
		okay=0
		missingReasons+="  - header file not found: $hdrPth${nl}"
	fi
	if [[ ! -e $matPth ]]; then
		okay=0
		missingReasons+="  - transform file not found: $matPth${nl}"
	fi

	# ---- discover repetition files (we're cd'd into subjBase, so globs are relative) ----
	repGlob=${repPattern//"$tok"/$x3}
	mapfile -t repFiles < <(ls -v $repGlob 2>/dev/null)   # unquoted on purpose: glob must expand
	nReps=${#repFiles[@]}
	if ((nReps == 0)); then
		okay=0
		missingReasons+="  - no repetition files matched: $subjBase/$repGlob${nl}"
	fi

	# ---- resolve ROI inputs for this subject, per roiMode ----
	roiDirPth="$subjBase/$roiDirName"
	if [[ $roiMode == files ]]; then
		if [[ ! -d $roiDirPth ]]; then
			okay=0
			missingReasons+="  - ROI directory not found: $roiDirPth${nl}"
		fi
	elif [[ $roiMode == seg ]]; then
		for segName in "${segNames[@]}"; do
			if [[ ! -e "$subjBase/$segName" ]]; then
				okay=0
				missingReasons+="  - segmentation file not found: $subjBase/$segName${nl}"
			fi
		done
	fi

	if ((okay == 0)); then
		msg="Skipping $x:${nl}${missingReasons}"
		printf '%s' "$msg"
		failLog+=$msg
		((nfail++))
		continue
	fi

	masterPth="$subjBase/$masterSeg"
	if [[ ! -e $masterPth ]]; then
		msg="Skipping $x: masterSeg file not found ($masterPth)${nl}"
		printf '%s' "$msg"
		failLog+=$msg
		((nfail++))
		continue
	fi

	# ---- transform each repetition from function space into anatomy space ----
	repAnatFiles=()
	for ((idx = 0; idx < nReps; idx++)); do
		rep_i=$((idx + 1))
		repFileBase=${repFiles[idx]}
		repAnatFile=${repFileBase%.nii.gz}_anat.nii.gz
		repAnatPth="$subjBase/$repAnatFile"
		repAnatFiles+=("$repAnatFile")

		if [[ ! -e $repAnatPth ]] || ((regenTransform == 1)); then
			t0=$(date +%s)
			3dAllineate -1Dmatrix_apply "$tformName" -source "$repFileBase" -master "$masterSeg" \
				-prefix "$repAnatFile" -final "$interpMethod" -overwrite
			t1=$(date +%s)
			dt=$((t1 - t0))
			totTime=$((totTime + dt))
			((compCount++))
			avgTime=$((totTime / compCount))

			# approximate remaining-work estimate: assumes remaining subjects have the same rep count
			repsLeftThisSubj=$((nReps - rep_i))
			subjsLeft=$((nsubj - subj_i))
			transformsLeft=$((repsLeftThisSubj + subjsLeft * nReps))
			estSec=$((avgTime * transformsLeft))
			hh=$((estSec / 3600))
			mm=$(((estSec % 3600) / 60))
			ss=$((estSec % 60))
			printf 'Transforming rep %d/%d -> anat space for %s (subject %d/%d)\n' "$rep_i" "$nReps" "$x" "$subj_i" "$nsubj"
			printf 'Est. time remaining: %02d:%02d:%02d\n\n' "$hh" "$mm" "$ss"
		fi
	done

	anyMissing=0
	for repAnatFile in "${repAnatFiles[@]}"; do
		[[ -e "$subjBase/$repAnatFile" ]] || anyMissing=1
	done
	if ((anyMissing == 1)); then
		msg="Transform failed for $x: one or more anat-space repetition files not produced${nl}"
		printf 'Warning: %s' "$msg"
		failLog+=$msg
		((nfail++))
		continue
	fi

	# ---- build asePths cell array from transformed repetition files ----
	asePthList=()
	for repAnatFile in "${repAnatFiles[@]}"; do
		asePthList+=("$subjBase/$repAnatFile")
	done
	asePths=$(make_cell "${asePthList[@]}")

	# ---- build ROI-related matlab args for this subject ----
	roiArgStr=""
	if [[ $roiMode == files ]]; then
		mapfile -t roiFiles < <(ls "$roiDirPth")
		if ((${#roiFiles[@]} == 0)); then
			msg="Skipping $x: no files found in $roiDirPth${nl}"
			printf '%s' "$msg"
			failLog+=$msg
			((nfail++))
			continue
		fi
		roiPthList=()
		for r in "${roiFiles[@]}"; do
			roiPthList+=("$roiDirPth/$r")
		done
		roiArgStr=",'ROI_Files',$(make_cell "${roiPthList[@]}")"
	elif [[ $roiMode == seg ]]; then
		segPthList=()
		for segName in "${segNames[@]}"; do
			segPthList+=("$subjBase/$segName")
		done
		roiArgStr=",'segFiles',$(make_cell "${segPthList[@]}"),'segLUTs',${segLUTsCell}"
	fi

	# ---- output location: one timestamped run folder per batch invocation, per subject ----
	runOutDir="$subjBase/$fitsSubdir/$runTimestamp"
	mkdir -p "$runOutDir"
	outPref="$runOutDir/ase_fit_${x3}"

	printf '%s\n' "runWithProgress(@() ASE_Fitting_V3(${B0},'${hdrPth}',${asePths}${roiArgStr},'outPref','${outPref}','verbose',true,'writeOpts',${writeOptsCell},'overwriteOpts',${overwriteOptsCell}${tauRangeArg}${tCritBufferArg}), '${x}', ${subj_i}, ${nsubj});" >> "$mfile"
	# record completed subject so a failed run can resume (append mode 'at', not 'wt')
	printf '%s\n' "fID = fopen('${failFile}', 'at');" >> "$mfile"
	printf '%s\n' "if fID ~=-1; fprintf(fID, '+${x}\n');end" >> "$mfile"
	printf '%s\n' "if fID ~=-1; fclose(fID);end" >> "$mfile"
	((nqueued++))
done
cd "$allSubjDir" || exit 1

echo "exit;" >> "$mfile"

if ((nqueued > 0)); then
	echo "Opening MATLAB..."
	"$matlabBin" -batch "run('${mfile}')"
else
	echo "No subjects queued; skipping matlab."
fi

rm -f "$mfile"

##########################################################################
echo ""
echo "=================================================="
echo "Run complete: $nqueued subject(s) queued for fitting"
if ((nfail == 0)); then
	echo "No failures."
	rm -f "$failFile"
else
	echo "$nfail failure(s) encountered:"
	printf '%s' "$failLog"

	if ((omitFailLog == 0)); then
		{
			echo ""
			echo "$nfail failure(s)"
			echo ""
			printf '%s' "$failLog"
		} >> "$failFile"
		echo "Failure details written to: $failFile"
	else
		rm -f "$failFile"
		echo "Failure details omitted (omitFailLog=1)"
	fi
fi
echo "=================================================="
