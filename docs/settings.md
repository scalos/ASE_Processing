# Settings reference

Every setting lives in `BashFiles/ase_config.json`. Precedence: script defaults < config < command-line flags.

- **Pipeline keys** are read by `asePreProcess.sh` (override with `--set key=value`) and `aseFit.sh` (dedicated flags, e.g. `--segNames`).
- **Fit keys** are read by `runASEFromConfig.m`. Override with `aseFit.sh --fit key=value`; values are JSON.
- **The `quality` block** merges key by key, e.g. `--fit 'quality={"enabled":true}'` keeps every other quality key. Unknown quality keys raise an error.
- **Every run records the full resolved configuration** in `fits/<timestamp>/config_resolved.json`.

"Changed" marks defaults that differ from the pre-November-2026 pipeline; the reasons are in [decisions.md](decisions.md).

## Pipeline: data layout and preprocessing

| key | default | meaning |
|---|---|---|
| `allSubjDir` | `/autofs/…/yen/ase` | folder holding all subject folders and `dataList` |
| `dataList` | `datalist_ase` | subject list: column 1 is the subject ID, column `echo1Col` the echo-1 series |
| `echo1Col`, `echoSpacing`, `nEchos` | 28, 2, 2 | echo-k series = echo-1 series + (k−1)·echoSpacing |
| `aseSubjDir` | `session_ase` | per-subject working folder |
| `rawDataDir`, `dcmIndexFile` | `raw_data/dcm`, `dcm.index.dat` | DICOM location and index |
| `scannerCode`, `commentScript` | 3, `Siemens_ase_XA30.sh` | decodes dTE / z-moment from the DICOM image comments |
| `dcm2niixMinFiles` | 500 | minimum DICOM count for a series to convert |
| `rawNii`, `pattern_dicomHdr`, `pattern_dteHdr` | `ase_%SERIES%.nii`, … | names of the converted series and their headers |
| `hdrFileName` | `ase_hdr.json` | combined header (TEs, dTEs, z-moments, `SliceThickness_mm`) |
| `aseProcName`, `motionParamName` | `ase_e%ECHO%_proc.nii.gz`, `ase_e%ECHO%_mc.1D` | motion-corrected data and motion parameters |
| `volregOpts` | `-Fourier -twopass -twodup` | 3dvolreg options |
| **`motionMode`** | **`echo1`** (changed) | `echo1`: 3dvolreg on echo 1, same transforms applied to all echoes. `perEcho`: register each echo separately (old behaviour). |
| `motionMatName`, `motionApplyInterp` | `ase_e1_mc_mat.aff12.1D`, `wsinc5` | saved echo-1 transforms; interpolation used when applying them to echo 2 |
| `baseName`, `baseFrame` | `ase_base.nii.gz`, 0 | motion and registration reference (echo-1 frame 0, the spin echo) |
| `baseAlName`, `regBaseName`, `allineateOpts` | …, `anat_fs_brain_al.nii.gz`, `-cost lpa+ …` | ASE → anatomy registration |
| `tformName`, `tformApplyName` | `ase_base_al_mat.aff12.1D`, `…_inv_mat…` | forward and inverse transforms |
| `transferSrcDir`, `transferNames` | `session_afni_nortc`, brain, mask, xvent, bstem | anatomy copied from AFNI `+orig` to NIfTI. Entries that a later stage reads (`regBaseName`, `segNames`, `qualityMaskNames`) are required. Other entries, e.g. bstem, are optional: a missing source is a warning, and the file is transferred on a later run once the source exists. |
| `segNames`, `segSuffix`, `segLUTs` | aparc+aseg, wmparc; `_ase`; `GM_WM_long.txt` | segmentations moved to ASE space (NN) for the ROI table; empty `segNames` means no ROI stats. In `aseFit.sh`, a missing `<seg>_ase` file is left out with a warning (its ROIs are absent from the table); the subject is skipped only when none is found. |
| `qualityMaskNames` | `anat_fs_xvent_mask_al.nii.gz` | tissue-fraction masks resampled to the ASE grid (linear); CSF = 1 − fraction |
| `tissuePosteriorName` | `""` (off) | optional SynthSeg posteriors (e.g. `anat_fs_al_posteriors.nii.gz`) for an alternative tissue fraction |
| `tissueFracName`, `posteriorNonTissueVols`, `posteriorNVols` | `anat_fs_tissueFrac_al.nii.gz`, `0,3,4,11,12,16,21,22`, 33 | SynthSeg tissue fraction = 1 − (background + CSF + ventricles). To use it, also set `qualityMaskNames` to `tissueFracName`. |
| `fitsSubdir`, `outPref` | `fits`, `ase_fit_results` | output folder (`<subject>/<aseSubjDir>/<fitsSubdir>/<timestamp>/`) and file prefix. `aseFit.sh --tag <label>` names the run folder `<timestamp>_<label>` so later steps can select it (`--fitTag`). |
| `matlabBin`, `matdir` | `matlab`, repo | MATLAB location |
| **`B0`** | null = header `B0_T` (changed) | field strength (T) for t_crit and the OEF scale. By default the **actual** field from each subject's header, `B0_T` = ImagingFrequency / 42.577478518 MHz/T, written by `asePreProcess.sh`. It used to be the nominal `MagneticFieldStrength` (3). Set a number (or `aseFit.sh --B0`) to override. Headers without `B0_T` are rejected: rerun `asePreProcess.sh --regenHdr`. |

**Re-running preprocessing stages:** `--regenNifti`, `--regenComments`, `--regenHdr`, `--regenMotion`, `--regenTransfer`, `--regenRegister`, `--regenTransform`. Each one also forces the stages downstream of it. Without them, existing outputs are kept.

## Fit: model and frames

| key | default | meaning |
|---|---|---|
| `Hct`, `Hct_correction` | 0.42, 0.85 | large-vessel Hct × small/large-vessel ratio |
| `dChi_0` | 1.8e-7 | Δχ between fully oxygenated and fully deoxygenated blood (SI, per unit Hct) |
| `StO2_max` | 0.6 | sets t_c = 1/δω(1 − StO2_max), and t_crit = 1.5·t_c. Used **only** for the regime split. **Flagged for revisit:** 0.6 is a typical, not a maximal, venous saturation. |
| **`tCrit_buffer`** | **`[3, 12]`** (changed from 3) | ms below / above t_crit that are excluded. A scalar means symmetric. |
| `tauRange` | `[]` (all) | restrict τ to [min max] ms |
| `tauUsePts` | `[]` | 98 × nEchos 0/1 mask of frames to use |
| **`zshimZeroFrameInFit`** | **false** (changed) | include the z-shim block's zero-moment frame (raw frame 96) in the regime fits |
| `averageReps` | false | average the repetitions before fitting (default: fit them stacked) |
| `dBCorr` | false | through-slice z-shim correction (settings in `quality.zshim_*`) |
| `physBd_R2`, `physBd_R2P`, `physBd_vCBV`, `physBd_OEF` | [0 Inf], [0 Inf], [0 1], [0 1] | plausible ranges: output maps are clipped to them, and the quality bit 0 test uses them |
| `physBdClipVal` | null | null = clip to the bound; a number = set out-of-range values to it |
| `normErr` | true | normalize the error map by mean signal |
| `minOverlapFrac` | 0.25 | minimum fraction of an ASE voxel inside a region for it to count as a member |
| `writeOpts`, `overwriteOpts` | maps, coeffs, stats | what to write and what may be overwritten |
| `writeRawMaps` | false | also write unclipped `_raw` maps |
| `warnSegIssues`, `verbose` | true, true | |

## Fit: `quality` block (`Aux_Mat_Functions/qualityDefaults.m`)

**Master switches.** All off by default, so default outputs never depend on them:

| key | default | meaning |
|---|---|---|
| `enabled` | false | build the voxel quality mask and apply it to ROI stats (turns on whatever diagnostics the enabled criteria need) |
| `repNoise` | false | write σ / SNR maps |
| `zshimEstimate` | false | write Gz maps without correcting |
| `fitExtensions` | false | write vCBV_quad, SEs, z, kZ and per-echo R2′ |
| `noiseWeighting` | false | weight ROI membership by echo-1 SNR² |
| `minVoxKept` | 5 | regions with fewer kept effective voxels get NaN |

**z-shim settings.** Used when `dBCorr` or a z-shim criterion is on:

| key | default | meaning |
|---|---|---|
| `zshim_GzRange`, `zshim_GzStep` | [−150 150], 1 | Gz search grid (µT/m) |
| `zshim_weightF2` | true | weight corrected frames by F² (false = old behaviour: unweighted, NaN below `zshim_Fmin`) |
| `zshim_Fcut` | 0.1 | frames with F below this are dropped |
| `zshim_Fmin` | 0.5 | "heavily attenuated" threshold for bit 8 |
| `zshim_jointGz` | false | one Gz fitted to both echoes. Available but off: it is biased where the echoes disagree. |

**Criteria** (bit numbers as in [pipeline.md](pipeline.md#4-corrections-and-diagnostics-all-optional)):

| key(s) | default | notes |
|---|---|---|
| `useBounds` | true | bit 0 |
| `useQuadSign`, `quadSign_kZmin` | true, 2 | bit 1: k/SE(k) > 2 |
| `useVcbvZ`, `vcbvZ_max` | true, 3 | bit 2: **placeholder; recalibrate** |
| `useSNR`, `snr_min` | true, 10 | bit 3 |
| `useCSF`, `csf_max` | true, 0.2 | bit 4: needs `qualityMaskFiles` |
| `useVein`, `vein_mode`, `vein_relK`, `vein_refTissueMin`, `vein_vCBVmax` | true, `relative`, 4, 0.9, 0.15 | bit 5. Relative mode: threshold = 4 × the median vCBV of tissue voxels (CSF < 10%, plausible fit). `vein_vCBVmax` is used only in `absolute` mode. |
| `useZshimEdge` | true | bit 6 |
| `useZshimEcho`, `zshimEcho_maxDiff` | true, 20 | bit 7: **placeholder; recalibrate** |
| `useZshimFrac`, `zshimFrac_max` | true, 0.3 | bit 8 |
| `useR2PEcho`, `r2pEcho_pMax` | **false**, 0.01 | bit 9: off, because it flags a systematic model effect rather than bad voxels |

## Group statistics (`aseGroupStats.sh` → `ASE_groupROIStats.m`)

Command-line flags only (no config section).

| flag | default | meaning |
|---|---|---|
| `--covTable` | `<allSubjDir>/cov_table_ase` | **Defines the subjects.** Header line (may start with `#`), one row per subject, first column = subject ID. `#` rows are excluded; `NA` = missing. The data list is not used. |
| `--select` (repeatable) | `group=ctrl,mdt` | row filters, all must hold: `col=a,b`, `col!=a`, `col>x`, `col>=x`, `col<x`, `col<=x`. `--selectAll` uses every row. |
| `--covariates` | `age` | covariate columns (`""` = none). Numeric ones are centred on the cohort mean. |
| `--categorical` | none | numeric columns to treat as categorical (e.g. `sex`); text columns always are (reference = first level in sorted order) |
| `--stats` | `OEF,R2P,vCBV` | ROI-table columns to model |
| `--minVox` | 0 | use a subject's ROI value only if its `nVox` ≥ this |
| `--minN` | 5 | minimum subjects for a model (and ≥ coefficients + 2) |
| `--fitTag` | none = latest run | each subject's latest run tagged `<label>`; no fallback to other runs |
| `--settingsCheck` | `error` | `error` / `warn` / `off` when the chosen runs' `config_resolved.json` "fit" sections differ (output-only keys ignored) |
| `--outPath` | `<allSubjDir>/group_stats` | results go to `<outPath>/<timestamp>[_<fitTag>]/` |
| `--paint` + `--templateDir`, `--templateNames`, `--LUTs`, … | off | paint the result with `paintStats.sh --table`; columns default to every `<stat>_adjMean` and `<stat>_<term>_t` (`--paintStats` to choose) |

Model per ROI and stat: `value ~ 1 + covariates` (`fitlm`). A subject's value is used when finite and above `minVox`; subjects with a missing covariate are left out. q = Benjamini–Hochberg FDR across ROIs, per stat and term.

## Painting (`paintStats.sh` → `ASE_paintROIStats.m`)

| flag | default | meaning |
|---|---|---|
| `--templateDir`, `--templateNames` | `/autofs/…/yen/templates`, required | template segmentation(s) |
| `--LUTs` | `GM_WM_long.txt` | **painting LUT(s).** Regions are placed in LUT order; a region that overlaps an earlier one is dropped whole (no trimming) and listed in the report. List only the regions to paint. |
| `--targetStats` | `OEF,R2P,vCBV` | table columns to paint |
| `--fitTag`, `--allRuns` | latest run | per-subject mode: which run(s) to paint |
| `--table`, `--outDir` | — | table mode: paint one table (e.g. `group_ROI_Stats.csv`) |
| `--fillValue` | **NaN** (changed) | value outside painted regions and for regions without a value (was 0, which looked like a real value) |
| `--omitCodes`, `--omitNames` | none | LUT entries not to paint |
