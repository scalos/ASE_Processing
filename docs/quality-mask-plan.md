# Action plan for Claude Code

Implement voxel-quality exclusion before ROI analysis in scalos/ASE\_Processing, in the phases below, one commit per phase. The other tabs in this doc hold the reasoning; this tab is the work order.

## Ground rules

- Default behaviour must not change: every new feature sits behind a config flag in a new `"quality"` block of the `"fit"` section, off by default until Phase 9 validation is done.
- Bash stays basic (no `mapfile`, no associative arrays), matching the existing scripts. MATLAB follows the existing `arguments`-block style.
- Never apply MATLAB's built-in `sinc` (normalized, sin(πx)/(πx)). Write the unnormalized form explicitly: `s = sin(x)./x; s(x==0) = 1;`.
- Every new map gets a manifest entry under `paths` and is written with `writeMap`/`writeVec`.
- Test data lives in `testData/` (gitignored). For now there is one local test subject; the user may add more later.

## Guardrails (non-negotiable)

**Git.**

- Before any edit, create a new branch from main (e.g. `quality-mask`) and do all work there.
- Editing existing files on the branch is allowed.
- Never merge into main, open-and-merge a PR, rebase main, force-push, or delete branches. Stop and ask the user before anything that touches main.
- Commit at the end of each phase with a message naming the phase.

**Data.**

- Work only with the local test subject in `testData/` (gitignored). Do not access the remote cluster.
- Never overwrite or delete existing data files (headers, processed volumes, previous fit runs). Before running a script, check what it writes: `asePreProcess.sh --regen*` overwrites in place, and `aseFit.sh` writes fit folders and a failure log into `allSubjDir`.
- Before the first run, copy the test subject to a scratch location (e.g. `testData/scratch/<subj>`) and point `--allSubjDir` there, so the original stays pristine.
- Ask the user before any step that would modify files outside the repo or outside the scratch copy.

**Purpose and cleanup.**

- The goal is to evaluate whether the code changes work, not to produce new datasets. Treat every output as disposable diagnostics.
- All outputs (fit runs, intermediate maps, logs, plots) go into scratch directories under `testData/scratch/`. Nothing is written next to the original subject data.
- Claude Code is responsible for cleaning up after itself: once the diagnostics for a step have been reviewed and summarized for the user, delete that step's scratch outputs. Keep only the Phase 0 baseline until the work is finished, then delete it too.
- Report findings as written summaries (numbers, short tables, file paths that were checked), so nothing needs to be kept on disk.
- Before ending a session, confirm that `testData/scratch/` is empty except for the baseline, and that the original subject folder is unchanged.

## Phase 0 — Baseline

1. Run the current pipeline (`aseFit.sh --singleSubjID <id>`) on the local test subject (in its scratch copy) and keep the outputs as the baseline: `_ROI_Stats.csv` and the maps.
2. Record runtime and peak memory per subject so the later phases can be compared.

## Phase 1 — Header metadata

1. `BashFiles/asePreProcess.sh`, `build_header`: read `SliceThickness` (mm) from each echo's dcm2niix sidecar with `json_num`, and write it into each echo block of `ase_hdr.json` as `"SliceThickness_mm"`.
2. `Aux_Mat_Functions/readASEHeader.m`: expose `info.SliceThickness_mm` (per echo), `info.zShimInds` (the trailing `nZshimming` frame indices) and `info.zMoments` (nZshim × nEchos, converted to T·s/m).
3. Missing `SliceThickness_mm` is allowed: error only when z-shim correction is enabled.
4. Regenerate the header for the test subjects with `asePreProcess.sh --regenHdr`.

## Phase 2 — Repetition noise

1. New `Aux_Mat_Functions/estimateRepNoise.m`. Input: the stacked data for one echo plus `info.repInds`. Output: σ map (X × Y × Z) per echo.
2. Method: d\_i = S\_rep1,i − S\_rep2,i over the nSweep points; subtract each voxel's mean d (drift); σ² = mean(d²)/2. For nReps > 2, pool per-point sample variances.
3. In `ASE_fitting.m`, compute σ from the unaveraged data. When `averageReps` is true, load the stack separately for this step.
4. Outputs: `_sigma_e1`, `_sigma_e2`, plus an SNR map (baseline signal / σ).

## Phase 3 — Z-shim gradient estimate and correction

1. New `CoreMatFiles/fitZshimGradient.m`. Per echo, grid search for Gz (default −150 to +150 µT/m, 1 µT/m step):
   1. Template f\_j(M) = |sinc(γ·Δz·(Gz\_j·t\_deph − M)/2)| at the 8 shim moments.
   2. A\_j = Σ y·f\_j / Σ f\_j²; R\_j = Σ (y − A\_j·f\_j)², vectorized over all voxels.
   3. Best j per voxel, with parabolic refinement. Flag grid-edge results and peaks outside the sampled moment range.
   4. t\_deph of the shim block = its entry in `taus_full` (44 ms for echo 1, 76 ms for echo 2 in the test protocol).
2. Outputs: `_Gz_e1`, `_Gz_e2`, the fitted amplitude, and the model-free attenuation ratio (zero-moment shim frame / A).
3. Correction: F = |sinc(γ·Δz·Gz·t\_deph/2)| per voxel and frame, with t\_deph = `taus_full`. Divide `ase_sorted` by F. Set values with F < `quality.zshim_Fmin` (default 0.5) to NaN voxelwise; the existing `frameMask` is per frame only, so this is a new voxel × frame mask.
4. Replace the `dBCorr` error in `ASE_fitting.m` with this path. The z-shim frames themselves stay excluded from the regime fits.
5. Unit test: synthesize a voxel with known Gz and check that the grid search recovers it within 2 µT/m at realistic noise.

## Phase 4 — Fit extensions

1. `fitASERegimes`: when the short design has a τ² column, also return k = `BS(:,:,:,2)` and vCBV\_quad = −0.3·R2′² / k (NaN where k ≥ 0).
2. Propagate σ (log space, σ/S per frame) through the OLS designs to get standard errors of vCBV\_int and vCBV\_quad. Output a z-map of their difference.
3. Per-echo long fits (only when nEchos = 2): fit each echo's long frames separately for R2′\_e1 and R2′\_e2. Get RSS for the joint and split models (`fitMultiRegressor` already returns rms and nValid, so RSS = rms²·n). F = (RSS\_joint − RSS\_split) / (RSS\_split/(n−4)).
4. Outputs: `_vCBV_quad`, `_vCBV_z`, `_R2P_e1`, `_R2P_e2`, `_R2P_Fstat`. Run before and after z-shim correction on the test subjects and keep both.

## Phase 5 — CSF fraction in ASE space

1. `asePreProcess.sh`: add a pipeline key `qualityMaskNames` (default `anat_fs_xvent_mask_al.nii.gz`). Transform these like the segmentations, but onto the ASE grid with linear interpolation so the result is a fraction (0–1), not a label. Output `<name>_ase.nii.gz`.
2. `aseFit.sh`: pass the transformed mask path(s) to MATLAB as a new reserved option, e.g. `qualityMaskFiles`, alongside `segFiles`.

## Phase 6 — Quality module

1. New `CoreMatFiles/computeVoxelQuality.m`. Inputs: the fit struct (raw, unclamped maps), σ maps, Gz outputs, the vCBV and R2′ consistency stats, the CSF fraction, and `opts.quality`. Outputs: logical `keep` (X × Y × Z) and a uint16 `flags` bitmask.
2. Each criterion has its own enable flag and threshold in `quality`. A disabled criterion never sets its bit.
3. Bit definitions:

| Bit | Criterion | Default threshold |
| --- | --- | --- |
| 0 | Implausible raw values (vCBV, OEF, R2′, R2 outside bounds) | config bounds |
| 1 | k ≥ 0 or vCBV\_quad ≤ 0 | — |
| 2 | vCBV\_int vs vCBV\_quad disagreement | abs(z) > 3 (to calibrate) |
| 3 | Low SNR | SNR < 10 (to calibrate) |
| 4 | CSF fraction | > 0.2 |
| 5 | Large vein (high vCBV\_int) | vCBV > 0.15 |
| 6 | Z-shim: Gz at grid or moment-range edge | — |
| 7 | Z-shim: echo 1 vs echo 2 Gz disagreement | to calibrate |
| 8 | Z-shim: too many frames below Fmin | > 30% of long frames |
| 9 | R2′ echo inconsistency | F test p < 0.01 (to calibrate) |

4. `keep` = no enabled bit set. Thresholds marked "to calibrate" are placeholders until Phase 9.

## Phase 7 — ROI integration

1. In `ASE_fitting.m`, when `stats` is requested and quality is enabled, force the voxelwise fit to run first even if `maps`/`coeffs` are not requested.
2. Apply the mask to the membership matrix before `regionMeans`: `M = M .* keep(:)'`. Hard exclusion only for now.
3. Optional soft weighting behind `quality.noiseWeighting` (default false): multiply M by 1/σ² in log space (≈ S²/σ²). Keep it off until compared against plain occupancy weights.
4. Add `_ROI_Stats.csv` columns: `nVox_kept`, `fracExcluded`, the excluded fraction per criterion bit, and the median of the voxel maps (`OEF_med`, `R2P_med`, `vCBV_med`, `R2_med`) over kept voxels.
5. Below `quality.minVoxKept` (default 5 effective voxels), set the region's values to NaN.

## Phase 8 — Outputs, config, viewer

1. Write `_qualityMask` and `_qualityFlags`, plus every intermediate map from Phases 2–4, and list them all in the manifest.
2. Add the `"quality"` block to `ase_config.json` with every flag and threshold. Check that `runASEFromConfig.m` handles the nested struct (`convertSpecial` already recurses) and that `config_resolved.json` records it.
3. `ASE_Viewer.m`: add entries to `mapDefs` for the flags, the Gz maps, σ/SNR, vCBV\_quad and the R2′ F map. Show excluded voxels in ROI mode, e.g. a different overlay colour.

## Phase 9 — Validation and calibration

1. Run with quality enabled on the test subjects and compare against the Phase 0 baseline.
2. Calibrate each "to calibrate" threshold from its distribution in clean deep white matter. Then check the combined false-exclusion rate in that tissue (target under about 10%), not each test alone.
3. Report per ROI: the shift in OEF, R2′ and vCBV caused by each criterion on its own (selection-bias check).
4. Visual check in the viewer: flagged voxels should cluster at CSF borders, sinuses and large veins, not randomly.
5. Summarize the findings for the user before turning any criterion on by default.
