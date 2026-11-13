# Voxel-quality mask: overnight work report

Branch: `quality-mask` · Started 2026-11-12 20:55 · Test subject: `bch046_250816`

You asked me to work through the night without checking in, so I did **not** stop for approval at
each phase gate the pasted prompt asks for. Instead, each phase has its own commit, and the
summary I would have sent you is written below. Questions that I would have asked are collected in
**Open questions**. Where something was ambiguous I took the conservative option and flagged it,
rather than guessing silently.

---

## 0. Pre-flight

### Guardrails, restated

- **Git:** all work happens on the `quality-mask` branch. I never merge, rebase, push or force-push
  to main, never open-and-merge a PR, and never delete branches. Anything that would touch main
  waits for you. Each phase ends with one commit that names the phase.
- **Data:** I use only the local subject in `testData/`, with no cluster access. Original data
  (headers, processed volumes, earlier fit runs) is never overwritten or deleted. Every script
  runs against a scratch copy, with `--allSubjDir` pointed at that copy. I check what a script
  writes before running it (`--regenHdr` overwrites in place, and `aseFit.sh` writes fits and a
  failure log into `allSubjDir`).
- **Purpose:** the goal is to check that the code changes work, not to produce datasets. Every
  output is a disposable diagnostic. I report numbers in this document and delete each phase's
  scratch outputs once they're summarized. Only the Phase 0 baseline stays until the end.
- **Defaults unchanged:** every new feature sits behind `fit.quality.*` flags that are off by
  default.

### Branch

`quality-mask` already existed and pointed at the same commit as `main` (`d235598`), with a clean
tree, so I used it as is. `git merge-base main quality-mask` = `d235598`.

### Scratch copy

- I copied `testData/bch046_250816/session_ase/` without its `fits/` folder, plus `dcm.index.dat`,
  the subject's `datalist_ase` row and `GM_WM_long.txt`, to `testData/scratch/qm/` (287 MB).
  I didn't copy `raw_data/` or `session_afni_nortc/`: no stage I run needs them.
- Every run uses `--allSubjDir testData/scratch/qm` and
  `--segLUTs testData/scratch/qm/GM_WM_long.txt`. The config's `segLUTs` default is a `/autofs`
  cluster path, so it has to be overridden.
- **Already present:** `testData/scratch/bch046_250816/` was there before I started (created about
  20:22–20:38, partial `session_ase`, no raw NIfTIs). I didn't create it, so I left it untouched.
  Decide whether to delete it.
- Before any run, I recorded a size+mtime fingerprint of every original file (1104 files). I
  re-check it at the end (see **Final checks**).

### Plan files

The prompt says to read `testData/plan/` and several reasoning files. Only
`testData/plans/Action plan.md` exists (the folder is `plans`, not `plan`). The reasoning docs
(vCBV consistency, repetition noise, z-shim, R2′ echo check) **are not on disk**, so I worked from
the action plan plus the code.

---

## Phase 0 — Baseline

**Command:**
```
/usr/bin/time -l bash BashFiles/aseFit.sh --allSubjDir testData/scratch/qm \
    --singleSubjID bch046_250816 --segLUTs $PWD/testData/scratch/qm/GM_WM_long.txt
```

**Runtime and memory:**
- 19.5 s wall time (MATLAB startup included).
- Peak memory 5.10 GB: `maximum resident set size` of the process tree.

**Data:**
- 128×128×30 voxels, 98 frames per echo.
- 2 reps × 45 sweep frames plus 8 z-shim frames.
- TE = 61.2 and 93.2 ms. `SliceThickness` = 4.8 mm.

**Regimes:**
- t_crit (×1.5) = 17.36 ms. Short regime 0–14 ms (30 frames). Long regime 21–76 ms.

**Baseline maps** (over all finite, nonzero voxels in the full FOV; no brain mask):

| map | n | median | p5 | p95 |
|---|---|---|---|---|
| OEF | 187,731 | 0.200 | 0.034 | 0.929 |
| R2′ (s⁻¹) | 350,069 | 3.39 | 0.13 | 32.2 |
| vCBV | 212,303 | 0.080 | 0.004 | 0.872 |
| R2 (s⁻¹) | 372,111 | 3.76 | 0.14 | 21.6 |
| imErr | 187,731 | 0.063 | 0.012 | 0.275 |

**Baseline ROI stats:**
- 160 regions. Median ROI values: OEF 0.268, R2′ 8.98, vCBV 0.125, R2 12.0; median nVox 275.
- 15 ROIs have NaN OEF. The max ROI OEF is 1 (clipped).
- vCBV 0.125 is high for tissue, which suggests contamination from large veins and CSF.

| ROI | OEF | R2′ | vCBV | R2 | nVox |
|---|---|---|---|---|---|
| Left-UnsegmentedWhiteMatter | 0.308 | 5.92 | 0.089 | 9.48 | 1679.8 |
| Left-Thalamus-Proper | 0.227 | 9.01 | 0.184 | 11.60 | 430.9 |
| ctx-lh-precuneus | 0.270 | 6.04 | 0.104 | 12.43 | 611.5 |
| ctx-lh-superiorfrontal | 0.250 | 7.16 | 0.133 | 11.61 | 1166.5 |

**Surprises:**
1. **Echo-1 zero-moment z-shim frame used as a long-regime sample.** It is frame 96
   (dTE = 44 ms, ZMoment = 0), and the existing code uses it because it masks only
   `zShim ~= 0`. The same applies to echo 2. This is probably fine (it is a real 44 ms frame),
   but it's worth knowing.
2. **Pre-existing fit differs.** The fit already in `bch046_250816/session_ase/fits/20261112_151816`
   has different ROI values from today's baseline. It was presumably made by the code before
   commit `d235598`. I didn't investigate further.

Kept: `testData/scratch/qm/bch046_250816/session_ase/fits/20261112_205823` (baseline).
Deleted: nothing.

---

## Phase 1 — Header metadata

**Changes:**
- `BashFiles/asePreProcess.sh`, `build_header`: reads `SliceThickness` from each echo's dcm2niix
  sidecar using `json_num`, and writes `"SliceThickness_mm"` into each echo block, right after
  `EchoTime_ms`. If the value is missing, the field is left out and a warning is added; this is
  not an error.
- `Aux_Mat_Functions/readASEHeader.m`: new fields `info.zShimInds` (trailing frame indices),
  `info.zMoments` (nZshim × nEchos, Tms/m × 1e-3 → T·s/m) and `info.SliceThickness_mm`
  (nEchos × 1, NaN if absent).
- `readASEHeader.m` also gets a small robustness fix. If the echo blocks have different optional
  fields, `jsondecode` returns a cell and `[e{:}]` failed. Missing fields are now padded.
- The "error only when z-shim correction is enabled" check is enforced in `ASE_fitting.m` in
  Phase 3.

**Ran:** `asePreProcess.sh --allSubjDir testData/scratch/qm --singleSubjID bch046_250816 --regenHdr`
(scratch only). Before running, I checked that with the existing outputs present only the header
stage runs. The log confirmed it ran only "Building header".

**Checks:**
- **Header diff:** a `jq` diff of the header before and after, ignoring `created`, shows exactly
  two added lines, `"SliceThickness_mm": 4.8` for each echo.
- **MATLAB values:** `SliceThickness_mm = [4.8 4.8]`, `zShimInds = 91:98`.
- **zMoments:** echo 1 = [3.819 3.055 2.291 1.527 0.763 0 −0.763 −1.527]·1e-6 T·s/m;
  echo 2 = [6.596 … −2.638]·1e-6.
- **Old header** (no field): `[NaN NaN]`, no error.
- **Mixed header** (field on echo 1 only): `[4.8 NaN]`, no error. Before the fix, this case
  crashed.

**Surprises:**
- The pre-existing warning "ZMoment values for echo 2 differ from echo 1" is expected: echo 2's
  shim moments are larger by about 76/44, consistent with a longer dephasing time. Phase 3 relies
  on this.
- The zero-moment shim frame is the 6th of 8, so the sampled moment range is asymmetric
  (−1.5 to +3.8 µT·s/m for echo 1).

**Deleted from scratch:** `prepASE_failures_*.txt`. The regenerated scratch header stays, because
it is the input for later phases.

---

## Phase 2 — Repetition noise

**Changes:**
- **New `Aux_Mat_Functions/estimateRepNoise.m`:** estimates σ as the residual of the additive
  rep × point model, with (nSweep−1)(nReps−1) degrees of freedom.
  - For 2 reps this equals Σ(d − d̄)² / (2(nSweep−1)), where d = rep1 − rep2. That is the plan's
    formula, except I divide by n−1 rather than n, because the mean drift d̄ has been subtracted.
    The difference is 2% at n = 45.
  - For more than 2 reps it pools the per-point variances after removing each rep's offset, as
    the plan asks.
  - It also returns S0, the rep-mean signal at the dTE = 0 sweep positions. The SNR map is S0/σ.
    "Baseline signal" was not defined further in the plan, so I chose the spin-echo frames.
- **New `Aux_Mat_Functions/qualityDefaults.m`:** the single source of the `quality` block
  defaults, covering every flag and threshold for Phases 2–9. Everything is off by default,
  unknown keys raise an error, and partial blocks are filled in.
- **`CoreMatFiles/ASE_fitting.m`:**
  - New `opts.quality`.
  - σ is computed from the unaveraged data. In stacked mode it is sliced from the stack; with
    `averageReps` each echo file is read again.
  - Outputs `_sigma_e1/_e2` and `_SNR_e1/_e2`, written by the new `writeQualityMaps` with
    manifest entries.
  - σ is needed when any of `quality.repNoise`, `fitExtensions`, `noiseWeighting`, or
    (`enabled` and `useSNR`) is set.
- **New `Tests/test_estimateRepNoise.m`:** 4 tests, all passing.
  - Known σ with rep drift: recovered within 3%.
  - 3 reps: recovered within 3%.
  - The 2-rep pooled form matches the plan formula exactly.
  - A single rep returns NaN with a warning.

**Runs** (scratch, about 12 s and 4.4–5.2 GB each):

| run | result vs baseline |
|---|---|
| defaults | all 9 maps and the ROI CSV bit-identical |
| `--fit 'quality={"repNoise":true}'` | existing outputs bit-identical; adds 4 maps |
| same + `averageReps=true` | σ bit-identical to the stacked run |

**Numbers:**
- **Intensity mask:** S0 > 0.2 × p99, a crude stand-in for a brain mask. I used it because
  ASE-grid brain masks don't exist yet.

| echo | median σ | median SNR | p5–p95 SNR | fraction SNR < 10 |
|---|---|---|---|---|
| e1 | 9.5 | 102 | 21–218 | 1.5% |
| e2 | 8.3 | 77 | 20–156 | 1.9% |

**Surprise:**
- In the background (S0 < 2% of p99), rep-difference σ is about 2. The Rayleigh estimate from the
  background mean is about 26. So the background is not thermal-noise-like. It is probably ghost
  or interpolation structure from Fourier motion correction, or scanner-side filtering.
- In tissue, σ is about 4–5× the background value. That means the rep-difference σ captures
  physiological and motion noise, not only thermal noise. That is what we want for the fit
  uncertainty, but **SNR < 10 will rarely trigger in tissue**, so the threshold needs calibrating
  in Phase 9.

**Deleted from scratch:** the 3 Phase 2 fit folders and their logs. The baseline is kept.

---

## Phase 3 — Z-shim gradient estimate and correction

**Changes:**
- **New `CoreMatFiles/fitZshimGradient.m`:**
  - Vectorized grid search: one matrix product Y·F gives R_j for all voxels and all 301 grid
    points.
  - Parabolic refinement, then the amplitude is recomputed at the refined Gz.
  - Outputs Gz, A, R, the attenuation ratio y(M=0)/A, and an `edge` bitmask
    (1 = grid edge, 2 = peak outside the sampled moments).
  - Uses the unnormalized `sin(x)./x` with `x==0 → 1`, never MATLAB's `sinc`.
- **New `CoreMatFiles/zshimAttenuation.m`:** F = |sinc(γΔz·Gz·t/2)|.
- **`ASE_fitting.m`:**
  - The `dBCorr` error is gone. Frames with nonzero z-moment are now always excluded from the
    regime fits.
  - `dBCorr=true` runs the correction: each frame is divided by F, where F uses that frame's own
    echo's Gz and t = |`taus_full`|.
  - Where F < `quality.zshim_Fmin`, the voxel–frame value becomes NaN. This is a new voxel × frame
    exclusion.
  - Voxels with no Gz estimate (no signal) are left uncorrected (F = 1).
  - `dBCorr` without `SliceThickness_mm` or without z-shim frames raises an error. The
    estimate-only path just warns and skips.
  - `quality.zshimEstimate` alone writes the maps without changing the fit.
  - New per-echo outputs `_Gz`, `_zshimAmp`, `_zshimAtten` and `_zshimEdge`, plus
    `_zshimFracLow`: the fraction of the used long frames with F < Fmin, which feeds bit 8. All
    are in the manifest.
- **New `Tests/test_fitZshimGradient.m`** (5 tests, all passing; Phase 2 tests still pass):
  - Gz ∈ {−30, −10, 0, 12.3, 40, 75} µT/m at peak SNR 80, using the real echo-1 moments, 44 ms
    and 4.8 mm: median and p90 error under 2 µT/m, amplitude within 2%.
  - Noiseless recovery is exact (to 0.1).
  - Both edge bits fire.
  - Background returns NaN.
  - The attenuation is exactly 0 at the first null of sin(x)/x.

**Runs:**
- Estimate only: all baseline outputs bit-identical, 9 maps added.
- `dBCorr=true`: 12.6 s, 6.1 GB.

**Brain-mask diagnostics.** These use a resampled FreeSurfer brain mask (fraction > 0.9,
81,335 voxels) and aseg WM/GM, made with `3dAllineate` from the existing inverse matrix and kept
in `testData/scratch/qm/diag/`.

| | e1 | e2 |
|---|---|---|
| median Gz in brain (µT/m) | −4.6 | −6.6 |
| median attenuation y(0)/A | 0.985 | 0.928 |
| grid-edge / outside-moment flags | 0.1% / 8.3% | 0.7% / 14.3% |

**Model check on the data:**
- The mean normalized z-shim profile over brain voxels at the largest moment is 0.262 (e1) and
  0.203 (e2).
- The |sinc| model with Δz = 4.8 mm predicts 0.26 and 0.21. The rect-profile model and the
  units are consistent.

**Echo convention check.** I fitted echo 2 under alternative conventions:

| echo 2 convention | median \|Gz₂−Gz₁\| | corr | edge |
|---|---|---|---|
| **M₂, t = 76 ms (plan)** | **3.3** | **0.50** | **14%** |
| M₁+M₂ (cumulative), 76 ms | 7.2 | 0.51 | 13.5% |
| M₂, t = 44 ms | 10.2 | 0.63 | 10.6% |
| −M₂ | 26.4 | −0.50 | 14% |

The plan's convention fits best. The **absolute sign of Gz is not identifiable** from a single
echo (M → −M only flips Gz, with identical residuals), but the relative sign between echoes is.
That doesn't matter for the correction, because |sinc| is even.

**Effect of correction** (brain medians, baseline → corrected):

| map | baseline | corrected |
|---|---|---|
| R2′ (s⁻¹) | 7.54 | **5.26** (WM 6.01 → 4.52) |
| vCBV | 0.143 | **0.104** |
| OEF | 0.220 | 0.224 |
| R2 (s⁻¹) | 12.44 | 12.66 |

R2′ and vCBV fall substantially and OEF is nearly unchanged, which is what you'd expect if
through-slice dephasing had been inflating R2′ and vCBV together.

**Surprises and concerns:**
1. **6.1% of brain voxels lose their fit after correction** (R2′ finite: 100% → 93.9%).
   - 14.6% of brain voxels have at least one long frame below Fmin, and 11.4% have more than 30%.
   - Echo-2 frames all have t ≥ 32 ms, so in high-Gz voxels the whole of echo 2 is NaN'd. Then
     the long design [1, −TE, −τ] is rank-deficient (one TE), and `fitMultiRegressor` returns NaN.
   - Bit 8 is meant to flag exactly these voxels, but they currently show up as missing values
     instead.
   - **Open question Q1.**
2. **ROI path with correction.** `regionMeans` averages each frame over the voxels that are still
   finite. After the voxel × frame NaNs, long-τ frames are averaged over a different
   (low-Gz-only) subset of voxels than short-τ frames. That is a selection bias inside the ROI
   curve. Phase 7's voxel mask helps, because bit-8 voxels are removed entirely. **Open
   question Q2.**
3. **Whole-FOV medians are misleading.** Without a brain mask, the run log reported echo-2 median
   Gz = −20 and 75% edge-flagged, because background voxels dominate. The verbose line is now
   just informational.

**Deleted from scratch:** the 2 Phase 3 fit folders and logs. Kept:
`testData/scratch/qm/diag/` (resampled brain mask and aseg), which later diagnostics need. It
will be deleted at the end.

---

## Phase 4 — Fit extensions

**Changes:**
- **`fitASERegimes`** (in `ASE_fitting.m`):
  - Also returns `k` = `BS(:,:,:,2)` and `vCBV_quad = −0.3·R2′²/k` (NaN where k ≥ 0) whenever
    the short design has a τ² column.
  - It also returns its design matrices, so the error propagation uses exactly the fitted design.
  - The ROI path gets the same fields; they aren't used there yet.
- **New `CoreMatFiles/vcbvConsistency.m`:** analytic standard errors of vCBV_int and vCBV_quad,
  and the z-score of their difference.
  - It uses per-voxel OLS projection rows P = (XᵀWX)⁻¹XᵀW, built with `pageinv`/`pagemtimes` in
    voxel chunks, so each voxel's own valid-frame set (z-shim NaNs included) is respected.
  - Log-space noise per frame is σ_echo/S, scaled by 1/√nReps for averaged frames.
  - It includes the **coupling through R2**: the short fit uses ln S + R2·TE with R2 from the
    long fit, which correlates vCBV_int and k with the long data.
  - vCBV_quad uses the delta method.
- **New `fitR2PerEcho`** (local to `ASE_fitting.m`):
  - Separate long fits per echo, ln S = c_e − R2′_e·τ.
  - RSS = rms²·n from `fitMultiRegressor`.
  - F = (RSS_joint − RSS_split)/(RSS_split/(n−4)) with df (1, n−4). The split model has 4
    parameters and the joint model 3, nested. The p-value is kept internally for bit 9.
- **Outputs:** `_vCBV_quad`, `_vCBV_z`, `_vCBV_SE`, `_vCBVquad_SE`, `_R2P_e1`, `_R2P_e2`,
  `_R2P_Fstat`, all in the manifest. They run when `quality.fitExtensions` is set, or when
  `quality.enabled` is set with any criterion that needs them.
- **New `Tests/test_vcbvConsistency.m`:** Monte Carlo over 20,000 voxels with the real protocol
  (2 reps, 2 echoes, σ = 4/3.5).

| quantity | Monte Carlo | analytic |
|---|---|---|
| SE(vCBV_int) | SD 0.01073 | **0.01076** |
| SE(vCBV_quad) | robust SD 0.0102 | **0.0099** |
| z | median −0.16, robust SD 0.92, \|z\|>3 in 0.18% | — |

The small bias and the narrower spread of z come from the delta-method approximation of a
ratio. The tolerances are set to that level, and all 10 repo tests pass.

**Runs** (each about 15 s, 6.6–7.0 GB peak, up from 5.1 GB):
- Extensions without correction: all baseline outputs bit-identical.
- Extensions with `dBCorr=true`.

**Real-data diagnostics** (brain mask; WM in brackets):

| | no correction | with dBCorr |
|---|---|---|
| vCBV_int median | 0.143 [0.118] | 0.104 [0.090] |
| vCBV_quad median | 0.182 [0.151] | 0.118 [0.101] |
| fraction with k < 0 (vCBV_quad defined) | 73.8% [82.0%] | 64.2% [74.8%] |
| SE_int / SE_quad median | 0.0044 / 0.037 | 0.0039 / 0.027 |
| z median, robust SD | −0.75, 1.29 | −0.54, 1.31 |
| \|z\| > 3 | 22.1% [18.5%] | 22.1% [20.3%] |
| R2′ joint / e1 / e2 median | 7.54 / 5.23 / 7.77 | 5.26 / 4.00 / 5.35 |
| F-test p < 0.01 | **78.8%** [82.9%] | **63.8%** [72.2%] |

**Surprises and concerns:**
1. **The R2′ echo F-test rejects in most of the brain.**
   - R2′ from echo 2 (τ = 52–76 ms) is consistently about 2 s⁻¹ higher than from echo 1
     (τ = 21–44 ms).
   - That fits a slope that is still increasing toward the asymptotic R2′. Echo-1's long frames
     start at 1.2 × t_crit, where the curve is not yet linear.
   - With about 130 frames, even small systematic misfit becomes very significant.
   - **At its placeholder threshold, bit 9 would exclude about 80% of WM.** It needs a
     calibrated threshold, or an effect-size criterion such as |R2′_e1 − R2′_e2| relative to its
     SE. **Open question Q3.**
2. **Bit 1 (k ≥ 0) would exclude about 26% of the brain** (18% of WM), and |z| > 3 another
   20%. vCBV_quad is noisy (SE about 0.03, roughly 8× SE_int), and its median exceeds vCBV_int.
   z is overdispersed (robust SD 1.3, not 1), which suggests a systematic model difference on
   top of noise. Calibration is needed (Phase 9).
3. Correction moves vCBV_quad toward vCBV_int (0.18 → 0.12 vs 0.10), and roughly halves the
   F-statistic. That is consistent with through-slice dephasing being part of the misfit.
4. **Memory:** peak rises by about 2 GB with all diagnostics on. Runtime is unchanged.

**Deleted from scratch:** the 2 Phase 4 fit folders and logs. The plan says to "keep both" the
pre- and post-correction runs. I kept them as the numbers above, in line with the
nothing-on-disk rule.

---

## Phase 5 — CSF fraction in ASE space

**Changes:**
- **`asePreProcess.sh`:**
  - New pipeline key `qualityMaskNames`, default `anat_fs_xvent_mask_al.nii.gz`. It's also a
    `--set` key and in the config.
  - The masks are added to the transfer list (deduplicated).
  - New stage 8, `stage_qualityMasks`:
    `3dAllineate -1Dmatrix_apply <tformApplyName> -master <baseName> -final linear` produces
    `<name><segSuffix>.nii.gz` on the ASE grid.
  - Planning: the stage runs when the output is missing, or on `--regenTransform`, re-registration
    or re-transfer, like the segmentation transforms.
- **`aseFit.sh`:**
  - New key and flag `qualityMaskNames`.
  - Masks that exist are passed as `'qualityMaskFiles'`. Missing ones only add a warning, so
    default runs never fail on them.
- **`runASEFromConfig.m`:** `qualityMaskFiles` is a reserved per-subject option. It's passed
  through and recorded in `config_resolved.json` (verified).
- **`ASE_fitting.m`:**
  - New `opts.qualityMaskFiles`. The geometry must match the ASE data (`isequalNiftiGeometry`,
    error on mismatch).
  - The masks are read **only when a quality feature is on**, so the default output can't change.
  - Writes `_csfFrac`.

**Mask semantics (decision, open question Q4).** `anat_fs_xvent_mask_al` is **not a CSF mask**.
It is 1 throughout the brain except the ventricles: the mean value is 0 in aseg 4/5/14/15/43/44
and 1 in WM, cortex, thalamus and aseg 24. So I treat each quality mask as a *tissue* fraction,
and CSF fraction = 1 − fraction. With several masks, the maximum CSF fraction is used. Outside the
brain the CSF fraction is 1.

**Runs:**
- `asePreProcess.sh` (no regen flags) on scratch: only the new stage ran, writing
  `anat_fs_xvent_mask_al_ase.nii.gz`.
- Default fit: bit-identical to baseline, even though `qualityMaskFiles` is now passed.
- `repNoise` fit: writes `_csfFrac`.

**CSF fraction > 0.2 by region** (aseg resampled to the ASE grid with NN):

| region | n | fraction flagged |
|---|---|---|
| ventricles (4/43) | 1,496 | 100% |
| WM (2/41) | 29,398 | 0.3% |
| cortex | 28,714 | 14.1% (edge partial volume against outside-brain) |
| thalamus | 914 | 1.6% |
| aseg "CSF" (24) | 73 | 21.9% |
| whole brain (mask > 0.9) | 81,335 | 18.6% |

**Notes:**
- The source mask is float, so the linear resampling gives true fractions. A byte mask (like
  `anat_fs_brain_al_mask`) would be rounded back to 0/1 by `3dAllineate` unless it's converted
  first. That matters if you add byte masks to `qualityMaskNames`.
- Sulcal CSF counts as tissue in this mask, so bit 4 catches ventricles and the brain's outer
  boundary, not sulci.

**Deleted from scratch:** the 2 fit folders, logs and the prep log. Kept: the new
`anat_fs_xvent_mask_al_ase.nii.gz` in the scratch session folder, which is an input for later
phases.

---

## Phase 6 — Quality module

**Changes:**
- **New `CoreMatFiles/computeVoxelQuality.m`:** `[keep, flags, bitInfo]`, with a uint16 bitmask
  using the plan's bits 0–9.
  - Each bit has its own `use*` flag and threshold in `quality` (`qualityDefaults.m`).
  - A disabled criterion never sets its bit.
  - An enabled criterion whose input is missing (for example, no quality mask) warns
    (`computeVoxelQuality:missingInput`), is reported as `Available=false` and sets nothing.
  - `keep` = no bit set.
- **Interpretation choices (flagged):**
  - **Bit 0:** raw (unclamped) R2/R2′/vCBV/OEF outside the `physBd_*` bounds, **or NaN fit
    values** (except OEF, which is NaN by construction where vCBV ≤ 0 or R2′ ≤ 0, cases already
    caught).
  - **Bit 1:** `~(vCBV_quad > 0)`, which covers k ≥ 0 because vCBV_quad is NaN there.
  - **Bit 3:** uses the **lowest SNR over echoes**.
  - **Bit 6:** either edge bit, in either echo.
  - Bits 2, 7 and 9 use the placeholder thresholds |z| > 3, |ΔGz| > 20 µT/m and p < 0.01.
    I picked 20 µT/m because the plan gives no number.
- **`ASE_fitting.m`:** with `quality.enabled`, it computes `qual.keep`/`qual.flags` from the raw
  maps and prints the bit table when verbose. **This phase has no effect on outputs** (ROI use is
  Phase 7, writing is Phase 8).
- **New `Tests/test_computeVoxelQuality.m`:** 4 tests, and all 18 repo tests pass.
  - 11 synthetic voxels, each breaking exactly one criterion: each sets exactly its own bit.
  - Disabling each criterion in turn never sets its bit.
  - A missing input warns and skips.
  - A NaN fit sets bit 0.

**Run:** `--fit 'quality={"enabled":true}'` runs end to end, with all inputs available
(CSF mask, σ, Gz, z, F).
- **Whole-FOV result:** 6,608 of 483,421 fitted voxels kept.
- **Not meaningful:** these counts are dominated by background voxels, which have no ROI
  membership anyway. The brain-restricted per-bit breakdown is in Phase 9, once Phase 8 writes
  the flags.

**Deleted from scratch:** the Phase 6 fit folder and log.

---

## Phase 7 — ROI integration

**Changes** (all in `ASE_fitting.m`):
- **Forced voxel fit:** with `stats` requested and `quality.enabled`, the voxelwise fit always
  runs. Maps and coeffs are still written only if requested.
- **Hard exclusion:** before `regionMeans`, `M = M .* keep(:)'`.
- **Noise weighting** (`quality.noiseWeighting`, default off): M is also multiplied by the
  echo-1 SNR² (≈ S²/σ², i.e. 1/var(log S)), normalized by its median over kept voxels.
  - The plan didn't say which echo or frame provides S and σ, so I used the baseline (dTE = 0)
    signal of echo 1.
- **`minVoxKept`:** regions with fewer kept effective voxels than `quality.minVoxKept` (5) get
  NaN.
- **New ROI columns** (only when quality is enabled; the default CSV is unchanged):
  - `nVox_kept` and `fracExcluded`, both occupancy-weighted.
  - `fracExcl_b<k>_<name>` for every enabled bit. Bits overlap, so they don't sum to the total.
  - `OEF_med`, `R2P_med`, `vCBV_med`, `R2_med`: medians of the raw voxel maps over each region's
    kept voxels.

**Runs and checks:**
- **Default:** maps and ROI CSV bit-identical to baseline.
- **`quality.enabled`, all `writeOpts`:** 15.7 s, 6.9 GB.
- **`quality.enabled` with `writeOpts=["stats"]`:** only the CSV, manifest and config are
  written, and the ROI table is **identical** to the full run. This confirms the forced voxel
  fit.
- **Voxel-median vs ROI-fit consistency** over kept voxels (medians across regions):

| | voxel median | ROI fit |
|---|---|---|
| OEF | 0.322 | 0.334 |
| R2′ | 4.58 | 4.78 |
| vCBV | 0.069 | 0.066 |

  These agree well. Mean-signal fitting and voxel medians give similar answers once the bad
  voxels are gone.
- **Noise weighting vs occupancy weights** (111 regions): median |ΔOEF| 0.007, |ΔR2′| 0.25,
  |ΔvCBV| 0.003. That's a small effect, so it stays off as planned.

**Results with every criterion at its placeholder threshold** (no dBCorr), for illustration only:
- **91.4% of ROI occupancy excluded.** 49 regions fall below `minVoxKept`, against 15 NaN-OEF
  regions at baseline.
- **Per-bit median excluded fraction per region:**

| bit | criterion | median fraction |
|---|---|---|
| 9 | r2pEcho | **0.80** |
| 5 | vein | 0.42 |
| 4 | csf | 0.20 |
| 1 | quadSign | 0.17 |
| 2 | vcbvZ | 0.17 |
| 6 | zshimEdge | 0.12 |
| 0 | implausible | 0.11 |
| 8 | zshimFrac | 0.08 |
| 7 | zshimEcho | 0.07 |
| 3 | lowSNR | 0.00 |

- **Shift across the 107 regions with values in both runs:** OEF 0.263 → 0.334, R2′ 7.42 → 4.78,
  vCBV 0.125 → 0.066.
- **These are not results.** They show why the defaults must stay off until Phase 9. Bit 9 and
  the vein threshold dominate, and the vein bit uses vCBV that hasn't been through z-shim
  correction.

**Deleted from scratch:** the 4 Phase 7 fit folders and logs.

---

## Phase 8 — Outputs, config, viewer

**Changes:**
- **`ASE_fitting.m`:**
  - Writes `_qualityMask` (1 = kept) and `_qualityFlags` (bitmask; exact integers in single
    precision).
  - Every Phase 2–4 intermediate map is written as well, each with a manifest entry. The
    manifest has 33 path entries with everything on, and no map on disk lacks an entry.
- **`BashFiles/ase_config.json`:**
  - New `fit.quality` block with all 26 flags and thresholds, every feature off. It matches
    `qualityDefaults()` exactly (checked in MATLAB).
  - New pipeline key `qualityMaskNames` (from Phase 5).
- **`runASEFromConfig.m`:**
  - `convertSpecial`/`encodeSpecial` already recurse, so the nested block works as is.
  - **Fix found while testing:** `--fit 'quality={...}'` used to *replace* the whole quality
    block, so `config_resolved.json` recorded only the override keys. Nested struct overrides
    are now **merged key by key**, and the block goes through `qualityDefaults`. That validates
    the keys and records the full 26-key block.
  - Verified with `--fit 'quality={"enabled":true,"snr_min":7}'`: the resolved config and the
    manifest both show `snr_min=7`, `csf_max=0.2`, 26 keys.
  - A typo (`"enabeld"`) fails cleanly: "Unknown quality option(s): enabeld", recorded in the
    failure log. The `aseFit.sh --help` text documents the merge.
- **`ASE_Viewer.m`:**
  - `loadManifest` now loads every map listed in `mapDefs` that the manifest has, instead of a
    hard-coded list.
  - New `mapDefs` entries: qualityMask, qualityFlags, vCBV_quad, vCBV_z, R2P_e1/e2, R2P_Fstat,
    Gz_e1/e2, zshimFracLow, sigma_e1/e2, SNR_e1/e2, csfFrac.
  - **ROI mode:** when a `qualityMask` is loaded, excluded member voxels are drawn **in cyan**
    over the red ROI overlay, and the ROI curve is computed from kept voxels only, as in the fit.
    Optional noise weighting is not reproduced in the viewer.

**Checks:**
- **Default run:** bit-identical to baseline.
- **`quality.enabled` + `dBCorr`:** flags are integers, max 1021, and `mask == (flags==0)`
  exactly.
- **Viewer** (driven from `matlab -batch`):
  - Loaded all 22 maps.
  - Opened the window and cycled through all 22 volume entries without error.
  - Selected `Left-Thalamus-Proper`, and the cyan excluded-voxel overlay was created.
  - I couldn't *look* at the rendering, so a visual check is still needed (Phase 9.4).
- **Side note:** manifest `paths` are stored without the `.nii` extension, so `isfile` on them
  fails. This was already the case before the branch (for example `CoeffsLong`); the viewer
  copes. Not changed.

**Deleted from scratch:** 4 fit folders, logs and the failure log from the typo test.
