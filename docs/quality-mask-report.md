# Voxel-quality mask: overnight work report

Branch: `quality-mask` · Started 2026-11-12 20:55 · Test subject: `bch046_250816`

You asked me to work through the night without checking in, so I did **not** stop for approval at
each phase gate the pasted prompt asks for. Instead, each phase has its own commit, and the
summary I would have sent you is written below. Questions that I would have asked are collected in
**Open questions**. Where something was ambiguous I took the conservative option and flagged it,
rather than guessing silently.

## TL;DR

- **All 10 phases (0–9) are implemented and committed** on `quality-mask`.
  - Every feature sits behind `fit.quality.*` (or `dBCorr`) and is off by default.
  - Default outputs are bit-identical to the baseline at every phase.
  - 14 MATLAB unit tests (`Tests/`) all pass.
- **Unit tests confirm three estimators:**
  - The z-shim gradient estimator recovers Gz within 2 µT/m at realistic SNR.
  - The repetition-noise estimator recovers σ within 3%.
  - The analytic vCBV SEs match Monte Carlo.
- **On real data, z-shim (dB) correction matters most.** Brain R2′ drops from 7.5 to 5.3 s⁻¹ and
  vCBV from 0.14 to 0.10, while OEF barely changes. The cost is about 6% of brain voxels losing
  their fit (Q1).
- **Calibration** (one subject, clean deep WM):
  - Two placeholder criteria don't work as specified. The **R2′ echo F-test** (bit 9) flags
    64–80% of the brain and biases OEF by +0.09. **k ≥ 0** (bit 1) is noise-dominated.
  - A calibrated set (|z| > 8.9, |ΔGz| > 12.7 µT/m, SNR < 10, vein > 0.59, with bits 1 and 9
    off) excludes **3.3% of clean deep WM** and about 43% of the brain.
  - The excluded voxels cluster at CSF borders and susceptibility regions, not randomly.
- **Nothing should be turned on by default yet.** 8 open questions are at the end.

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
- ~~Sulcal CSF counts as tissue in this mask, so bit 4 catches ventricles and the brain's outer
  boundary, not sulci.~~ **Corrected 2026-11-16:** this was wrong. The mask is the FreeSurfer
  aseg labels minus the ventricles, so sulcal CSF (aseg = 0) is 0 in it and counts as CSF. See Q4
  under Decisions.

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

---

## Phase 9 — Validation and calibration

**Code change in this phase:** the verbose z-shim summary lines in `ASE_fitting.m` now cover
only voxels with real signal (echo-1 shim amplitude > 0.2 × p99, which covers 94% of the brain
mask). Before, background voxels dominated them. They now read e1/e2 median Gz −5.1/−6.6 µT/m,
which agrees with the brain-mask analysis. Defaults are still bit-identical, and all 14 unit
tests pass.

### Setup

**Runs:**
- All criteria on at placeholder thresholds, with and without `dBCorr`.
- One run per criterion on its own, plus a reference with every criterion off. All use `dBCorr`
  and stats only.
- Calibrated "set C" (see below).
- Everything on with `averageReps=true`, to exercise the other frame layout: it runs, Gz is
  identical, and the excluded occupancy is about the same (88.0% vs 88.6%).

**Clean deep WM** (the calibration tissue):
- aseg 2/41 resampled to the ASE grid, eroded with a 3×3×3 cube.
- CSF fraction = 0 and brain fraction > 0.9.
- 3,233 voxels.
- Because **only one subject is available**, the thresholds are calibrated and evaluated on the
  same data. Expect them to be optimistic.

### 9.2 Thresholds

**Per-bit exclusion at the placeholder thresholds**, as % of voxels:

| bit | criterion | clean deep WM | WM | GM | deep GM | brain |
|---|---|---|---|---|---|---|
| 0 | implausible | 1.0 | 12.8 | 22.5 | 21.9 | 19.7 |
| 1 | quadSign | 5.0 | 25.2 | 34.2 | 18.7 | 35.8 |
| 2 | vcbvZ (\|z\|>3) | 9.3 | 15.2 | 16.6 | 14.4 | 14.2 |
| 3 | lowSNR (<10) | 0.0 | 0.1 | 0.8 | 0.5 | 0.9 |
| 4 | csf (>0.2) | 0.0 | 0.3 | 12.9 | 2.1 | 18.6 |
| 5 | vein (>0.15) | 6.6 | 19.6 | 30.7 | 39.3 | 29.3 |
| 6 | zshimEdge | 0.8 | 9.4 | 20.3 | 12.7 | 15.8 |
| 7 | zshimEcho (>20) | 0.6 | 5.1 | 13.9 | 13.2 | 10.9 |
| 8 | zshimFrac (>0.3) | 0.6 | 6.4 | 14.8 | 11.8 | 11.4 |
| 9 | r2pEcho (p<0.01) | **81.1** | 72.1 | 60.3 | 54.0 | 63.8 |
| — | **kept** | 15.3 | 13.7 | 12.2 | — | 10.4 |

This table uses `dBCorr`. Without correction the vein rate in clean deep WM is 13.9% and
r2pEcho is 91%.

**Distributions in clean deep WM** (with `dBCorr`):
- **|vCBV_z|:** p50 0.88, p90 2.98, p95 4.33, p99 **8.87**. It is NaN in 5% (k ≥ 0).
- **min-echo SNR:** p1 58, p5 77, p50 133.
- **|Gz_e1 − Gz_e2|:** p50 1.6, p95 4.3, p99 **12.7** µT/m.
- **R2′ F statistic** (df 1,150): p50 **25**, p99 156. The p-value at p99 is about 0.
- **vCBV_int:** p50 0.079, p95 0.164, p99 0.594, a heavy tail.
- **|R2′_e1 − R2′_e2| / R2′:** p50 0.29, p90 0.47, p99 0.88.

**Calibration verdicts.** The rule was a threshold at p99 of clean deep WM, so each criterion
alone falsely excludes about 1%.
- **Bit 2 (vcbvZ):** threshold **|z| > 8.9** (7.5 without correction), against the
  placeholder 3. z is overdispersed in real tissue: robust SD about 1.3, plus a heavy tail.
- **Bit 7 (zshimEcho):** threshold **|ΔGz| > 12.7 µT/m**, against 20.
- **Bit 3 (lowSNR):** **can't be calibrated this way.**
  - Clean WM has no low-SNR voxels. Its p1 (SNR 58) would flag 56% of deep GM, which has
    naturally lower signal.
  - SNR < 10 already flags only 0.9% of the brain. **Recommendation: keep 10.** The criterion
    is nearly inert on this subject.
- **Bit 9 (r2pEcho):** **the p-value form is unusable.**
  - F is large everywhere (median 25 in the cleanest tissue). R2′ from echo 1 is consistently
    about 30% lower than from echo 2: a systematic regime/model effect (Phase 4), not a
    voxel defect.
  - A p-test also penalizes high-SNR tissue the most.
  - The calibrated equivalent (F > 156) flags about 0%.
  - **Recommendation:** disable it, or redefine it as an effect size (for example the relative
    echo difference against the clean-WM p99 of 0.88). **Q3.**
- **Bit 5 (vein):** a fixed 0.15 is inappropriate here, because vCBV_int runs high (clean-WM
  median 0.08–0.10). p99 of clean WM gives 0.59, which only catches extreme voxels; p95 would
  be 0.16. **Q5.**
- **Bit 1 (quadSign):** flags 5% of clean WM and 25–36% of WM and brain, because k is noisy
  (SE_quad is about 8× SE_int).
  - It is the largest *sole cause* of false exclusion in clean WM.
  - **Recommendation:** disable, or flag only k significantly ≥ 0 (k/SE_k > 2). **Q6.**

**Combined false exclusion in clean deep WM** (target < 10%), with `dBCorr`:

| set | clean deep WM | WM | GM | deep GM | brain |
|---|---|---|---|---|---|
| A: calibrated z/ΔGz, SNR 10, r2pEcho off, vein 0.15 | 12.3% | 47% | 70% | 68% | 68% |
| B: A + vein at p99 (0.594) | 7.8% | 38% | 59% | 41% | 60% |
| **C: B + quadSign off** | **3.3%** | **22%** | **44%** | **32%** | **43%** |

Without correction: A 18.1%, B 6.3%, C 2.8%.

**Set C** meets the target. Its per-bit rates (clean deep WM / brain):

| criterion | clean deep WM | brain |
|---|---|---|
| implausible | 1.0% | 19.7% |
| vcbvZ | 1.0% | 5.8% |
| lowSNR | 0% | 0.9% |
| csf | 0% | 18.6% |
| vein | 1.0% | 4.1% |
| zshimEdge | 0.8% | 15.8% |
| zshimEcho | 1.0% | 14.2% |
| zshimFrac | 0.6% | 11.4% |

Set C config:
```json
"quality": {"enabled": true, "useQuadSign": false, "useR2PEcho": false,
            "vcbvZ_max": 8.9, "snr_min": 10, "zshimEcho_maxDiff": 12.7, "vein_vCBVmax": 0.59}
```
with `dBCorr: true`.

### 9.3 Per-ROI selection bias

**Method:** each criterion on its own, with calibrated thresholds, against the reference with
`dBCorr` on and every criterion off. The numbers are the median across regions [IQR] of the
change in the ROI-fit values. "lost" = regions that go from finite OEF to NaN.

| run | lost | median excluded | ΔOEF | ΔR2′ | ΔvCBV |
|---|---|---|---|---|---|
| *dBCorr vs baseline* | 8 | 0 | +0.012 [+0.002, +0.047] | −2.22 [−5.2, −1.0] | −0.046 [−0.107, −0.018] |
| implausible | 0 | 0.17 | −0.031 [−0.111, −0.013] | +0.19 | +0.020 |
| quadSign | 10 | 0.29 | +0.014 | +0.15 | −0.003 |
| vcbvZ | 0 | 0.05 | −0.011 | +0.10 | +0.007 |
| lowSNR | 0 | 0.00 | 0.000 | 0.00 | 0.000 |
| csf | 1 | 0.20 | −0.010 | −0.28 | −0.001 |
| vein | 14 | 0.01 | 0.000 | −0.01 | 0.000 |
| zshimEdge | 2 | 0.12 | −0.009 | −0.05 | +0.001 |
| zshimEcho | 0 | 0.11 | −0.014 | −0.14 | +0.002 |
| zshimFrac | 0 | 0.08 | −0.011 | −0.02 | +0.004 |
| **r2pEcho (placeholder)** | **26** | 0.63 | **+0.093** [+0.048, +0.186] | −0.80 | −0.039 |
| **set C (combined)** | 7 | 0.44 | **−0.031** [−0.108, −0.014] | −0.30 | +0.010 |

**What the table shows:**
- **r2pEcho** has the largest selection bias, roughly +0.09 OEF. That is another reason to keep
  it off.
- **"implausible"** has the next largest, −0.03 OEF. It removes voxels with failed or
  out-of-range fits, and those are not random with respect to OEF.
- The z-shim and CSF criteria each move OEF by about 0.01, and lowSNR does nothing.
- **"vein"** loses 14 regions while excluding only about 1% of occupancy. Removing a small
  region's few high-vCBV voxels can push its ROI-fit vCBV to ≤ 0, which makes OEF NaN. Worth a
  look. **Q5.**
- **Set C** gives OEF for 147 of 160 regions, against 137 with dBCorr alone and 145 at
  baseline. Median OEF across regions: baseline 0.268, dBCorr 0.276, set C 0.231.

**Example ROIs** (baseline / dBCorr / set C):

| ROI | OEF | R2′ | vCBV | excluded (set C) |
|---|---|---|---|---|
| L UnsegmentedWM | 0.308 / 0.377 / 0.256 | 5.92 / 4.33 / 5.71 | 0.089 / 0.053 / 0.103 | 31% |
| L Thalamus | 0.227 / 0.258 / 0.214 | 9.01 / 5.92 / 5.71 | 0.184 / 0.106 / 0.124 | 34% |
| L Caudate | 0.370 / **0.686** / 0.290 | 15.1 / 10.2 / 11.2 | 0.189 / 0.069 / 0.179 | 52% |
| ctx-lh-precuneus | 0.270 / 0.273 / 0.256 | 6.04 / 5.66 / 5.34 | 0.104 / 0.096 / 0.097 | 31% |
| ctx-lh-lateralorbitofrontal | **0.748** / 0.756 / 0.206 | 26.8 / 19.3 / 14.1 | 0.166 / 0.118 / 0.317 | **94%** |
| ctx-lh-inferiortemporal | NaN / NaN / 0.186 | 16.3 / 6.8 / 7.2 | 0.000 / 0.000 / 0.181 | **91%** |
| L Hippocampus | 0.262 / 0.303 / 0.199 | 19.1 / 9.5 / 5.5 | 0.337 / 0.145 / 0.129 | 79% |

Regions over the sinuses and petrous bone (orbitofrontal, inferior temporal, hippocampus) lose
80–94% of their voxels. Their implausible values (OEF 0.75) become plausible, but they then rest
on very few voxels. Consider a stricter `minVoxKept`, or reporting `fracExcluded` alongside.

### 9.4 Spatial check

I rendered a montage of 8 axial slices with set C flags coloured by category over the echo-1
S0, and looked at it. It was then deleted.

**What it showed:**
- **Z-shim flags (red)** sit in the inferior temporal lobes and the orbitofrontal / gyrus
  rectus region. Those are the classic air–tissue susceptibility areas next to the mastoids and
  the frontal sinus.
- **CSF flags (blue)** trace the ventricles and the cortical surface.
- **Vein flags (yellow)** are sparse: lateral cerebellum (near the transverse sinus) and one
  frontal blob near the frontal horn.
- **Fit/consistency flags (green)** are scattered along the cortical ribbon.
- **Deep WM is largely clean.**
- **Unexpected:** a bilateral paracentral **red (z-shim) cluster in the most superior slice**.
  It could be a slab-edge or slice-profile effect, or the superior sagittal sinus. Worth looking
  at in the viewer.

**Clustering numbers:**
- 61% of voxels next to CSF or outside the brain are flagged, against 26% of interior voxels.
- In-plane neighbours of flagged voxels are 73% flagged, against 43% expected at random and 16%
  around kept voxels. The flags cluster; they are not random.
- By slice (inferior → superior), the flagged fraction is 93, 84, 72, 53, … 25–39 in the
  mid-brain … then 54, 68, 82, 93. It is highest at both slab ends.
- Category shares of flagged voxels: CSF 44%, z-shim 37%, fit/consistency 16%, vein 4%.

**A visual check in `ASE_Viewer` is still needed.** I drove it in batch mode (Phase 8), but
only you can judge the overlay by eye.

### 9.5 Summary and recommendation

**Do not turn anything on by default yet.** In order of priority:

1. **The plumbing works:**
   - Default outputs are bit-identical to baseline at every phase.
   - All diagnostics have unit tests (14, all passing).
   - Runtime is unchanged (about 12–16 s per subject).
   - Peak memory goes from 5.1 GB to about 7.0 GB with everything on.
2. **Z-shim correction (`dBCorr`) is the most consequential change.**
   - It brings R2′ down by about 2.2 s⁻¹ and vCBV by about 0.05 per ROI, with OEF about the same
     (+0.01).
   - It makes vCBV_quad agree better with vCBV_int, and halves the R2′ echo F statistic.
   - Cost: it loses fits in about 6% of brain voxels (Q1), and has the ROI-curve subset issue
     (Q2).
   - It deserves its own decision before the mask.
3. **For the mask, set C is a reasonable starting point:**
   - Calibrated z and ΔGz thresholds, SNR at 10, vein at clean-WM p99.
   - quadSign and r2pEcho off.
   - Clean deep WM false exclusion is 3.3%.
   - It still removes about 43% of brain voxels, mostly CSF borders and susceptibility regions,
     which is what it should do.
   - **Validate on more subjects** before trusting the thresholds.
4. **Redesign, don't just recalibrate, bits 1 and 9** (Q3, Q6).

---

## Open questions

Each of these is a decision I made in your absence (or deliberately left alone). All are easy to
change.

- **Q1 — Fits lost after dB correction.**
  - When F < Fmin wipes out every echo-2 long frame in a voxel, the long design loses its TE
    contrast and the fit becomes NaN. That is about 6% of brain voxels, which then get bit 0.
  - **Options:**
    - (a) Keep this (the current behaviour).
    - (b) In those voxels, fall back to a single-echo long design.
    - (c) Divide by F but don't NaN frames, relying on bit 8 instead.
- **Q2 — ROI curves with voxel × frame NaNs.**
  - `regionMeans` averages each frame over whichever voxels are still finite, so long-τ frames
    come from a lower-Gz subset than short-τ frames.
  - **Options:** drop from the ROI any voxel that has a NaN'd used frame (stricter), or keep the
    per-frame NaNs only for the voxel fits.
- **Q3 — Bit 9 (R2′ echo consistency).**
  - The p-value test flags 64–80% of the brain, because echo-1 R2′ is consistently about 30%
    lower: the slope is still in the transition region (echo-1 long frames start at 1.2 ×
    t_crit).
  - **Options:** disable it, switch to an effect size (relative difference > clean-WM p99), or
    start echo 1's long range later.
- **Q4 — `anat_fs_xvent_mask_al` semantics.** It is brain-minus-ventricles, so I used
  CSF = 1 − mask. Should we use a real CSF partial-volume map instead (for example aseg CSF
  labels plus sulcal CSF)? Bit 4 currently catches ventricles and the outer brain edge, but not
  sulci. *(Premise corrected 2026-11-16: xvent does include sulcal CSF; see Decisions, Q4.)*
- **Q5 — Vein threshold.**
  - A fixed 0.15 doesn't fit this subject, whose vCBV_int is high overall.
  - Should it be relative (a percentile of clean WM, or a multiple of the regional median)?
  - Also, the vein bit alone loses 14 small regions: their ROI vCBV becomes ≤ 0.
- **Q6 — Bit 1 (k ≥ 0).** k is noisy, so this flags 25–36% of WM and brain. Should it require
  significance (k/SE_k > 2, from the SEs computed in Phase 4)?
- **Q7 — Smaller choices:**
  - σ uses n−1 degrees of freedom.
  - SNR baseline = dTE = 0 frames.
  - Bit 3 uses the minimum over echoes.
  - Noise weighting uses echo-1 SNR².
  - The dB correction uses each echo's own Gz rather than an echo-averaged Gz.
  - The bit 7 placeholder was 20 µT/m.
  - `qualityMaskNames` uses `segSuffix` for its output name.
- **Q8 — Housekeeping:**
  - The plan's reasoning documents aren't in `testData/plans/`; only the action plan is.
  - `testData/scratch/bch046_250816/` predates this session. I left it untouched. Delete it?

## Final checks (2026-11-12 ~21:35)

- **Branch:** `quality-mask`, one commit per phase (Phase 0–9). Nothing merged, pushed or
  rebased. `main` is untouched.
- **Original data:** unchanged. The size+mtime fingerprint of all 1,104 files under
  `testData/bch046_250816`, `datalist_ase`, `GM_WM_long.txt`, the old failure log and `plans/`
  matches the pre-run snapshot. The md5 of `session_ase/*.json` and `*.gz` matches too.
- **`testData/scratch/` contents:**
  - **`qm/`:** the scratch copy of the subject, which holds the regenerated header and the
    resampled `anat_fs_xvent_mask_al_ase.nii.gz`. Its `fits/` holds **only the Phase 0 baseline**
    (`20261112_205823`), plus `phase0_run.log`.
  - **`bch046_250816/`:** the pre-existing folder (Q8).
  - Per the guardrails, the baseline stays until you sign off. Delete `testData/scratch/qm/`
    when you're done.
- **Dropbox note:**
  - `testData/` is inside Dropbox. Once, a deleted scratch fit folder reappeared, presumably
    restored by sync; I deleted it again.
  - The Phase 9 cleanup happened about 21:34.
  - If anything reappears under `testData/scratch/qm/bch046_250816/session_ase/fits/` besides
    `20261112_205823`, it's safe to delete.
- **Tools left outside the repo** (not committed):
  - `~/.claude/overnight_qm/` holds the diagnostic MATLAB scripts (`diag/`), the original-data
    fingerprint and the restart watchdog's log. Safe to delete.
  - The launchd restart agent was **unloaded and removed** when the work finished (about 21:37),
    so no further automatic restarts will happen.

## Decisions (review session, 2026-11-16)

### Q1 — F² weighting of z-shim-corrected frames (decided: implement)

**Change:**
- With `dBCorr`, each corrected voxel-frame gets weight F² in every log-space fit. Dividing by
  F scales the noise by 1/F, so F² is the matching inverse-variance weight.
- Frames are set to NaN only below `zshim_Fcut` (0.1), where the sinc model itself is
  unreliable. `zshim_Fmin` (0.5) now only defines "heavily attenuated" for bit 8.
- `zshim_weightF2=false` restores the old behaviour (unweighted, NaN below Fmin).
- **Where the weights apply:**
  - `fitMultiRegressor` takes an optional weight array. Unit weights are bit-identical to OLS,
    and weighted fits match `lscov`.
  - The weights flow through `fitASERegimes`, `fitR2PerEcho` and `vcbvConsistency`. The SEs use
    the weighted projection and the measured (attenuated) signal, σ/(F·S).
- **ROI fit:** each region-mean frame gets weight Σm²/Σ(m²/F²) over the voxels valid in that
  frame.
- **New config keys:** `zshim_weightF2` (true) and `zshim_Fcut` (0.1). The config block now has
  28 keys.

**Tests:**
- New `Tests/test_fitMultiRegressor.m`.
- New Monte Carlo test with attenuated frames:
  - The analytic SE is 0.01452, against an MC SD of 0.01451.
  - The weighted fit has SD 0.0145, against 0.0181 unweighted.
- All 17 tests pass.

**Scratch subject** (brain = xvent tissue fraction > 0.9):
- Default outputs are still bit-identical to the Phase 0 baseline.
- **Voxels that lost their fit:** under legacy dBCorr this is 5.5% of brain. All of them get
  R2′ and vCBV back with F² weighting.
- **Recovered voxels against their neighbours:**
  - R2′ is 21 vs 16 s⁻¹ and vCBV 0.32 vs 0.25.
  - OEF is 0.22 vs 0.20.
  - vCBV SE is about 2× larger.
  - These voxels sit in susceptibility regions, and bit 8 still flags them.
- **Voxels with no frame below F = 0.5** are essentially unchanged: median |ΔOEF| is 0.0004.
- Brain medians are unchanged: R2′ 5.03 → 5.11, vCBV 0.101 → 0.101, OEF 0.221 → 0.222.
- **ROI stats:**
  - Median across regions, F² minus legacy: ΔOEF +0.001 [IQR −0.009, +0.013].
  - Finite-OEF regions drop from 137 to 132. High-Gz voxels now re-enter the region means at
    full membership, which is the Q2 question.
  - L Caudate OEF is 0.370 at baseline, 0.686 with legacy dBCorr and 0.612 with F².

### Q2 — ROI curves: one fixed weight per voxel (decided: option c)

**Change:**
- With `dBCorr` and F² weighting, each voxel's ROI membership is multiplied by one weight, used
  the same way in every frame.
- The weight is F² at the voxel's worst used frame (long or short regime).
- If any used frame was dropped (F < `zshim_Fcut`), the weight is 0. That keeps the voxel set
  identical across frames.
- Every frame of a region curve therefore averages the same voxels in the same proportions.
- The per-frame ROI fit weights from Q1 stay.
- There is no new config key: this comes with `zshim_weightF2`.

**Scratch subject:**
- Default outputs are still bit-identical to the Phase 0 baseline.
- **Weights:** the total ROI weight is 69% of occupancy, and 13.6% of occupancy has weight 0.
- **Finite-OEF regions:** 145 at baseline, 137 with legacy dBCorr, 132 with Q1 alone, and
  **153** with Q1 plus Q2(c).
  - Against legacy, 17 regions gain an OEF, mostly orbitofrontal, inferior temporal, entorhinal
    and the corpus callosum.
  - One region loses its OEF (L Amygdala).
- **Clipped regions:** OEF ≥ 1 in 7 regions with legacy, 2 with (c). OEF > 0.6 in 13 regions
  with legacy, 6 with (c).
- **Median across regions:** OEF is 0.276 with legacy and 0.257 with (c), and R2′ is 5.94 for
  both. vCBV is 0.091 with legacy and 0.102 with (c).

| ROI | OEF base / legacy / c | R2′ base / legacy / c | vCBV base / legacy / c |
|---|---|---|---|
| L Caudate | 0.370 / 0.686 / **0.395** | 15.1 / 10.2 / 9.0 | 0.189 / 0.069 / 0.105 |
| L Thalamus | 0.227 / 0.258 / 0.233 | 9.01 / 5.92 / 5.76 | 0.184 / 0.106 / 0.114 |
| L UnsegmentedWM | 0.308 / 0.377 / 0.307 | 5.92 / 4.33 / 4.19 | 0.089 / 0.053 / 0.063 |
| ctx-lh-precuneus | 0.270 / 0.273 / 0.272 | 6.04 / 5.66 / 5.59 | 0.104 / 0.096 / 0.095 |
| ctx-lh-lateralorbitofrontal | 0.748 / 0.756 / 0.319 | 26.8 / 19.3 / 19.1 | 0.166 / 0.118 / 0.277 |
| ctx-lh-inferiortemporal | NaN / NaN / 1.000 | 16.3 / 6.8 / 8.0 | 0.000 / 0.000 / 0.022 |
| L Hippocampus | 0.262 / 0.303 / 0.268 | 19.1 / 9.5 / 8.3 | 0.337 / 0.145 / 0.143 |

**Interpretation:** the legacy L Caudate value of 0.686 came from the curve-shape bias that (c)
removes. Legacy dBCorr pulled the region's vCBV down to 0.069.

### Q3 — Bit 9 off (a), and the regime-boundary investigation (c)

**(a)** `useR2PEcho` now defaults to false. The per-echo R2′ maps and the F statistic are still
written as diagnostics.

**(c) Code change:** `tCrit_buffer` accepts a scalar (symmetric, unchanged default 3 ms) or
`[below above]`, e.g. `--fit 'tCrit_buffer=[3,12]'`. The 1.5·t_c boundary itself is not touched.

**Model simulation:**
- Noise-free Yablonskiy–Haacke cylinder signal, ln S = −ζ·f_s(δω·τ) − R2·TE, on this
  protocol's τ/TE sampling. It goes through the pipeline's regime split and designs.
- **The lower buffer (short regime) has no effect:** −1, −3 and −6 ms give identical results.
  All the bias comes from long-regime frames that are still in the transition.
- **The echo gap appears only when the true t_c is longer than the assumed one.**
  - The assumed t_c is 11.57 ms, from StO2_max = 0.6, i.e. (1−Y) ≥ 0.4.
  - At OEF 0.35–0.45 the gap is under 3%.
  - At lower OEF the gap and the biases grow:

| true OEF (t_c) | buffer +3: echo gap / ΔvCBV / ΔOEF | buffer +12 | buffer +15 |
|---|---|---|---|
| 0.25 (18.5 ms) | 16% / −13% / +13% | 9% / −6% / +6% | 7% / −4% / +5% |
| 0.20 (23.1 ms) | 24% / −27% / +27% | 16% / −19% / +17% | 14% / −17% / +15% |
| 0.18 (25.7 ms) | 27% / −35% / +36% | 19% / −27% / +25% | 17% / −25% / +23% |

- R2 is essentially unbiased (< 1%), and R2′ is biased by only 2–12%. The bias goes mostly into
  vCBV, through the long-regime intercept, and from there into OEF.

**Real data** (dBCorr with Q1/Q2, lower buffer 3, brain = xvent fraction > 0.9):

| buffer above | echo-1 long frames/rep | R2′ e1 / e2 | echo gap (median) | F median | vCBV | OEF | SE(vCBV) | ROI: nOEF / OEF / vCBV |
|---|---|---|---|---|---|---|---|---|
| +3 (current) | 24 | 3.71 / 5.24 | 28.5% | 21.7 | 0.101 | 0.222 | 0.0039 | 153 / 0.257 / 0.102 |
| +6 | 21 | 3.88 / 5.24 | 25.4% | 13.3 | 0.106 | 0.214 | 0.0042 | 154 / 0.249 / 0.107 |
| +9 | 18 | 4.14 / 5.24 | 20.9% | 7.3 | 0.110 | 0.207 | 0.0045 | 154 / 0.242 / 0.110 |
| +12 | 15 | 4.28 / 5.24 | 18.5% | 4.1 | 0.112 | 0.203 | 0.0047 | 154 / 0.238 / 0.112 |
| +15 | 12 | 4.36 / 5.29 | 17.6% | 2.0 | 0.116 | 0.197 | 0.0051 | 154 / 0.232 / 0.116 |

- **The echo gap does not depend on |Gz| below 25 µT/m.**
  - It is 25.5–28.4% at +3 and 16.2–17.8% at +12.
  - It is higher only above 25 µT/m (40% and 28%).
  - So incomplete z-shim correction isn't the cause.
- **The data follow the model at OEF ≈ 0.18–0.20, which agrees with the measured OEF.**
  - Gap in low-Gz voxels: 25.5% → 16.2%.
  - Model at OEF 0.20: 23.8% → 16.0%.
  - The echo-1/echo-2 R2′ discrepancy is therefore what static-dephasing theory predicts when the
    long regime starts before the asymptote.
- **Interpretation:**
  - The 1.5·t_c criterion is right, but t_c is computed for (1−Y) = 0.4. This tissue has
    (1−Y) ≈ 0.2, so t_c ≈ 23 ms and 1.5·t_c ≈ 35 ms.
  - Echo 1's long frames (21–44 ms) are mostly before that.
  - With the current boundary the model implies OEF is overestimated by about 25% and vCBV
    underestimated by the same amount in this subject.
  - The buffer can reduce this but not remove it, because echo 1 stops at 44 ms.
  - One subject, idealized model (cylinders, static dephasing, no diffusion).

**Decision (option 1): `tCrit_buffer` default is now `[3, 12]`.**
- It is set in `BashFiles/ase_config.json` and as the `ASE_fitting` default, so the two agree.
- **This is the first change to default outputs on this branch.** Default runs no longer match
  the Phase 0 baseline.
- Scratch subject, default run (no dBCorr), against baseline:
  - Finite-OEF regions: 145 → 143.
  - Median across regions: OEF 0.268 → 0.249, R2′ 8.98 → 9.23, vCBV 0.125 → 0.140.

**Simulation script:** `Scripts/regimeBoundarySim.m`, parameterized: OEF, ζ, buffers,
StO2_max and protocol. It reproduces the tables above.

**StO2_max** (only used to compute t_crit; noise-free simulation, ζ = 0.06):

| StO2_max (1.5·t_c) | buffer | short / echo-1 long frames | ΔOEF at true OEF 0.20 / 0.30 / 0.40 / 0.50 |
|---|---|---|---|
| 0.6 (17.4 ms) | [3,3] | 15 / 24 | +27 / +7 / +5 / +6 % |
| 0.6 (17.4 ms) | [3,12] | 15 / 15 | +17 / +2 / +6 / +7 % |
| 0.7 (23.1 ms) | [3,3] | 21 / 18 | +20 / +3 / +4 / +5 % |
| 0.8 (34.7 ms) | [3,3] | 32 / 7 | +10 / 0 / +4 / 0 % |
| 0.8 | [3,12] | 32 / 0 | no fit (no echo-1 long frames, R2 not identifiable) |

What the table shows:
- **Raising StO2_max reduces the low-OEF bias.**
- **The longer short regime costs little in the model.**
- **The real cost is precision:** fewer echo-1 long frames weaken the R2 identification. This
  simulation has no noise, so it doesn't show that cost. It needs a Monte Carlo before any
  change.
- **A floor of about 5% at high OEF remains** in every setting. There, the short frames
  (≤ 14 ms) run past the true t_c (9–12 ms).

**StO2_max with noise** (`regimeBoundarySim(..., 'nMC',5000)`):
- Gaussian noise; each τ frame acquired 2 times.
- **SNR:** echo-1 SNR at τ = 0 is 100 (typical tissue in this subject) or 50.
- **Columns:** voxel-level OEF bias / robust SD / RMSE, all in % of truth.
- **The model checks out:** at SNR 1e5 it reproduces the noise-free numbers.

SNR 100, ζ = 0.04:

| setting | true OEF 0.20 | 0.30 | 0.40 | vCBV SD @0.20 |
|---|---|---|---|---|
| 0.6 [3,3] (old) | +27 / 15 / 31 | +7 / 11 / 13 | +5 / 12 / 13 | 12 |
| **0.6 [3,12] (default)** | +18 / 15 / 23 | +3 / 11 / 12 | +5 / 14 / 15 | 14 |
| 0.7 [3,3] | +20 / 13 / 24 | +3 / 10 / 11 | +4 / 12 / 13 | 13 |
| 0.7 [3,6] | +18 / 14 / 22 | +2 / 11 / 11 | +5 / 13 / 14 | 14 |
| 0.8 [3,3] | +11 / 14 / 18 | 0 / 13 / 13 | +5 / 18 / 19 | 18 |

At SNR 50 and at ζ = 0.06, the ranking is the same, and the SD scales roughly as 1/SNR and 1/ζ.
Failed fits are under 1% everywhere.

**Conclusion:** keep StO2_max = 0.6, and flag it for revisiting. A comment next to the option
in `ASE_fitting.m` marks it.
- **At voxel level, no setting wins clearly.**
  - StO2_max = 0.7 is about equivalent to the new buffer default.
  - 0.8 trades about 9 points of bias at OEF 0.2 for about 50% more SD at OEF 0.4.
  - Net RMSE differences are a few points.
- **At ROI level the picture is different.** Averaging over hundreds of voxels removes most of
  the noise, so bias dominates. There, 0.8 [3,3] would cut the low-OEF bias from 18% to 11% at
  no cost in accuracy at high OEF.
- **What would justify revisiting:**
  - A cohort with systematically low OEF or low Hct.
  - Group comparisons across OEF levels. The bias depends on OEF, so it compresses
    differences: true 0.20 vs 0.40 reads as about 0.235 vs 0.42 with the current default.
- **No setting removes the bias.** About 10% at low OEF and 4–5% at high OEF remains with every
  regime split. The real fix is fitting the full static-dephasing model.

### Q4 — CSF mask (decided: keep xvent; SynthSeg posteriors available as an opt-in)

**Correction:** `anat_fs_xvent_mask_al` = aseg-labelled brain minus the ventricles. Checked at
1 mm on the scratch subject:
- It is 1 in 100% of aseg cortex voxels.
- It is 0 in all 140,722 voxels that are inside the brain mask but unlabelled by aseg, which are
  mostly sulcal CSF.
- So CSF = 1 − mask already includes sulcal CSF as FreeSurfer defines it. The overnight claim
  that it doesn't was wrong.

**Alternative tested:** SynthSeg posteriors, `session_afni_nortc/anat_fs_al_posteriors.nii.gz`.
- **Format:** 33 maps (background + 32 labels). Labels were matched against
  `anat_fs_al_vol.csv`, e.g. map 17 = "csf" (225,770 mm³).
- **Tissue fraction** = 1 − (background + CSF + all ventricles).
- **In the same aseg-unlabelled voxels:** the posterior tissue fraction averages 0.50.
  SynthSeg's cortex extends into FreeSurfer's sulcal CSF.

| | xvent | SynthSeg posterior |
|---|---|---|
| cortex with CSF fraction > 0.2 (ASE grid) | 14.9% | 9.0% |
| WM / deep GM / ventricles | 0.5% / 4.4% / 97.7% | 0.5% / 5.1% / 94.7% |
| median CSF exclusion per cortical ROI (set C + dBCorr) | 26% | 17% |
| ROI OEF, posterior − xvent | — | +0.002 cortex, +0.001 WM parcels |
| regions with OEF | 145 | 145 |

**Decision:** keep xvent as the default. It is consistent with the FreeSurfer segmentation that
defines the ROIs and needs no extra input. ROI values barely depend on the choice.

**Opt-in SynthSeg stage** (`asePreProcess.sh` stage 7b, off by default):
- **Enable it with:**
  - `--set tissuePosteriorName=anat_fs_al_posteriors.nii.gz`
  - `--set qualityMaskNames=anat_fs_tissueFrac_al.nii.gz`
  - `qualityMaskNames` set the same way for `aseFit.sh`.
- **Keys:** `tissueFracName`, `posteriorNonTissueVols` (0-based maps, default
  `0,3,4,11,12,16,21,22`) and `posteriorNVols` (33; checked).
- **Behaviour:**
  - It runs from the subject folder, because AFNI rejects paths with spaces.
  - It is rebuilt only when missing or on `--regenTransfer`.
- **Tested on scratch:** the default run is a no-op, and the opt-in run builds and resamples the
  fraction.

### vCBV scale vs An & Lin 2003 (discussion during Q5)

**Question:** pipeline vCBV (9–16%) is about 3× the published ASE values (An & Lin 2003, MRM
50:708). Is the difference in acquisition or in analysis?

**The paper**, read from `../literature/An_2003.pdf`:
- **Acquisition:** 3T Allegra, 3 mm isotropic, single echo (TE 62 ms), 3×3 low-pass filter.
- **Flow suppression:** diffusion gradients, b = 13–120 s/mm² (Venc 0.3–1 cm/s).
- **Fit:** linear long-regime fit over dephasing times 20–40 ms. vCBV = extrapolated line at 0
  minus ln(mean spin echo).
- **ROIs:** hand-drawn GM/WM on 3 slices, artifact regions excluded.
- **Values from the figures** (approximate):
  - vCBV with flow suppression: WM 3.0–3.4%, GM 4.0–4.2%. Without suppression (b = 0):
    3.6% / 5.0%.
  - OEF 35–40%.
  - R2′: WM 2.6, GM 3.5 s⁻¹.
- **So:** by the paper's own numbers, intravascular signal adds only about 0.5–1 point of vCBV.

**Refit of this subject with the paper's design:** `Scripts/anLin2003DesignFit.m`. Masks are
clean deep WM, WM and cortex, each with tissue fraction > 0.9.

| analysis | vCBV cWM / WM / ctx | R2′ cWM / WM / ctx | OEF cWM / WM / ctx |
|---|---|---|---|
| paper design (echo 1, 20–40 ms), raw | 4.7 / 5.8 / 8.0% | 3.4 / 4.0 / 5.7 | 0.34 / 0.32 / 0.30 |
| + z-shim correction (F², as dBCorr) | 4.1 / 4.5 / 6.0% | 3.0 / 3.1 / 4.2 | 0.34 / 0.31 / 0.29 |
| paper design, window 30–44 ms, raw | 7.4 / 9.1 / 11.8% | 4.1 / 5.0 / 6.9 | 0.26 / 0.25 / 0.24 |
| pipeline defaults ([3,12], 2 echoes), no dBCorr | 11.1 / 13.3 / 16.1% | 5.2 / 6.2 / 8.5 | 0.22 / 0.21 / 0.21 |
| pipeline defaults, dBCorr | 8.8 / 10.0 / 11.3% | 4.3 / 4.7 / 5.8 | 0.23 / 0.21 / 0.21 |

The 3×3 smoothing changes nothing material.

**Conclusions:**
- **With the paper's analysis, this dataset gives values close to the paper's.** The roughly
  3× gap comes mostly from the analysis, not the acquisition.
- **Correction of an earlier inference:** AnLab's current code agreeing with the pipeline meant
  the data were responsible. It doesn't, because that code uses a multi-echo design like ours,
  not the 2003 single-echo fit.
- **The dephasing-time window drives absolute vCBV and OEF.**
  - vCBV goes about 4% (20–40 ms) → 6–7% (30–44 ms) → 9–11% (pipeline, out to 76 ms, both
    echoes).
  - OEF goes 0.34 → 0.25 → 0.22.
  - ln S keeps steepening across 20–76 ms. This is the same effect as the Q3 echo gap.
- **Which window is right isn't decidable within the regime-split method:**
  - Static-dephasing theory (`regimeBoundarySim`) says early windows underestimate vCBV and
    overestimate OEF when the true t_c is long.
  - The continued steepening into echo 2 suggests an extra decay that the model doesn't include
    (macroscopic in-voxel gradients, diffusion), and that would bias late windows the other way.
  - This strengthens the case for a full-model fit (open project).
- **Acquisition notes for future scans** are real but secondary:
  - No flow suppression adds about +15–25% vCBV, per the paper's Fig. 4.
  - 4.8 mm slices give more through-slice dephasing. With the paper's design, z-shim
    correction lowered vCBV from 4.7% to 4.1%.

### Q5 — Vein threshold (decided: relative, option b)

**The "regions lost to the vein bit" from Phase 9 aren't a threshold problem.**
- With the vein bit enabled alone, susceptibility regions lost their OEF. Those regions are full
  of failed voxel fits: voxel-median vCBV is −0.5 to −0.8 in orbitofrontal and inferior temporal.
- Trimming only the high-vCBV tail leaves a set dominated by negative values, so the region's
  vCBV goes ≤ 0.
- With bit 0 (implausible values) also on, as in any real configuration, both tails are removed.
- So the vein bit should not be used without bit 0.

**Change:**
- New `quality.vein_mode`, default `"relative"`: flag vCBV > `vein_relK` (4) × the median
  vCBV over tissue voxels.
- **Tissue voxels** for that median: csfFrac < 1 − `vein_refTissueMin` (0.9), with all four raw
  maps finite and within `physBd_*`.
- The threshold scales with the subject's vCBV level, and with the analysis choices that set it
  (window, dBCorr).
- **Inputs and failure handling:**
  - Relative mode needs csfFrac, i.e. the quality masks. Without it the bit is reported as
    unavailable and skipped.
  - With fewer than 100 reference voxels it warns.
- **Absolute mode** (`"absolute"`) keeps `vein_vCBVmax`.
- `computeVoxelQuality` returns the reference and threshold (4th output), and verbose runs print
  them.
- **Tests:** 3 new; 20 of 20 pass.

**Scratch subject** (set C criteria + dBCorr, against the same run without the vein bit):
- **Threshold:** 0.48 (median tissue vCBV 0.120 over 51,807 voxels).
- **Flagged:** clean deep WM 2.0%, WM 5.5%, cortex 9.8%, deep GM 9.4%.
- **Vein-excluded fraction per region:** median 3.7%, p90 28%.
- **ROI effect:** OEF median +0.0006 [IQR 0, +0.008], vCBV −0.003, R2′ −0.09.
- **5 regions lose their OEF.** All five fall below `minVoxKept` (5 effective voxels): left
  entorhinal and pars orbitalis cortex, right entorhinal, pars orbitalis and frontal-pole WM.
  These regions were already 95%+ excluded by the z-shim and CSF criteria.
