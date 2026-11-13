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
