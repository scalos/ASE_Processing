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
