# Motion correction and ASE contrast: problem and plan

Status: **open, not started** (written 2026-11-20). Nothing in the pipeline has changed because of this page.

## The problem

`asePreProcess.sh` (stage 4, `motionMode = echo1`) runs `3dvolreg -Fourier -twopass` on all 98 echo-1 frames against raw frame 0 (dTE = 0). The other echoes reuse those transforms. The old cluster script (`clusterFiles/batch_preprocess`) did the same for echo 1, so both pipelines share this issue.

3dvolreg finds the rigid transform that makes each frame match the base image by least squares. That assumes each frame differs from the base only by head position. In ASE the contrast changes from frame to frame **by design**:
- the signal decays with τ;
- dropout near the sinuses and skull base grows with τ;
- the z-shim block changes the through-slice dephasing directly.

The sweep is acquired in an interleaved order (dTE = 0, 44, 35, 26, 17, 8, 43, …), so neighbouring frames differ strongly in contrast. 3dvolreg partly reads those contrast changes as head motion.

## Evidence (bch046_250816, echo 1)

**Sweep frames.** After removing slow drift (a quadratic in time), dTE explains this much of each motion parameter's frame-to-frame variation:

| parameter | variance explained by dTE | change from dTE 0 to 44 ms |
|---|---|---|
| dP | 83% | +0.62 mm |
| yaw | 68% | −0.17° |
| dS | 61% | −0.27 mm |
| pitch | 32% | +0.16° |
| roll | 9% | +0.03° |
| dL | 5% | −0.02 mm |

A real head does not move in step with the dTE schedule. The affected directions (anterior–posterior, superior–inferior, pitch) match dropout spreading in the inferior frontal and temporal regions.

**Z-shim block (frames 90–97).** These frames differ only in the applied z-moment, yet the parameters span a 2.5° pitch range and a 1.2 mm dS range. Across the whole 90-frame sweep the ranges are 0.37° and 0.51 mm.

**Cross-echo check.** Before 2026-11-16 echo 2 was registered separately, across TE, against the echo-1 base. That gave a constant offset of about 0.45 mm dS, 0.45° pitch and 0.4 mm dP (and +5% vCBV). This was fixed by `motionMode = echo1`, see [decisions.md](decisions.md). The same mechanism drives the within-echo bias above.

Plot: `docs/figures/preproc_registration/reg5_motion.png` (made by `Scripts/comparePreprocRegistration.m`).

## Reading reg5_motion.png (old vs new preprocessing, bch046)

**Top six panels: motion parameters per frame.** The x-axis is the acquired frame (0–97, scan order); the grey band is the z-shim block (frames 90–97). Each panel is one of 3dvolreg's six rigid-body parameters (dS, dL, dP = shifts toward superior, left, posterior), relative to raw frame 0.

| marker | what it is |
|---|---|
| blue line | new pipeline: one 3dvolreg run over all of echo 1; the same transforms are applied to echo 2 |
| orange circles / squares | old pipeline, echo 1: 3dvolreg run separately on its two dTE-sorted halves (split 1 / split 2) |
| grey circles / squares | old pipeline, echo 2: its own 3dvolreg, registered to the echo-1 base |

1. **Orange on blue: old and new echo-1 motion are identical** (RMS difference ≤ 0.002 ° or mm). 3dvolreg registers each frame to the base on its own, so the old repetition split (or any frame order) does not change a frame's transform.
2. **Grey off blue: the old echo-2 registration was wrong.** Echo 2 comes from the same excitation 32 ms later, so its true motion equals echo 1's. The old echo-2 estimates sit at a near-constant offset (dS ≈ −0.5 mm, pitch ≈ +0.4°, dP ≈ +0.4 mm, roll ≈ +0.1°): an error from matching a TE 93 ms image (more decay and dropout) to a TE 61 ms base. The new pipeline applies the echo-1 transforms to echo 2 instead (decision 2026-11-16, worth +5% vCBV).
3. **The blue 5-frame sawtooth is contrast, not motion.** It follows the interleaved dTE order (table above). Both pipelines share it.
4. **The z-shim block has the largest apparent motion**: up to 1.9° pitch and −1.4 mm dS (blue), and 2.7° and −2.4 mm (grey). Same mechanism: the z-moments change the inferior dropout.
5. **The slow trends are probably the real motion, and they are small**: roll rises to about 0.2° by frame 25 and drifts back, dL steps up about 0.1 mm near frame 17, yaw drifts slowly. For this subject the real motion (under about 0.2° and 0.2 mm) is smaller than the contrast-driven sawtooth, so here motion correction mostly adds contrast-driven shifts. Check whether that holds for subjects who move more before choosing an option below.

**Bottom panel: old vs new motion-corrected images.** For each frame: RMS(old − new) / RMS(new) inside the brain, in %. It is a whole-image summary dominated by edges and high-contrast areas, not a per-voxel signal change.
- **Echo 1: about 0.01%.** The images are effectively identical.
- **Echo 2: a steady ~4.5% across the sweep.** This mixes the transform offset in point 2 with an interpolation difference (old: 3dvolreg Fourier; new: 3dAllineate wsinc5). The offset probably dominates (half a millimetre at edges with 1.7 mm voxels costs several percent), but this figure does not separate the two. Re-applying the old echo-2 transforms with wsinc5 would.
- **A spike to 33% in the z-shim block (frames 90–93)**: the two pipelines' transforms disagree most there (the old echo-2 registration's large spurious moves against the echo-1 ones). Neither is right; it shows these frames are the most sensitive to motion handling. They are used only to estimate Gz, never in the decay fits.

**Conclusion for old vs new:** echo-1 motion correction is unchanged, and the new pipeline removes a ~0.5 mm / 0.4° cross-echo misalignment in echo 2. The dTE-locked sawtooth and the z-shim spikes are present in both pipelines and are the subject of this plan.

## Why it matters even though it is sub-voxel

The errors are small: 0.6 mm against 1.72 × 1.72 × 4.8 mm voxels. But they are **correlated with τ, the variable the fit regresses on**. A random registration error would only add noise. A τ-correlated shift adds bias wherever the image changes steeply in space (CSF borders, cortex, the inferior frontal lobe): the shift creates a fake τ-dependence that the fit takes up as R2′ or short-regime curvature.

**Effect on z-shim.** The z-shim frames get spurious transforms before Gz is fitted from them (see the z-shim discussion and `fitZshimGradient.m`).

## Options

| # | approach | idea | effort |
|---|---|---|---|
| 0 | current | 3dvolreg on all echo-1 frames against frame 0 | — |
| 1 | same-contrast estimates, interpolated in time | Estimate motion only between frames with matching contrast: each repetition-2 frame against the repetition-1 frame with the same dTE, and/or only short-τ frames. Interpolate transforms for the other frames from their neighbours in time (3.8 s per frame; real motion is slow). Z-shim frames take the nearest sweep frame's transform. | low (AFNI + a little scripting) |
| 2 | weighted registration | Add a weight mask that excludes the dropout regions (inferior frontal and temporal) | low; reduces the bias but does not remove it |
| 3 | register to the model | Fit the model, predict each frame at its own τ, register each frame to its prediction, refit; repeat 2–3 times | high; the most principled |
| 4 | no correction within the run | Skip motion correction; flag subjects whose real motion (from option 1) is large | none |

## How to evaluate (bch046 first, then a few subjects with more motion)

For options 0, 1, 2 and 4, compare:
1. **Whether the parameters still follow dTE.** Repeat the table above; a good method should be near 0%.
2. **Fit quality:** the RMS residuals (`RMS_L`, `RMS_S`) over the brain and near CSF borders.
3. **The maps:** R2′, vCBV and OEF, especially in the inferior frontal and temporal lobes and at CSF borders.
4. **Z-shim:** Gz maps and how well the two echoes' Gz agree (bit 7) with each option's z-shim transforms.
5. **Between repetitions:** with real motion removed, frames with the same dTE in the two repetitions should match better (the repetition-difference noise σ).

Option 3 only if options 1 and 2 leave a clear bias.
