# Method decisions and open questions

This page lists the method decisions in short form. Full evidence (tables, simulations, figures) is in [quality-mask-report.md](quality-mask-report.md), the detailed lab log, under the section named in the last column. Test subject: `bch046_250816` (healthy volunteer). Every number below comes from this one subject.

## Decisions

| date | decision | setting / default | why | report section |
|---|---|---|---|---|
| 2026-11-12 | Voxel quality mask, noise, z-shim and fit-extension diagnostics added, all optional | `quality.*`, `dBCorr` (off) | The overnight build; default outputs unchanged | Phases 0–9 |
| 2026-11-16 | z-shim-corrected frames are weighted by F² instead of dropped at F < 0.5 | `zshim_weightF2` = true, `zshim_Fcut` = 0.1 | Dividing by F scales the noise by 1/F. Dropping frames lost about 6% of brain fits. | Q1 |
| 2026-11-16 | With dBCorr, each ROI average uses one fixed weight per voxel (F² at its worst frame) | (with `dBCorr`) | Every frame then averages the same voxel mix, so the curve shape isn't biased. Regions with OEF: 137 → 153. | Q2 |
| 2026-11-16 | Per-echo R2′ F-test (bit 9) off | `useR2PEcho` = false | It flagged 64–81% of brain because of a systematic echo-1 / echo-2 slope difference, not bad voxels | Q3 |
| 2026-11-16 | The long regime starts later; the short regime is unchanged | `tCrit_buffer` = [3, 12] | Static-dephasing simulation: the echo gap appears when the true t_c is longer than assumed. Later start: OEF bias at true OEF 0.2 is +27% → +17%. Matches AnLab's 30 ms cutoff. | Q3 |
| 2026-11-16 | StO2_max kept at 0.6, **flagged for revisit** | `StO2_max` = 0.6 | Noise simulation: no clear winner at the voxel level. A higher value would reduce the ROI-level bias at low OEF. | Q3 (StO2_max) |
| 2026-11-16 | CSF fraction from xvent (FreeSurfer); SynthSeg posteriors optional | `qualityMaskNames` = xvent, `tissuePosteriorName` = "" | xvent already treats sulcal CSF as CSF, consistent with the FreeSurfer ROIs. SynthSeg moves ROI OEF by only about 0.002. | Q4 |
| 2026-11-16 | Vein threshold relative to the subject's tissue vCBV | `vein_mode` = relative, `vein_relK` = 4 | The absolute vCBV scale depends on the analysis (window, dBCorr). A fixed 0.15 flagged half of GM. | Q5 |
| 2026-11-16 | Curvature check (bit 1) requires significance | `quadSign_kZmin` = 2 | Plain k ≥ 0 was mostly noise. Significant cases are real shifted-echo artifacts (see open questions). | Q6 |
| 2026-11-16 | Joint (both-echo) Gz fit available but off | `zshim_jointGz` = false | More precise where the linear through-slice model holds, but biased where the echoes disagree (susceptibility regions). Revisit with thinner slices. | Q7 |
| 2026-11-16 | Echo 2 motion-corrected with echo 1's transforms | `motionMode` = echo1 | Both echoes come from the same excitation (32 ms apart). Registering echo 2 separately across TE added a spurious 0.4 mm / 0.4° offset and +5% vCBV. | Follow-ups: echo-2 motion |
| 2026-11-16 | Z-shim zero-moment frame (raw 96) left out of the regime fits | `zshimZeroFrameInFit` = false | It is acquired minutes after the sweep. Removing it shrinks the echo gap by 2 points, with negligible effect on the maps. Matches AnLab. | Follow-ups: frame 96 |
| 2026-11-16 | AnLab's robust regression and smoothing not adopted | — | Small effects; would conflict with the F² weights and the analytic SEs, and blur the quality flags | AnLab fitting steps not adopted |

## Open questions and known limitations

1. **Recalibrate the quality thresholds on several subjects** (next step). Bits 2 (vcbvZ) and 7 (zshimEcho) still have placeholder defaults. The working set-C values predate the November-16 changes.
2. **Absolute vCBV / OEF depend on the fitting window.** On this subject, An & Lin's 2003 design (echo 1, 20–40 ms) gives vCBV about 4–6% and OEF about 0.30. The pipeline gives about 10–11% and 0.21. ln S keeps steepening through 76 ms, so the regime split can't settle which is right. **Candidate project: fit the full static-dephasing model.** Tools: `Scripts/regimeBoundarySim.m`, `Scripts/anLin2003DesignFit.m`.
3. **Shifted-echo artifact.** In susceptibility regions the signal peaks at τ ≈ 17 ms, not at the spin echo, probably because of in-plane field gradients. Bit 1 catches it. A correction (fitting around the shifted peak) would be a separate project.
4. **StO2_max.** Revisit for cohorts with low OEF or low Hct, or for group comparisons: the OEF bias depends on OEF, so it compresses group differences.
5. **Acquisition notes for future scans** (to raise with the PI):
   - There are no flow-crushing gradients. By An & Lin's Fig. 4, that adds about 15–25% to vCBV.
   - The 4.8 mm slices increase through-slice dephasing. Thinner slices may also make the joint Gz fit worthwhile.
6. **ROI summary.** ROI-fit values (region-averaged signal) and voxel medians differ by about 20% in OEF. Decide which one to report.
