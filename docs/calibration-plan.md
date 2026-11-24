# Calibration plan: quality thresholds and correction defaults across the o2co2 cohort

**Status:** in progress. See [Current status and where to pick up](#current-status-and-where-to-pick-up) below.

**Goal:** turn the voxel quality mask from a one-subject prototype into a default the group analyses can rely on. Everything in [decisions.md](decisions.md) so far comes from one subject (`bch046_250816`), and "set C" was calibrated and evaluated on that same subject, before the November-16 fit changes. This plan recalibrates on the cohort, with subjects held out for validation, and ends with a decision on what to turn on by default.

## Current status and where to pick up

*Updated 2026-11-18, end of day.*

**2026-11-20 to 11-24 (preprocessing review):**
- Cohort preprocessing finished; the only failure is bch045 (the comment script can't find series 116 in `dcm.index.dat`: its ASE comes from a second exam on 08-21 that the index doesn't list). Fixed 2026-11-24: the comment script now reads the whole series folder (see decisions.md).
- **Adopted 2026-11-24: the anatomy comes from the FreeSurfer recon** (`anatSource = fsrecon`, 12-parameter registration). See [anat-route-report.md](anat-route-report.md).
- **Calibration is on hold** until the cohort has been preprocessed again with the new anatomy stages.

**Next, on the cluster (before anything below):**
1. Rerun `Scripts/linkO2co2Ase.sh --src <o2co2_ase>`: it adds the `<subj>/fsrecon` links (existing links are kept).
2. Make FreeSurfer available (`mri_convert` on PATH); `asePreProcess.sh` checks it at startup.
3. Run `asePreProcess.sh` on the cohort. Conversion, comments, header and motion correction are kept; stages 5–8 run for every subject (new anatomy, new registration). Expect a note per subject that the registration settings changed or that the anatomy is new.
4. Read the newest `prepASE_failures_*.txt`. bch045 should now pass, with a note that folder 116 holds DICOM series 16. It never got past stage 2, so the cohort run does stages 2–8 for it with no extra flags. The z-moment warning no longer appears for every subject.

The cluster steps below (Phase 0 QC onwards) follow once this pass is done.

**Done:**
- `o2co2_ase_proc` built on the cluster from links (`Scripts/linkO2co2Ase.sh`). Brainstem masks are missing for every subject; preprocessing and fitting now warn instead of failing.
- AFNI on the cluster switched to 26.0.10 (the 2018 build lacked the `lpa+` cost). `asePreProcess.sh` now prints the AFNI build at startup and stops before any subject if its `3dAllineate` doesn't know the cost in `allineateOpts`.
- Cohort preprocessing restarted with 26.0.10 (running overnight).
- Decisions 1–4 recorded (section 5).
- `Scripts/calibrateQuality.m` (Phase 3) written and run end to end on `bch046_250816`. The two Phase 1 fit commands were checked on that subject as well: both runs complete (about 75 MB each).

**Next, on the cluster:**
1. Confirm preprocessing finished: read the newest `prepASE_failures_*.txt` in `o2co2_ase_proc`. Expect one warning per subject about the missing brainstem source. Anything else needs a look.
2. Check the startup line says `AFNI: /usr/pubsw/packages/AFNI/26.0.10 ...`. If preprocessing was restarted without `--regenNifti`, conversion and motion correction are still from AFNI 18.0.22: decide whether to rerun them (section 6 log).
3. Phase 0 QC: motion and a registration look per subject; fill in the QC table.
4. Phase 1: the two diagnostic runs (`calDiag-dB`, `calDiag-noDB`), then check `fitASE_failures_*.txt`.
5. Phase 2: export with the tar loop and copy to `testData/cohort/` here.

**Next, here, once the export arrives:**
1. Choose the calibration / validation split from the QC table (healthy cohort, balanced by site, group and age).
2. Run `calibrateQuality(rootDir, calibSubjects=[...])` on both runs.
3. Two things the test run showed, to settle before reading the thresholds:
   - **Clean-WM mask size.** The script's clean deep WM on `bch046` has 1,671 voxels; the Phase 9 mask had 3,233. Both erode by a 3×3×3 cube, but the script selects WM by occupancy ≥ 0.9 on the ASE grid, while Phase 9 resampled the aseg labels to the ASE grid and required brain fraction > 0.9. The stricter occupancy cut is the likely cause; not yet checked. Decide on one definition (and whether a through-slice erosion of 4.8 mm slices is wanted).
   - **p99 is fragile to small clusters.** On `bch046`, 23 clean-WM voxels in slices 15–16 have echo-2 Gz fits far from echo 1's (median +27 vs −90 µT/m), so the zshimEcho p99 is 73 µT/m, against a p95 of 5.6. The same subject gave 12.7 in Phase 9 with the larger mask. Look at that cluster in the viewer, and consider a more robust rule (e.g. p99 with the cluster's slices checked, or a multiple of a robust spread).

## 1. Scope

**In scope:**

| item | current state | what we decide |
|---|---|---|
| `dBCorr` (z-shim correction) | off by default; the most consequential change on the test subject | on or off by default, decided before the mask |
| bit 2 vcbvZ, `vcbvZ_max` | placeholder 3 (set C: 8.9, outdated) | threshold |
| bit 7 zshimEcho, `zshimEcho_maxDiff` | placeholder 20 µT/m (set C: 12.7, outdated) | threshold |
| bit 1 quadSign, `quadSign_kZmin` | 2, calibrated on one subject (noise rate 2.3%) | confirm |
| bit 5 vein, `vein_relK` | relative, k = 4, one subject | confirm |
| bit 3 lowSNR, `snr_min` | 10, nearly inert on one subject | confirm or drop |
| bits 4, 6, 8 (csf, zshimEdge, zshimFrac) | fixed rules | confirm the rates are sensible across subjects |
| bit 9 r2pEcho | off (it flags a model effect) | stays off; confirm the echo-1/echo-2 R2′ gap is systematic across subjects |
| `minVoxKept`, reporting `fracExcluded` | 5 effective voxels | value; whether group stats should also filter on `fracExcluded` |
| `quality.enabled` default | off | on or off, with which settings |

**Out of scope here** (open questions 2, 4 and 6 in [decisions.md](decisions.md)): the fitting window and the full static-dephasing model, `StO2_max`, and ROI-fit versus voxel-median reporting. They change the values the thresholds are measured on, so they should either be settled first or wait until after this calibration. The plan assumes they wait: we calibrate with the current fit defaults and rerun the calibration if those defaults change later.

## 2. Method (carried over from Phase 9, extended to the cohort)

**Calibration tissue: clean deep white matter**, per subject, on the ASE grid:
- aseg labels 2 / 41 (cerebral WM), membership from `segToAseMembership` (occupancy ≥ 0.25), eroded with a 3×3×3 cube;
- CSF fraction = 0 and brain fraction > 0.9.

**Rule per criterion:** threshold at the p99 of the criterion in clean deep WM, so each criterion on its own falsely excludes about 1% of tissue that should be kept.

**Combined target:** with every enabled criterion, false exclusion in clean deep WM below 10% for every subject, not just on average.

**What changes from the one-subject calibration:**
1. **Hold-out.** Subjects are split into a calibration set and a validation set. Thresholds come from the calibration set only; every rate reported for the decision comes from the validation set. Also run leave-one-subject-out to see how much a threshold depends on single subjects.
2. **Pooled vs per-subject.** For each criterion, compare the per-subject clean-WM p99 values. If they agree (say within ±20%), use one pooled threshold. If they spread widely, the criterion should be relative to the subject (as the vein criterion already is), e.g. a multiple of the subject's clean-WM median or robust SD.
3. **Healthy cohort.** Calibrate on the healthy cohort only, `group=ctrl` and `mdt` (19 subjects with data: 11 ctrl, 8 mdt), so disease doesn't move the reference. Then check the false-exclusion rate in the remaining subjects (AD `mdh199`, ipv `mdh206` and `mdh209`, covid `bch001`). A criterion that excludes much more in a patient group is either picking up pathology (bad: biases group comparisons) or real artifacts (fine), and we look at those cases in the viewer.
4. **Outcomes beyond false exclusion,** evaluated at the ROI level:
   - selection bias: per-ROI change in OEF, R2′ and vCBV caused by each criterion alone and by the combined set;
   - coverage: number of ROIs with finite OEF per subject, and per-ROI `fracExcluded`;
   - between-subject spread: CV of ROI OEF across the healthy cohort, with and without the mask. A useful mask should not increase it, and ideally reduces it.

## 3. Phases

### Phase 0. Cohort QC (cluster, before fitting)

1. Read the `prepASE_failures_*.txt` log. Subjects expected to need attention:
   - all of them: missing `anat_fs_bstem_al+orig` (now a warning only);
   - `bch045_250823` and `bch046_250816`: the datalist uses the first of two ASE runs (`ase_116`, `ase_022`). The second runs (`ase_124`, `ase_026`) are commented out; they reflect acquisition problems, not test–retest scans, and are not used.
2. Motion: from `ase_e1_mc.1D`, max displacement and max frame-to-frame movement per subject. Flag subjects above about 1 mm / 1° for a closer look, not automatic exclusion.
3. Registration: look at `ase_base_al` over `anat_fs_brain_al` for every subject in `aseView.sh` or AFNI. A bad registration would corrupt the segmentations, CSF fraction and ROI membership, and every calibration number with them.
4. Record the result in a table (subject, group, motion, registration OK, notes) in section 6 of this file.

### Phase 1. Diagnostic fit runs (cluster)

Two runs per subject, identical except for `dBCorr`. Both are analysed in full, and the `dBCorr` default is decided from the two side by side. They write every diagnostic map, keep raw (unclipped) maps, and compute the mask at the current thresholds. The thresholds themselves don't matter here: every threshold can be re-applied afterwards from the diagnostic maps.

```bash
Q='quality={"enabled":true,"repNoise":true,"zshimEstimate":true,"fitExtensions":true}'
bash BashFiles/aseFit.sh --allSubjDir <o2co2_ase_proc> --tag calDiag-dB   --fit dBCorr=true  --fit writeRawMaps=true --fit "$Q"
bash BashFiles/aseFit.sh --allSubjDir <o2co2_ase_proc> --tag calDiag-noDB --fit dBCorr=false --fit writeRawMaps=true --fit "$Q"
```

Check the `fitASE_failures_*.txt` log after each run and resolve any failures before exporting.

### Phase 2. Export (cluster → this machine)

Two parts. Part B is what lets us iterate quickly: with the fit inputs here, the selection-bias runs in Phase 4 run locally in about 15 s per subject instead of a cluster round trip each time.

**A. The two diagnostic runs**, per subject: the latest `*_calDiag-dB` and `*_calDiag-noDB` run folders (all files: maps, diagnostics, `_ROI_Stats.csv`, `_manifest.json`, `config_resolved.json`). About 80 MB per run uncompressed; the NIfTIs compress well.

**B. The fit inputs**, per subject, from `session_ase/`:
- `ase_e1_proc.nii.gz`, `ase_e2_proc.nii.gz`, `ase_hdr.json`, `ase_e1_mc.1D`
- `ase_base.nii.gz`
- `anat_fs_aparc_aseg_al_ase.nii.gz`, `anat_fs_wmparc_al_ase.nii.gz`, `anat_fs_xvent_mask_al_ase.nii.gz`

About 100 MB per subject, so roughly 2.5 GB for the cohort. Plus, once: `GM_WM_long.txt`, `cov_table_ase`, `datalist_ase`, and the preprocessing and fit failure logs.

```bash
cd <o2co2_ase_proc>
out=$HOME/ase_calib_export.tar
files=(cov_table_ase datalist_ase prepASE_failures_*.txt fitASE_failures_*.txt)
for s in $(grep -v '^#' cov_table_ase | awk 'NF {print $1}'); do
    [[ -d $s/session_ase ]] || continue
    for t in calDiag-dB calDiag-noDB; do
        r=$(ls -d $s/session_ase/fits/*_$t 2>/dev/null | sort | tail -n 1)
        [[ -n $r ]] && files+=("$r")
    done
    for f in ase_e1_proc.nii.gz ase_e2_proc.nii.gz ase_hdr.json ase_e1_mc.1D ase_base.nii.gz \
             anat_fs_aparc_aseg_al_ase.nii.gz anat_fs_wmparc_al_ase.nii.gz anat_fs_xvent_mask_al_ase.nii.gz; do
        [[ -e $s/session_ase/$f ]] && files+=("$s/session_ase/$f")
    done
done
tar -cf "$out" "${files[@]}" && ls -lh "$out"
# GM_WM_long.txt lives outside this folder (segLUTs); add it separately.
```

**Where it goes here:** `testData/cohort/` (gitignored), keeping the `<subj>/session_ase/...` layout so `aseFit.sh --allSubjDir testData/cohort` works unchanged. If 3 GB is too much for Dropbox to sync, use a folder outside Dropbox and point to it.

### Phase 3. Threshold calibration (here, from the exported diagnostic maps)

`Scripts/calibrateQuality.m` (written; `help calibrateQuality` for the options) does, for each subject and each run:
1. builds the clean-deep-WM mask (and WM, cortex, deep GM, brain masks for context);
2. extracts each criterion's statistic: |vCBV_z|, |Gz_e1 − Gz_e2|, kZ, min-echo SNR, vCBV relative to tissue median, CSF fraction, edge and frac flags, R2′ F and the relative echo R2′ difference;
3. writes one long table: subject, group, run, tissue, criterion, p50 / p90 / p95 / p99, n voxels.

**Outputs, for the decision:**
- per criterion, the per-subject clean-WM p99 (calibration set), and the pooled threshold;
- the spread of per-subject p99 → pooled or relative (section 2, point 2);
- per-bit and combined false exclusion per subject, on the validation set, for each candidate set;
- the same rates in WM, cortex, deep GM, brain, and per group;
- leave-one-subject-out stability of each threshold;
- everything with and without `dBCorr`.

Thresholds are applied to the saved diagnostic maps here (no refitting), using the same logic as `computeVoxelQuality.m`. To make sure the two agree, the script should call `computeVoxelQuality` itself on the loaded maps rather than reimplement the rules.

### Phase 4. ROI-level evaluation (local fits)

For the leading one or two candidate sets, plus each criterion alone and a reference with every criterion off, refit locally on the exported inputs (stats only, so each run takes seconds). Then compare, across subjects:
- per-ROI ΔOEF / ΔR2′ / ΔvCBV against the reference (selection bias), as in Phase 9.3 but as a cohort median [IQR];
- ROIs lost (finite → NaN) and ROIs gained;
- between-subject CV of ROI OEF in the healthy cohort;
- `fracExcluded` per ROI: which ROIs lose most of their voxels in most subjects (sinus and petrous regions on the test subject). This decides `minVoxKept` and whether group stats should drop ROIs above some `fracExcluded`.

### Phase 5. Visual check (you)

In `aseView.sh`, on 4–6 subjects chosen to cover the groups and the extremes (highest motion, most and least excluded): flags should cluster at CSF borders, sinuses, the petrous bone and large veins, not in deep tissue. Also check the superior paracentral z-shim cluster seen on the test subject: does it appear in every subject (slab edge) or only some?

### Phase 6. Decision and defaults

1. Summarize Phases 3–5 for you, with a recommendation for each item in section 1.
2. After you decide: update `ase_config.json` and `qualityDefaults.m`, add rows to [decisions.md](decisions.md), update [pipeline.md](pipeline.md) and [settings.md](settings.md), and add tests if any criterion's logic changes (for example a relative threshold).
3. On the cluster, rerun the cohort with the chosen defaults and a tag (e.g. `--tag calib1`) for the group analysis.

## 4. Where the work happens

| step | where | who |
|---|---|---|
| Phase 0 QC, Phase 1 runs, export | cluster | you |
| Phase 3 calibration script and analysis | this machine (MATLAB R2025b) | me |
| Phase 4 local refits | this machine | me |
| Phase 5 visual check | either | you |
| Phase 6 decisions | — | you, from my summary |
| final cohort run | cluster | you |

## 5. Decisions (2026-11-18)

1. **Calibration subjects:** the healthy cohort, ctrl + mdt. Split into calibration and validation halves after the Phase 0 QC, balanced by site (bch / mdh), group and age, with leave-one-subject-out on top (section 2).
2. **Targets:** about 1% false exclusion per criterion (clean-WM p99) and < 10% combined in clean deep WM, for every validation subject.
3. **`dBCorr`:** both runs are calibrated and evaluated in full; the default is chosen from the two results together.
4. **Second ASE runs** of `bch045` / `bch046`: not used (acquisition problems, not test–retest).

## 6. Log

| date | step | result |
|---|---|---|
| 2026-11-18 | plan drafted | — |
| 2026-11-18 | AFNI on the cluster | Registration failed (`-cost lpa+` unknown): the login PATH used `/usr/pubsw/packages/AFNI/new` = 18.0.22 (2018). Switched to `/usr/pubsw/packages/AFNI/26.0.10` (`~/.cshrc` lines 38 and 84). Cohort preprocessing restarted with 26.0.10. |
| 2026-11-18 | `calibrateQuality.m` | Written; tested on `bch046` (scratch copy, both diagnostic runs). Two findings: clean-WM mask smaller than Phase 9's; zshimEcho p99 driven by one 23-voxel cluster. See the status section. |
| 2026-11-18 | decisions 1–4 (section 5) | healthy cohort = ctrl + mdt; both `dBCorr` settings; second runs not used |
| 2026-11-24 | anatomy route comparison | `Scripts/compareAnatRoutes.m` on bch046: FreeSurfer + 12-parameter registration matches today's geometry (median 0.56 mm), removes two resampling steps and the BOLD dependency; rigid registration misses a 3–4% anterior–posterior compression. See [anat-route-report.md](anat-route-report.md). |
