# Outputs reference

## Preprocessing (`<subject>/session_ase/`)

| file | from stage | content |
|---|---|---|
| `ase_<series>.nii`, `dcm_hdr_<series>.json` | 1 | converted echo series (98 frames) and dcm2niix sidecar |
| `dcm_hdr_<series>_dtes.txt` | 2 | per frame: DICOM file, frame index, acquisition time (ms since midnight), dTE / z-moment / rep / echo |
| `ase_hdr.json` | 3 | combined header: TEs, dTEs, z-moments, `SliceThickness_mm` per echo |
| `ase_base.nii.gz` | 4 | echo-1 frame 0 (spin echo): motion and registration reference |
| `ase_e<k>_proc.nii.gz`, `ase_e<k>_mc.1D` | 4 | motion-corrected echoes and motion parameters (with `motionMode=echo1`, the echo-2 `.1D` is a copy of echo 1's) |
| `ase_e1_mc_mat.aff12.1D` | 4 | per-frame echo-1 transforms (`motionMode=echo1`) |
| `anat_fs_*_al.nii.gz` | 5 | transferred anatomy and segmentations |
| `ase_base_al*.nii.gz`, `ase_base_al_mat.aff12.1D`, `ase_base_al_inv_mat.aff12.1D` | 6 | ASE → anatomy registration |
| `anat_fs_aparc_aseg_al_ase.nii.gz`, `anat_fs_wmparc_al_ase.nii.gz` | 7 | segmentations in ASE space (on the 1 mm anatomical grid; ROI membership is sub-sampled from these) |
| `anat_fs_tissueFrac_al.nii.gz` | 7b | optional SynthSeg tissue fraction |
| `anat_fs_xvent_mask_al_ase.nii.gz` (or `…tissueFrac…_ase`) | 8 | tissue fraction on the ASE grid, used for the CSF criterion |

## Fit run (`<subject>/session_ase/fits/<timestamp>/`, prefix `ase_fit_results`)

**Always written:**
- **`_manifest.json`:** the frame map (τ, TE, echo, source frame), which frames were used, t_crit, all fit options including the resolved `quality` block, and the path of every file below.
- **`config_resolved.json`:** the full configuration that produced the run.

**Maps** (`writeOpts` contains `maps`):

| file | content |
|---|---|
| `_OEF`, `_R2P`, `_vCBV`, `_R2` | parameter maps, clipped to `physBd_*` (R2, R2′ in 1/s) |
| `_<map>_raw` | unclipped versions (`writeRawMaps`) |
| `_imErr` | RMS misfit of the full Bessel ASE model over all used frames (no z-shim correction), normalized by mean signal when `normErr` |

**Coefficients** (`coeffs`):

| file | content |
|---|---|
| `_CoeffsLong` | 4D: [intercept, R2, R2′] (2 echoes) |
| `_CoeffsShort` | 4D: [intercept, k] |
| `_RMS_L`, `_RMS_S` | log-space fit RMS per regime |

**Quality and correction diagnostics.** Written with `maps` when the matching feature is on:

| file | feature | content |
|---|---|---|
| `_sigma_e<k>`, `_SNR_e<k>` | repNoise | repetition-noise σ (signal units) and spin-echo SNR |
| `_Gz_e<k>`, `_zshimAmp_e<k>`, `_zshimAtten_e<k>`, `_zshimEdge_e<k>` | z-shim | per-echo Gz (µT/m), z-shim profile amplitude, measured attenuation y(M=0)/A, edge bits (1 = grid edge, 2 = peak outside the sampled moments) |
| `_GzJoint` | `zshim_jointGz` | joint Gz used for the correction |
| `_zshimFracLow` | z-shim | fraction of used long frames with F < `zshim_Fmin` |
| `_vCBV_quad`, `_vCBV_SE`, `_vCBVquad_SE`, `_vCBV_z` | fitExtensions | curvature-based vCBV, standard errors, z of the difference (k < 0 only) |
| `_kZ` | fitExtensions | k / SE(k) (bit 1) |
| `_R2P_e1`, `_R2P_e2`, `_R2P_Fstat` | fitExtensions | per-echo long-regime R2′ and the joint-vs-split F statistic |
| `_csfFrac` | quality masks | CSF fraction (1 − tissue fraction) |
| `_qualityMask`, `_qualityFlags` | `quality.enabled` | 1 = kept; 16-bit flag map (bit b = criterion b) |

## ROI table (`_ROI_Stats.csv`; `stats`, requires segmentations)

One row per region of the LUT. **The values come from fitting each region's averaged signal, not from averaging the voxel maps.** These differ: in GM the ROI-fit OEF is about 0.27 while voxel medians are about 0.22. Use the `_med` columns for voxel medians.

| columns | content |
|---|---|
| `OEF`, `R2P`, `vCBV`, `R2` | ROI-fit values, clipped to `physBd_*` |
| `B_L_1..3`, `RMS_L`, `nFrames_L`, `B_S_1..2`, `RMS_S`, `nFrames_S` | long / short coefficients, fit RMS and frames used |
| `nVox` | effective voxel count (occupancy-weighted) |
| `<map>_raw` | unclipped ROI values (`writeRawMaps`) |
| `nVox_kept`, `fracExcluded`, `fracExcl_b<k>_<name>` | with `quality.enabled`: kept voxels, total excluded fraction, and the fraction each criterion flags (overlapping) |
| `OEF_med`, `R2P_med`, `vCBV_med`, `R2_med` | with `quality.enabled`: medians of the voxel maps over kept voxels |
