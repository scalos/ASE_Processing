function q = qualityDefaults(q)
% qualityDefaults  Fill the "quality" option block of ASE_fitting with defaults.
%   q = qualityDefaults()      -> all defaults (every feature off)
%   q = qualityDefaults(q)     -> user struct with missing fields filled in
% Unknown fields raise an error (catches typos in the config).
%
% Every feature is OFF by default, so the default pipeline output is unchanged.
% Thresholds marked "to calibrate" are placeholders until validation.
    arguments
        q = struct()
    end
    if isempty(q), q = struct(); end
    assert(isstruct(q) && isscalar(q), 'quality options must be a struct (JSON object).');

    d = struct();
    % ---- master switch: build the voxel mask and apply it to the ROI stats ----
    d.enabled          = false;   % compute keep/flags and exclude voxels from ROI means
    % ---- diagnostics (computed automatically when needed by an enabled criterion) ----
    d.repNoise         = false;   % Phase 2: sigma / SNR maps from the repetitions
    d.zshimEstimate    = false;   % Phase 3: per-voxel Gz estimate from the z-shim frames
    d.fitExtensions    = false;   % Phase 4: vCBV_quad, vCBV z-map, per-echo R2'
    % ---- z-shim gradient search / correction (correction itself: opts.dBCorr) ----
    d.zshim_GzRange    = [-150 150];  % uT/m
    d.zshim_GzStep     = 1;           % uT/m
    d.zshim_Fmin       = 0.5;         % 'heavily attenuated' long frame (F < Fmin), counted for bit 8;
                                      %   also the NaN cutoff when zshim_weightF2 is false (legacy)
    d.zshim_weightF2   = true;        % dBCorr: weight each corrected voxel-frame by F^2 in the fits
    d.zshim_Fcut       = 0.1;         % dBCorr with F^2 weights: frames with F below this -> NaN
    % ---- ROI integration ----
    d.noiseWeighting   = false;   % weight ROI membership by S^2/sigma^2 (off until compared)
    d.minVoxKept       = 5;       % effective voxels; fewer -> region values NaN
    % ---- criteria (bit = position in quality flags, 0-based) ----
    d.useBounds        = true;    % bit 0: implausible raw values (physBd_* bounds)
    d.useQuadSign      = true;    % bit 1: k >= 0 or vCBV_quad <= 0
    d.useVcbvZ         = true;    % bit 2: |z(vCBV_int - vCBV_quad)| > vcbvZ_max
    d.vcbvZ_max        = 3;       %   to calibrate
    d.useSNR           = true;    % bit 3: SNR < snr_min
    d.snr_min          = 10;      %   to calibrate
    d.useCSF           = true;    % bit 4: CSF fraction > csf_max (needs qualityMaskFiles)
    d.csf_max          = 0.2;
    d.useVein          = true;    % bit 5: vCBV_int > vein_vCBVmax
    d.vein_vCBVmax     = 0.15;
    d.useZshimEdge     = true;    % bit 6: Gz at grid edge / peak outside sampled moments
    d.useZshimEcho     = true;    % bit 7: |Gz_e1 - Gz_e2| > zshimEcho_maxDiff
    d.zshimEcho_maxDiff = 20;     %   uT/m, to calibrate
    d.useZshimFrac     = true;    % bit 8: fraction of long frames with F < Fmin > zshimFrac_max
    d.zshimFrac_max    = 0.3;
    d.useR2PEcho       = false;   % bit 9: per-echo R2' F-test p < r2pEcho_pMax (off: flags a systematic
                                  %   echo-1 vs echo-2 R2' difference, not bad voxels; see report Q3)
    d.r2pEcho_pMax     = 0.01;    %   to calibrate

    fn = fieldnames(q);
    unknown = setdiff(fn, fieldnames(d));
    if ~isempty(unknown)
        error('Unknown quality option(s): %s', strjoin(unknown, ', '));
    end
    for i = 1:numel(fn)
        d.(fn{i}) = q.(fn{i});
    end
    q = d;
end
