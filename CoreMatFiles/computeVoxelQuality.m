function [keep, flags, bitInfo, info] = computeVoxelQuality(rawMaps, qual, physBds, q)
% computeVoxelQuality  Voxel exclusion mask and per-criterion bit flags.
%
%   rawMaps : struct with raw (unclamped) voxel maps R2, R2P, vCBV, OEF (X x Y x Z)
%   qual    : struct of quality diagnostics (any subset of):
%               SNR (X x Y x Z x nEchos), csfFrac, Gz (x nEchos), zshimEdge (x nEchos),
%               zshimFracLow, vCBV_quad, vCBV_z, kZ (k/SE_k), R2P_pval
%   physBds : struct with [lo hi] bounds per map (fields R2, R2P, vCBV, OEF)
%   q       : resolved quality options (qualityDefaults)
%
%   keep    : logical X x Y x Z, true where no enabled criterion is set
%   flags   : uint16 X x Y x Z bitmask; bit b (0-based) = criterion b (see bitTable below)
%   bitInfo : table with Bit, Name, Enabled, Available, nFlagged
%   info    : struct with the vein reference (veinRef: median tissue vCBV, veinThreshold,
%             veinRefN: voxels in the reference); NaN when not computed
%
% Vein (bit 5), q.vein_mode:
%   "relative" : vCBV > vein_relK * median(vCBV) over tissue voxels, i.e. csfFrac <
%                1 - vein_refTissueMin with all four raw maps finite and inside physBds.
%                Scales with the subject's (and the analysis') vCBV level; needs csfFrac.
%   "absolute" : vCBV > vein_vCBVmax.
%
% A disabled criterion never sets its bit. An enabled criterion whose input is missing is
% reported as Available = false (with a warning) and sets nothing. NaN diagnostics never set a
% bit, except bit 0, where a NaN fit value counts as implausible.
    arguments
        rawMaps struct
        qual struct
        physBds struct
        q struct
    end
    sz = size(rawMaps.vCBV, 1:3);
    flags = zeros(sz, 'uint16');
    names = bitTable();
    nBits = numel(names);
    enabled = false(nBits,1);  avail = false(nBits,1);  nFlag = zeros(nBits,1);

    function setBit(b, mask)
        mask = reshape(logical(mask), sz);
        flags(mask) = bitor(flags(mask), uint16(2^b));
        nFlag(b+1) = nnz(mask);
    end
    function tf = has(varargin)
        tf = all(cellfun(@(f) isfield(qual,f) && ~isempty(qual.(f)), varargin));
    end

    % bit 0: implausible raw values (or no fit)
    enabled(1) = q.useBounds;  avail(1) = true;
    if q.useBounds
        bad = false(sz);
        for f = {'R2','R2P','vCBV','OEF'}
            nm = f{1};
            if ~isfield(rawMaps,nm) || isempty(rawMaps.(nm)), continue; end
            x = rawMaps.(nm);
            if strcmp(nm,'OEF')
                % OEF is NaN by construction where vCBV<=0 or R2P<=0 (already caught there)
                bad = bad | x < physBds.(nm)(1) | x > physBds.(nm)(2);
            else
                bad = bad | isnan(x) | x < physBds.(nm)(1) | x > physBds.(nm)(2);
            end
        end
        setBit(0, bad);
    end
    % bit 1: short-regime curvature significantly of the wrong sign: k/SE(k) > quadSign_kZmin
    % (k should be < 0; noise alone gives k >= 0 in many voxels, so plain k >= 0 is not used)
    enabled(2) = q.useQuadSign;  avail(2) = has('kZ');
    if enabled(2) && avail(2), setBit(1, qual.kZ > q.quadSign_kZmin); end
    % bit 2: vCBV_int vs vCBV_quad disagreement
    enabled(3) = q.useVcbvZ;  avail(3) = has('vCBV_z');
    if enabled(3) && avail(3), setBit(2, abs(qual.vCBV_z) > q.vcbvZ_max); end
    % bit 3: low SNR (lowest SNR over echos)
    enabled(4) = q.useSNR;  avail(4) = has('SNR');
    if enabled(4) && avail(4), setBit(3, min(qual.SNR,[],4) < q.snr_min); end
    % bit 4: CSF partial volume
    enabled(5) = q.useCSF;  avail(5) = has('csfFrac');
    if enabled(5) && avail(5), setBit(4, qual.csfFrac > q.csf_max); end
    % bit 5: large vein (raw intercept vCBV), absolute or relative to the subject's tissue vCBV
    info = struct('veinRef', NaN, 'veinThreshold', NaN, 'veinRefN', 0);
    enabled(6) = q.useVein;
    relVein = strcmpi(q.vein_mode, 'relative');
    avail(6) = ~relVein || has('csfFrac');
    if enabled(6) && avail(6)
        if relVein
            ok = qual.csfFrac < 1 - q.vein_refTissueMin;
            for f = {'R2','R2P','vCBV','OEF'}
                x = rawMaps.(f{1});
                ok = ok & isfinite(x) & x >= physBds.(f{1})(1) & x <= physBds.(f{1})(2);
            end
            info.veinRefN = nnz(ok);
            info.veinRef = median(rawMaps.vCBV(ok));
            info.veinThreshold = q.vein_relK * info.veinRef;
            if info.veinRefN < 100 || ~(info.veinRef > 0)
                warning('computeVoxelQuality:veinReference', ...
                    'Vein reference from %d tissue voxels (median vCBV %g): threshold may be unreliable.', ...
                    info.veinRefN, info.veinRef);
            end
        else
            info.veinThreshold = q.vein_vCBVmax;
        end
        if isfinite(info.veinThreshold), setBit(5, rawMaps.vCBV > info.veinThreshold); end
    end
    % bit 6: z-shim estimate at the grid edge / outside the sampled moments (any echo)
    enabled(7) = q.useZshimEdge;  avail(7) = has('zshimEdge');
    if enabled(7) && avail(7), setBit(6, any(qual.zshimEdge > 0, 4)); end
    % bit 7: echo 1 vs echo 2 Gz disagreement
    enabled(8) = q.useZshimEcho;  avail(8) = has('Gz') && size(qual.Gz,4) >= 2;
    if enabled(8) && avail(8), setBit(7, abs(qual.Gz(:,:,:,1) - qual.Gz(:,:,:,2)) > q.zshimEcho_maxDiff); end
    % bit 8: too many long frames attenuated below Fmin
    enabled(9) = q.useZshimFrac;  avail(9) = has('zshimFracLow');
    if enabled(9) && avail(9), setBit(8, qual.zshimFracLow > q.zshimFrac_max); end
    % bit 9: R2' inconsistent between echos
    enabled(10) = q.useR2PEcho;  avail(10) = has('R2P_pval');
    if enabled(10) && avail(10), setBit(9, qual.R2P_pval < q.r2pEcho_pMax); end

    missing = find(enabled & ~avail);
    for i = missing(:)'
        warning('computeVoxelQuality:missingInput', ...
            'Quality criterion bit %d (%s) is enabled but its input is unavailable; skipped.', i-1, names{i});
    end
    keep = flags == 0;
    bitInfo = table((0:nBits-1).', names(:), enabled, avail, nFlag, ...
        'VariableNames', {'Bit','Name','Enabled','Available','nFlagged'});
end

function names = bitTable()
    names = {'implausible','quadSign','vcbvZ','lowSNR','csf','vein', ...
             'zshimEdge','zshimEcho','zshimFrac','r2pEcho'};
end
