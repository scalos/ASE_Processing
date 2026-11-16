function [dTEs, zShimming, TEs, info] = readASEHeader(hdr_File)
% Reads the combined ASE header (ase_hdr.json, written by the prep script).
%
%   dTEs      : ndTE x 1            dTE per frame (ms), acquisition order (not sorted)
%   zShimming : ndTE x nEchos       z-moment per frame (Tms/m), same frame order as dTEs
%   TEs       : nEchos x 1          echo times (ms), ordered by echo number
%   info      : struct (optional 4th output) with B0 (actual field, T; from B0_T or the imaging
%               frequency), B0_nominal, ImagingFrequency, MaxDTE_ms,
%               StepDTE_ms, nSweep, nReps, nZshimming, repInds, zShimInds,
%               zMoments (nZshim x nEchos, T*s/m), SliceThickness_mm (nEchos x 1,
%               NaN if absent from the header), series, file
%
% The header lists frames for ALL reps in one list per echo. Each NIfTI holds one rep.
% Split convention:
%   nSweep = MaxDTE/StepDTE + 1; rows 1..nSweep*nReps are the sweeps (one block per rep);
%   the remaining rows are extra (z-shim) frames, divided equally across reps in order.
% Every rep must then have identical dTE values (errors otherwise); z-moment may differ.
    arguments
        hdr_File {mustBeTextScalar}
    end
    hdr_File = char(hdr_File);
    assert(isfile(hdr_File), 'Header file not found: %s', hdr_File);

    h = jsondecode(fileread(hdr_File));
    required = {'MagneticFieldStrength','ImagingFrequency','echos'};
    for k = 1:numel(required)
        assert(isfield(h, required{k}), 'Header %s is missing field "%s".', hdr_File, required{k});
    end

    e = h.echos;
    if iscell(e)
        % echo entries with different optional fields decode as a cell: pad missing fields
        allF = {};
        for k = 1:numel(e), allF = union(allF, fieldnames(e{k}), 'stable'); end
        for k = 1:numel(e)
            for f = setdiff(allF, fieldnames(e{k}))'
                e{k}.(f{1}) = [];
            end
            e{k} = orderfields(e{k}, allF);
        end
        e = [e{:}];
    end
    e = e(:);
    nEchos = numel(e);
    if isfield(h,'nEchos')
        assert(h.nEchos == nEchos, ...
            'Header lists nEchos = %d but contains %d echo entries (%s).', h.nEchos, nEchos, hdr_File);
    end

    req = {'echo','EchoTime_ms','nFrames','MaxDTE_ms','StepDTE_ms','dTE_ms','ZMoment'};
    for k = 1:numel(req)
        assert(all(isfield(e, req{k})), ...
            'Header %s: echo entries are missing field "%s".', hdr_File, req{k});
    end

    [~, eo] = sort([e.echo]);
    e = e(eo);
    TEs = [e.EchoTime_ms].';

    % ---- per-echo consistency (dTE must match across echos) ----
    dte1 = e(1).dTE_ms(:);
    N = numel(dte1);
    tol = 1e-9;
    for i = 1:nEchos
        d = e(i).dTE_ms(:);  z = e(i).ZMoment(:);
        assert(numel(d) == e(i).nFrames, ...
            'Echo %d: nFrames (%d) does not match number of dTE values (%d).', e(i).echo, e(i).nFrames, numel(d));
        assert(numel(z) == numel(d), ...
            'Echo %d: ZMoment has %d entries but dTE has %d.', e(i).echo, numel(z), numel(d));
        assert(numel(d) == N && all(abs(d - dte1) <= tol), ...
            'dTE values for echo %d differ from echo %d.', e(i).echo, e(1).echo);
    end
    %echos match so:
    dTEs = dte1;
    zShimming = zeros(N,nEchos);
    for i = 1:nEchos
        zShimming(:,i) = e(i).ZMoment(:);
    end
    mx = e(1).MaxDTE_ms;
    st = e(1).StepDTE_ms;
    nSweep = round(mx/st)+1;
    %assume reps divide evenly with zshimming filling remainder
    nReps = round(numel(dTEs)/nSweep);
    nZshimming = mod(numel(dTEs),nSweep);
    repInds = zeros(nReps,nSweep);
    for rep = 1:nReps
        repInds(rep,:) = (nSweep)*(rep-1)+1:(nSweep)*(rep);
    end
    
    if nargout > 3
        info = struct();
        % actual field strength: B0_T written by asePreProcess.sh, else derived from the imaging
        % frequency (MagneticFieldStrength is the nominal value, e.g. 3 for a 2.89 T Prisma)
        if isfield(h,'B0_T') && ~isempty(h.B0_T)
            info.B0 = h.B0_T;
        else
            info.B0 = h.ImagingFrequency / 42.577478518;
        end
        info.B0_nominal       = h.MagneticFieldStrength;
        info.ImagingFrequency = h.ImagingFrequency;
        info.nEchos           = nEchos;
        info.nReps            = nReps;
        info.MaxDTE_ms        = mx;
        info.StepDTE_ms       = st;
        info.nSweep           = nSweep;
        info.nZshimming       = nZshimming;
        info.repInds          = repInds;
        % trailing z-shim frames (raw-file frame indices, same for every echo)
        info.zShimInds        = nSweep*nReps + (1:nZshimming);
        % z-shim moments, nZshim x nEchos, converted from Tms/m (header units) to T*s/m
        info.zMoments         = zShimming(info.zShimInds,:) * 1e-3;
        % slice thickness per echo (mm); NaN where the header does not provide it
        % (only required when z-shim correction is enabled)
        info.SliceThickness_mm = NaN(nEchos,1);
        if isfield(e,'SliceThickness_mm')
            for i = 1:nEchos
                st_i = e(i).SliceThickness_mm;
                if isnumeric(st_i) && isscalar(st_i), info.SliceThickness_mm(i) = st_i; end
            end
        end
        info.series = cell(nEchos,1);
        info.file  = cell(nEchos,1);
        for i = 1:nEchos
            if isfield(e,'series'), info.series{i} = char(string(e(i).series)); end
            if isfield(e,'file'),  info.file{i}  = cellstr(e(i).file); end
        end
    end
end