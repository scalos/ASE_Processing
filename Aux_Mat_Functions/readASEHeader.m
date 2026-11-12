function [dTEs, zShimming, TEs, info] = readASEHeader(hdr_File)
% Reads the combined ASE header (ase_hdr.json, written by the prep script).
%
%   dTEs      : ndTE x 1            dTE per frame (ms), acquisition order (not sorted)
%   zShimming : ndTE x nEchos       z-moment per frame (Tms/m), same frame order as dTEs
%   TEs       : nEchos x 1          echo times (ms), ordered by echo number
%   info      : struct (optional 4th output) with B0, ImagingFrequency, MaxDTE_ms,
%               StepDTE_ms, nSweep, nExtraPerRep, repRows, sortIdx, series, files, ...
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
    if iscell(e), e = [e{:}]; end
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
        info.B0               = h.MagneticFieldStrength;
        info.ImagingFrequency = h.ImagingFrequency;
        info.nEchos           = nEchos;
        info.nReps            = nReps;
        info.MaxDTE_ms        = mx;
        info.StepDTE_ms       = st;
        info.nSweep           = nSweep;
        info.nZshimming       = nZshimming;
        info.repInds          = repInds;
        info.series = cell(nEchos,1);
        info.file  = cell(nEchos,1);
        for i = 1:nEchos
            if isfield(e,'series'), info.series{i} = char(string(e(i).series)); end
            if isfield(e,'file'),  info.file{i}  = cellstr(e(i).file); end
        end
    end
end