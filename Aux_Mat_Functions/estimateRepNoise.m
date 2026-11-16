function [sigma, S0] = estimateRepNoise(echoData, repInds, opts)
% estimateRepNoise  Voxelwise noise level from repeated ASE sweeps of one echo.
%
%   echoData : X x Y x Z x ndTE, unaveraged frames of ONE echo (acquisition order,
%              as stored in the processed file)
%   repInds  : nReps x nSweep frame indices (info.repInds from readASEHeader);
%              column i of every row is the same dTE point
%   opts.baseInds : optional sweep positions (1..nSweep) used for the baseline signal
%              S0 (default: all positions; callers pass the dTE = 0 positions)
%
%   sigma : X x Y x Z, noise SD in signal units (NaN where it cannot be estimated)
%   S0    : X x Y x Z, mean signal over reps at opts.baseInds (for SNR = S0./sigma)
%
% Method (two reps): d_i = S_rep1,i - S_rep2,i over the nSweep points; the voxel's mean d
% (drift / rep offset) is removed and sigma^2 = sum((d - mean(d)).^2) / (2*(nSweep-1)).
% For nReps > 2 the same idea is pooled over points: the residual of the two-way
% (rep x point) additive model, i.e. per-point sample variances after removing each rep's
% offset, with (nSweep-1)*(nReps-1) degrees of freedom. For nReps = 2 the two forms are
% identical.
    arguments
        echoData {mustBeNumeric}
        repInds (:,:) {mustBeInteger, mustBePositive}
        opts.baseInds = []
    end
    [nReps, nSweep] = size(repInds);
    sz = size(echoData); sz(end+1:4) = 1;
    assert(max(repInds(:)) <= sz(4), ...
        'repInds reference frame %d but the data has %d frames.', max(repInds(:)), sz(4));
    if nReps < 2 || nSweep < 2
        warning('estimateRepNoise:tooFew', ...
            'Need >= 2 reps and >= 2 sweep points to estimate noise (got %d x %d); returning NaN.', nReps, nSweep);
        sigma = NaN(sz(1:3));
        S0 = NaN(sz(1:3));
        return
    end

    nVox = prod(sz(1:3));
    Y = reshape(double(echoData), nVox, sz(4));
    % nVox x nSweep x nReps
    S = zeros(nVox, nSweep, nReps);
    for r = 1:nReps
        S(:,:,r) = Y(:, repInds(r,:));
    end
    clear Y

    ptMean  = mean(S, 3);                        % per-point mean over reps
    repOff  = mean(S - ptMean, 2);               % per-rep offset (drift), nVox x 1 x nReps
    resid   = S - ptMean - repOff;
    dof     = (nSweep - 1) * (nReps - 1);
    sigma   = sqrt(sum(resid.^2, [2 3]) / dof);

    baseInds = opts.baseInds;
    if isempty(baseInds), baseInds = 1:nSweep; end
    S0 = mean(ptMean(:, baseInds), 2);

    bad = any(~isfinite(S), [2 3]);
    sigma(bad) = NaN;  S0(bad) = NaN;
    sigma = reshape(sigma, sz(1:3));
    S0    = reshape(S0,    sz(1:3));
end
