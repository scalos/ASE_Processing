function [err, nFrames] = computeASEError(ydata, taus, TEs, R2, R2P, vCBV, C, opts)
% RMS residual of the full ASE model (Bessel integral) against the raw data,
% in linear signal units, over all valid frames. No dB/sinc correction.
%
%   ydata : X x Y x Z x Nfr x nReps   raw signal (NaN / <=0 = excluded)
%   taus  : Nfr x 1  tau in ms (NaN = frame excluded, e.g. z-shim frames)
%   TEs   : Nfr x 1  TE of each frame in ms
%   R2    : X x Y x Z  (1/s), or [] for single echo (absorbed into C)
%   R2P   : X x Y x Z  (1/s)   raw, unclamped
%   vCBV  : X x Y x Z          raw, unclamped
%   C     : X x Y x Z  amplitude, exp(B_short(:,:,:,1))
    arguments
        ydata
        taus (:,1)
        TEs (:,1)
        R2
        R2P
        vCBV
        C
        opts.mask = []
        opts.chunkSize (1,1) = 5e4
        opts.nu (1,1) = 1000          % u-grid points (matches original)
        opts.nGrid (1,1) = 4000       % fc lookup-table points
        opts.xCap (1,1) = 200         % cap on table range (delta_omega*tau)
        opts.normalize (1,1) logical = false   % divide by mean valid signal
    end

    sz = size(ydata);  sz(end+1:5) = 1;
    Nv = prod(sz(1:3));  Nfr = sz(4);  Nr = sz(5);
    assert(numel(taus)==Nfr && numel(TEs)==Nfr, 'taus/TEs must have one entry per frame.');
    if isempty(R2), R2 = zeros(sz(1:3)); end

    % per-column vectors, ordered like reshape(ydata,Nv,Nfr*Nr): frame fastest, then rep
    tau = repmat(abs(taus(:))/1000, Nr, 1).';      % 1 x N, seconds
    TE  = repmat(TEs(:)/1000,       Nr, 1).';
    frameOK = isfinite(tau) & isfinite(TE);
    tau = tau(frameOK);  TE = TE(frameOK);

    Y2 = reshape(ydata, Nv, Nfr*Nr);
    Y2cols = find(frameOK);

    % delta_omega = R2P/vCBV  (equals Y_coef*OEF_raw)
    dw = R2P(:) ./ vCBV(:);
    valid = R2P(:) > 0 & vCBV(:) > 0 & isfinite(dw) & isfinite(R2(:)) & isfinite(C(:)) & C(:) > 0;
    if ~isempty(opts.mask), valid = valid & logical(opts.mask(:)); end

    err = NaN(Nv,1);  nF = zeros(Nv,1);
    if ~any(valid) || ~any(frameOK)
        err = reshape(err, sz(1:3));  nFrames = reshape(nF, sz(1:3));  return
    end

    % ---- one-time lookup table for fc(x) ----
    du = 1/opts.nu;  u = (0:opts.nu-1)'*du;
    w  = (2+u).*sqrt(1-u)./(u.^2 + eps);
    xmax = min(max(dw(valid))*max(tau), opts.xCap);
    xg   = linspace(0, xmax, opts.nGrid);
    fcg  = (w.' * (1 - besselj(0, 1.5*u*xg))) * du/3;       % 1 x nGrid

    % ---- chunked evaluation ----
    idx = find(valid);
    for s = 1:opts.chunkSize:numel(idx)
        v  = idx(s : min(s+opts.chunkSize-1, numel(idx)));
        Y  = Y2(v, Y2cols); %#ok
        ok = isfinite(Y) & Y > 0;
        Y(~ok) = 0;

        x  = dw(v) .* tau;                                  % nv x N
        fc = interp1(xg, fcg, x, 'linear', 'extrap');
        Sfit = C(v) .* exp(-R2(v).*TE - vCBV(v).*fc);

        res = (Sfit - Y) .* ok;
        n   = sum(ok, 2);
        e   = sqrt(sum(res.^2, 2) ./ n);
        if opts.normalize, e = e ./ (sum(Y,2) ./ n); end
        e(n == 0) = NaN;

        err(v) = e;  nF(v) = n;
    end
    err     = reshape(err, sz(1:3));
    nFrames = reshape(nF,  sz(1:3));
end