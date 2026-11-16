function out = vcbvConsistency(aseLong, aseShort, XL, XS, TEsShort, fitRes, sigLog, opts)
% vcbvConsistency  Standard errors of vCBV_int and vCBV_quad and the z-score of their difference.
%
%   aseLong, aseShort : X x Y x Z x NL / NS signal used by fitASERegimes (NaN/<=0 = excluded)
%   XL, XS            : long / short OLS design matrices exactly as used in fitASERegimes
%                       (long: [1 -TE -tau] for 2+ echos or [1 -tau] for one echo;
%                        short: [1 tau^2])
%   TEsShort          : NS x 1 TE of each short frame (s); used for the R2 coupling
%   fitRes            : output of fitASERegimes (needs BL, BS, R2P, vCBV, k, vCBV_quad)
%   sigLog.sigma      : X x Y x Z x nEchos noise SD (signal units, per single acquisition)
%   sigLog.echoL/S    : NL/NS x 1 echo index of each frame
%   sigLog.facL/S     : NL/NS x 1 noise factor per frame (1/sqrt(nAvg) for averaged frames)
%   sigLog.wL/wS      : optional X x Y x Z x NL/NS fit weights (F^2 from the z-shim correction),
%                       exactly as passed to fitMultiRegressor. With weights the projection is
%                       P = (X'WX)^-1 X'W, and since the data were divided by F the measured
%                       signal is F*S, so the log-space SD is sigma*fac/(sqrt(w)*S).
%   opts.chunk        : voxels per chunk (default 2e4)
%
% Noise model: independent frames, log-space SD = sigma_echo*fac / S (S = measured signal).
% Both estimates are linear (vCBV_int) or linearized (vCBV_quad, delta method) functions of the
% log data, so with per-voxel OLS projection P = (X'WX)^-1 X'W (W = valid frames):
%   vCBV_int  = P_L(1,:) yL - P_S(1,:) (yS + R2*TE_S),   R2 = P_L(2,:) yL   (2+ echos)
%   vCBV_quad = -0.3 R2'^2 / k,  R2' = P_L(end,:) yL,  k = P_S(2,:) (yS + R2*TE_S)
% Var(.) = sum_f g_f^2 sigLog_f^2 for the gradient g of each quantity w.r.t. the log data.
%
% out: SE_int, SE_quad, z (= (vCBV_int - vCBV_quad)/SE(diff)) where k < 0, and SE_k (SE of the
%      short-regime curvature k, incl. the R2 coupling) and kZ = k/SE_k wherever k is finite;
%      all X x Y x Z
    arguments
        aseLong
        aseShort
        XL (:,:) double
        XS (:,:) double
        TEsShort (:,1) double
        fitRes struct
        sigLog struct
        opts.chunk (1,1) double = 2e4
    end
    sz = size(aseLong, 1:3);  nVox = prod(sz);
    NL = size(XL,1);  NS = size(XS,1);  pL = size(XL,2);
    assert(size(XS,2) == 2, 'vcbvConsistency needs a short design with a tau^2 column.');
    twoEcho = pL == 3;
    iR2P = pL;                                     % R2' is the last long coefficient

    out.SE_int  = NaN(sz);  out.SE_quad = NaN(sz);  out.z = NaN(sz);
    out.SE_k = NaN(sz);  out.kZ = NaN(sz);
    R2P  = fitRes.R2P(:);  k = fitRes.k(:);
    dif  = fitRes.vCBV(:) - fitRes.vCBV_quad(:);
    vox  = find(isfinite(fitRes.vCBV(:)) & isfinite(R2P) & isfinite(k));   % SE_k for all of these
    useQ = isfinite(dif) & k < 0;                                           % vCBV_quad defined
    if isempty(vox), return; end

    % column scaling for conditioning (as in fitMultiRegressor)
    cL = sqrt(mean(XL.^2,1)); cL(cL==0) = 1;  XLs = XL./cL;
    cS = sqrt(mean(XS.^2,1)); cS(cS==0) = 1;  XSs = XS./cS;

    YL = reshape(aseLong,  nVox, NL);
    YS = reshape(aseShort, nVox, NS);
    wL = getW(sigLog, 'wL', nVox, NL);  wS = getW(sigLog, 'wS', nVox, NS);
    sig = reshape(sigLog.sigma, nVox, []);
    SE_int = NaN(nVox,1);  SE_q = NaN(nVox,1);  SE_d = NaN(nVox,1);  SE_kk = NaN(nVox,1);

    ws = warning('off','MATLAB:nearlySingularMatrix');
    cleanup = onCleanup(@() warning(ws));
    warning('off','MATLAB:singularMatrix');
    for s = 1:opts.chunk:numel(vox)
        idx = vox(s:min(s+opts.chunk-1, numel(vox)));  nv = numel(idx);
        [PL, sL] = proj(YL(idx,:), XLs, cL, sig(idx,:), sigLog.echoL, sigLog.facL, rowsOf(wL,idx));   % pL x NL x nv
        [PS, sS] = proj(YS(idx,:), XSs, cS, sig(idx,:), sigLog.echoS, sigLog.facS, rowsOf(wS,idx));   % 2 x NS x nv
        PL1 = squeeze3(PL(1,:,:));  PLr = squeeze3(PL(iR2P,:,:));                    % nv x NL
        PS1 = squeeze3(PS(1,:,:));  PS2 = squeeze3(PS(2,:,:));                       % nv x NS
        if twoEcho
            PL2 = squeeze3(PL(2,:,:));
            c1 = sum(PS1 .* TEsShort.', 2);                     % P_S(1,:)*TE_S
            c2 = sum(PS2 .* TEsShort.', 2);                     % P_S(2,:)*TE_S
            aL = PL1 - c1.*PL2;
        else
            PL2 = zeros(nv,NL);  c2 = zeros(nv,1);
            aL = PL1;
        end
        aS = -PS1;
        r = R2P(idx);  kk = k(idx);
        dqR = -0.6*r./kk;  dqk = 0.3*r.^2./kk.^2;
        qL = dqR.*PLr + dqk.*c2.*PL2;
        qS = dqk.*PS2;
        vInt = sum(aL.^2.*sL,2) + sum(aS.^2.*sS,2);
        vQ   = sum(qL.^2.*sL,2) + sum(qS.^2.*sS,2);
        vD   = sum((aL-qL).^2.*sL,2) + sum((aS-qS).^2.*sS,2);
        SE_int(idx) = sqrt(vInt);  SE_q(idx) = sqrt(vQ);  SE_d(idx) = sqrt(vD);
        % k = P_S(2,:)(yS + R2*TE_S): d k/d yS = P_S(2,:), d k/d yL = (P_S(2,:)*TE_S) * P_L(2,:)
        SE_kk(idx) = sqrt(sum((c2.*PL2).^2.*sL,2) + sum(PS2.^2.*sS,2));
    end
    SE_int(~useQ) = NaN;  SE_q(~useQ) = NaN;  SE_d(~useQ) = NaN;           % as before: k < 0 only
    out.SE_int  = reshape(SE_int, sz);
    out.SE_quad = reshape(SE_q,   sz);
    z = dif ./ SE_d;  z(~isfinite(z)) = NaN;
    out.z = reshape(z, sz);
    kZ = k ./ SE_kk;  kZ(~isfinite(kZ)) = NaN;
    out.SE_k = reshape(SE_kk, sz);
    out.kZ   = reshape(kZ, sz);
end

function [P, s2] = proj(Y, Xs, colScale, sig, echoIdx, fac, w)
% Per-voxel (weighted) LS projection rows (unscaled coefficients) and log-space noise variances.
    [nv, N] = size(Y);
    W  = isfinite(Y) & Y > 0;
    if isempty(w)
        Xw = Xs .* reshape(double(W).', N, 1, nv);         % N x p x nv
        A  = pagemtimes(Xw, 'transpose', Xw, 'none');       % p x p x nv  (W is 0/1)
        P  = pagemtimes(pageinv(A), 'none', Xw, 'transpose');   % p x N x nv (scaled coeffs)
    else
        w(~W) = 0;
        Xm = Xs .* reshape(double(W).', N, 1, nv);
        Xw = Xs .* reshape(w.', N, 1, nv);                  % W*X
        A  = pagemtimes(Xm, 'transpose', Xw, 'none');       % X'WX
        P  = pagemtimes(pageinv(A), 'none', Xw, 'transpose');   % (X'WX)^-1 X'W
    end
    P  = P ./ colScale(:);                                  % undo column scaling
    S  = Y;  S(~W) = NaN;
    sd = sig(:, echoIdx(:).') .* fac(:).' ./ S;              % nv x N log-space SD
    if ~isempty(w), sd = sd ./ sqrt(w); end                 % measured signal = F*S, w = F^2
    s2 = sd.^2;  s2(~W | ~isfinite(s2)) = 0;
end

function w = getW(sigLog, name, nVox, N)
    if isfield(sigLog, name) && ~isempty(sigLog.(name))
        w = reshape(double(sigLog.(name)), nVox, N);
        w(~isfinite(w)) = 1;
    else
        w = [];
    end
end

function M = squeeze3(A)
% 1 x N x nv -> nv x N
    M = reshape(A, size(A,2), size(A,3)).';
end

function r = rowsOf(w, idx)
    if isempty(w), r = []; else, r = w(idx,:); end
end
