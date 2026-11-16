function T = regimeBoundarySim(opts)
% regimeBoundarySim  Bias of the ASE regime fits against the full static-dephasing model.
%
% Generates noise-free signal from the Yablonskiy-Haacke cylinder model
%     ln S = ln C - R2*TE - zeta*f_s(dw*tau),
%     f_s(x) = (1/3) * int_0^1 (2+u) sqrt(1-u) (1 - J0(1.5 x u)) / u^2 du
% on the protocol's tau/TE sampling, then fits it exactly as ASE_fitting does (regime split
% at 1.5*t_crit -/+ tCrit_buffer, long design [1 -TE -tau], short design [1 tau^2] on
% ln S + R2*TE, vCBV = long intercept - short intercept, OEF = R2'/(vCBV*oefScale)), plus
% the per-echo long-regime R2' fits behind the bit-9 echo test.
%
% "OEF" here is (1 - Y_venous), as in ASE_fitting; dw = oefScale*OEF, true t_c = 1/dw.
% t_crit uses StO2_max as in ASE_fitting: t_c(assumed) = 1/(oefScale*(1 - StO2_max)).
%
% Example (the 2026-11-16 analysis in docs/quality-mask-report.md, Q3):
%   regimeBoundarySim('oef',[0.18 0.20 0.25 0.35 0.45], 'buffers',[3 3; 3 6; 3 9; 3 12; 3 15; 6 3; 1 3])
%   regimeBoundarySim('StO2_max',0.8, 'buffers',[3 3])
%
% Noise (nMC > 0): Monte-Carlo with Gaussian noise of SD 1 on the magnitude signal, scaled so
% that echo-1 SNR at tau = 0 is opts.SNR, and every tau frame acquired nReps times (as the
% stacked, unaveraged fit). The same noise realisations are used for every buffer/StO2_max
% at a given (zeta, OEF), so settings are compared pairwise. Reported per quantity:
%   med*  : median relative error (%)        sd*  : robust SD (IQR/1.349) relative to truth (%)
%   rmse* : sqrt(med^2 + sd^2) (%)           failOEF : fraction with no finite OEF
%   regimeBoundarySim('nMC',5000,'SNR',100,'StO2_max',0.7,'buffers',[3 3])
%
% Output table: one row per (zeta, OEF, buffer) with fitted values and relative errors (%).
% Needs CoreMatFiles on the path (fitMultiRegressor).
    arguments
        opts.oef (1,:) double = [0.20 0.25 0.35 0.45]
        opts.zeta (1,:) double = 0.06
        opts.buffers (:,2) double = [3 3; 3 6; 3 9; 3 12; 3 15]   % [below above] t_crit, ms
        opts.StO2_max (1,1) double = 0.6
        opts.R2 (1,1) double = 12                                  % 1/s
        opts.TE (1,:) double = [61.2 93.2]                         % ms, per echo
        opts.dTE (:,1) double = (0:44)'                            % ms, per-echo tau sweep
        opts.B0 (1,1) double = 3
        opts.Hct (1,1) double = 0.42
        opts.Hct_correction (1,1) double = 0.85
        opts.dChi_0 (1,1) double = 0.18e-6
        opts.verbose (1,1) logical = true
        opts.nMC (1,1) double = 0                                  % 0 = noise-free only
        opts.SNR (1,1) double = 100                                % echo-1, tau = 0, per acquisition
        opts.nReps (1,1) double = 2
        opts.seed (1,1) double = 1
    end
    gamma = 2.67502e8;
    oefScale = (4/3)*pi*gamma*opts.dChi_0*opts.Hct*opts.Hct_correction*opts.B0;
    tcrit = round(1.5*1000/(oefScale*(1-opts.StO2_max)), 3);          % ms, as in ASE_fitting
    fs = @(x) arrayfun(@(xx) integral(@(u) (2+u).*sqrt(1-u).*(1-besselj(0,1.5*xx*u))./u.^2, 0, 1)/3, x);

    % echo e frames: tau = dTE + (TE_e - TE_1), as taus_full in ASE_fitting
    nE = numel(opts.TE);  nD = numel(opts.dTE);
    tau = zeros(nD*nE,1);  TE = tau;  echo = tau;
    for e = 1:nE
        k = (e-1)*nD + (1:nD);
        tau(k) = (opts.dTE + opts.TE(e) - opts.TE(1))/1000;
        TE(k) = opts.TE(e)/1000;  echo(k) = e;
    end
    if opts.verbose
        fprintf('Assumed t_c (StO2_max=%.2f) = %.2f ms; boundary 1.5 t_c = %.2f ms\n', ...
            opts.StO2_max, tcrit/1.5, tcrit);
    end

    rows = {};
    for zeta = opts.zeta
    for oef = opts.oef
        dw = oefScale*oef;  R2P = zeta*dw;
        lnS = log(500) - opts.R2*TE - zeta*fs(dw*tau);
        Y = @(m) reshape(exp(lnS(m)),1,1,1,[]);
        if opts.verbose
            fprintf('\nzeta=%.3f OEF=%.2f: true t_c=%.1f ms (1.5 t_c=%.1f ms), R2''=%.2f 1/s\n', ...
                zeta, oef, 1000/dw, 1500/dw, R2P);
            fprintf('  buf[-,+]  nShort nLong(e1) | R2''e1  R2''e2  gap%% | R2''    R2     vCBV   OEF  (rel. error %%)\n');
        end
        for b = 1:size(opts.buffers,1)
            lo = max(0, tcrit - opts.buffers(b,1));  hi = tcrit + opts.buffers(b,2);
            sh = tau*1000 <= lo;  lg = tau*1000 > hi;
            BL = squeeze(fitMultiRegressor(Y(lg), [ones(nnz(lg),1), -TE(lg), -tau(lg)]));
            modS = exp(lnS(sh) + BL(2)*TE(sh));
            BS = squeeze(fitMultiRegressor(reshape(modS,1,1,1,[]), [ones(nnz(sh),1), tau(sh).^2]));
            vC = BL(1) - BS(1);  oefF = BL(3)/(vC*oefScale);
            r = NaN(1,2);
            for e = 1:min(nE,2)
                m = lg & echo == e;
                if nnz(m) >= 2
                    Be = squeeze(fitMultiRegressor(Y(m), [ones(nnz(m),1), -tau(m)]));  r(e) = Be(2);
                end
            end
            gap = 100*(r(2)-r(1))/r(2);
            err = 100*([BL(3)/R2P, BL(2)/opts.R2, vC/zeta, oefF/oef] - 1);
            if opts.verbose
                fprintf('  [%2g %2g]   %4d   %4d      | %5.2f  %5.2f  %5.1f | %+6.1f %+6.1f %+6.1f %+6.1f\n', ...
                    opts.buffers(b,:), nnz(sh), nnz(lg & echo==1), r(1), r(2), gap, err);
            end
            mc = NaN(1,10);
            if opts.nMC > 0
                mc = noiseMC(lnS, tau, TE, sh, lg, opts, zeta, oef, oefScale, R2P);
                if opts.verbose
                    fprintf('            noise (SNR %g, %d reps): OEF med %+5.1f sd %5.1f rmse %5.1f fail %4.1f%% | vCBV med %+5.1f sd %5.1f | R2 sd %4.1f | R2'' sd %4.1f\n', ...
                        opts.SNR, opts.nReps, mc(1), mc(2), mc(3), 100*mc(4), mc(5), mc(6), mc(8), mc(10));
                end
            end
            rows(end+1,:) = {zeta, oef, 1000/dw, opts.buffers(b,1), opts.buffers(b,2), nnz(sh), ...
                nnz(lg & echo==1), r(1), r(2), gap, err(1), err(2), err(3), err(4), ...
                mc(1), mc(2), mc(3), mc(4), mc(5), mc(6), mc(7), mc(8), mc(9), mc(10)}; %#ok<AGROW>
        end
    end
    end
    T = cell2table(rows, 'VariableNames', {'zeta','OEF','tc_true_ms','bufBelow','bufAbove', ...
        'nShort','nLongEcho1','R2P_e1','R2P_e2','echoGapPct','errR2P','errR2','errVCBV','errOEF', ...
        'medOEF','sdOEF','rmseOEF','failOEF','medVCBV','sdVCBV','medR2','sdR2','medR2P','sdR2P'});
end

function mc = noiseMC(lnS, tau, TE, sh, lg, opts, zeta, oef, oefScale, R2P)
% Monte-Carlo fit of noisy data; same realisations for every call with the same seed/OEF/zeta.
    rng(opts.seed + round(1e4*oef) + round(1e6*zeta));
    S = exp(lnS - lnS(1)) * opts.SNR;                       % echo-1 tau=0 -> SNR (sigma = 1)
    nF = numel(S);  nR = opts.nReps;
    Y = repmat(S.', opts.nMC, nR) + randn(opts.nMC, nF*nR);  % frames: [rep1 frames, rep2 frames, ...]
    Y = reshape(Y, opts.nMC, 1, 1, []);
    rep = @(v) repmat(v, nR, 1);
    tauR = rep(tau);  TER = rep(TE);  shR = rep(sh);  lgR = rep(lg);
    BL = fitMultiRegressor(Y(:,:,:,lgR), [ones(nnz(lgR),1), -TER(lgR), -tauR(lgR)]);
    modS = Y(:,:,:,shR) .* exp(BL(:,:,:,2) .* reshape(TER(shR),1,1,1,[]));
    BS = fitMultiRegressor(modS, [ones(nnz(shR),1), tauR(shR).^2]);
    vC = BL(:,:,:,1) - BS(:,:,:,1);  r2p = BL(:,:,:,3);  r2 = BL(:,:,:,2);
    O = r2p ./ (vC*oefScale);  O(vC <= 0 | r2p <= 0) = NaN;     % as in fitASERegimes
    st = @(x, t) [100*(median(x(:),'omitnan')/t - 1), 100*iqr(x(:))/1.349/t];
    so = st(O, oef);  sv = st(vC, zeta);  s2 = st(r2, opts.R2);  sp = st(r2p, R2P);
    mc = [so(1), so(2), hypot(so(1), so(2)), mean(~isfinite(O(:))), sv, s2, sp];
end
