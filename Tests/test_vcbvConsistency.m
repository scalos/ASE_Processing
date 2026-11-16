function tests = test_vcbvConsistency
% Monte-Carlo check of the analytic SEs in vcbvConsistency (2-echo protocol of the test subject).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
end

function testSEMatchesMonteCarlo(tc)
    rng(7);
    nMC = 20000;
    dte = repmat((0:44).',2,1);                       % 2 reps
    tau = [dte; dte + 32]/1000;  TE = [61.2*ones(90,1); 93.2*ones(90,1)]/1000;
    echo = [ones(90,1); 2*ones(90,1)];
    short = tau <= 0.014;  long = tau > 0.020;
    zeta = 0.04; R2 = 12; R2P = 6; dw = R2P/zeta; c = log(400);
    lnS = c - R2*TE - R2P*tau + zeta;                 % long-regime asymptote
    lnS(short) = c - R2*TE(short) - 0.3*zeta*(dw*tau(short)).^2;
    S = exp(lnS);
    sigma = [4; 3.5];                                 % per echo, signal units
    Y = S.' + sigma(echo).' .* randn(nMC, numel(S));
    Y = reshape(Y, nMC,1,1,[]);
    yL = Y(:,:,:,long);  yS = Y(:,:,:,short);
    XL = [ones(nnz(long),1), -TE(long), -tau(long)];
    [BL,~,~] = fitMultiRegressor(yL, XL);
    TES = TE(short);
    modS = yS .* exp(BL(:,:,:,2) .* reshape(TES,1,1,1,[]));
    XS = [ones(nnz(short),1), tau(short).^2];
    [BS,~,~] = fitMultiRegressor(modS, XS);
    fr.vCBV = BL(:,:,:,1) - BS(:,:,:,1);  fr.R2P = BL(:,:,:,3);  fr.k = BS(:,:,:,2);
    fr.vCBV_quad = -0.3*fr.R2P.^2./fr.k;  fr.vCBV_quad(~(fr.k<0)) = NaN;
    sl = struct('sigma', repmat(reshape(sigma,1,1,1,2), nMC,1,1,1), ...
                'echoL', echo(long), 'facL', ones(nnz(long),1), 'echoS', echo(short), 'facS', ones(nnz(short),1));
    o = vcbvConsistency(yL, yS, XL, XS, TES, fr, sl);
    sdInt = std(fr.vCBV(:));  sdQ = std(fr.vCBV_quad(:),'omitnan');
    fprintf('vCBV_int: MC SD %.4g, analytic SE %.4g | vCBV_quad: MC SD %.4g (robust %.4g), SE %.4g | z: mean %.3f SD %.3f\n', ...
        sdInt, median(o.SE_int(:)), sdQ, iqr(fr.vCBV_quad(:))/1.349, median(o.SE_quad(:),'omitnan'), mean(o.z(:),'omitnan'), std(o.z(:),'omitnan'));
    tc.verifyEqual(median(o.SE_int(:)), sdInt, 'RelTol', 0.05);
    tc.verifyEqual(median(o.SE_quad(:),'omitnan'), iqr(fr.vCBV_quad(:))/1.349, 'RelTol', 0.15);
    zr = iqr(o.z(:))/1.349;
    fprintf('z: median %.3f robust SD %.3f, frac |z|>3 = %.4f (Gaussian 0.0027)\n', median(o.z(:),'omitnan'), zr, mean(abs(o.z(:))>3));
    % vCBV_quad is a ratio (nonlinear), so z is slightly biased and not exactly N(0,1);
    % tolerances reflect the delta-method approximation, not exactness.
    tc.verifyEqual(zr, 1, 'AbsTol', 0.15);
    tc.verifyLessThan(abs(median(o.z(:),'omitnan')), 0.25);
    tc.verifyLessThan(mean(abs(o.z(:))>3), 0.01);
end

function testWeightedSEMatchesMonteCarlo(tc)
% Frames attenuated by a through-slice factor F (as with a large Gz), divided by F and fitted
% with F^2 weights: the analytic SE must still match the Monte-Carlo spread, and the weighted
% fit must be more precise than the unweighted one on the same corrected data.
    rng(11);
    nMC = 20000;
    dte = repmat((0:44).',2,1);
    tau = [dte; dte + 32]/1000;  TE = [61.2*ones(90,1); 93.2*ones(90,1)]/1000;
    echo = [ones(90,1); 2*ones(90,1)];
    short = tau <= 0.014;  long = tau > 0.020;
    zeta = 0.04; R2 = 12; R2P = 6; dw = R2P/zeta; c = log(400);
    lnS = c - R2*TE - R2P*tau + zeta;
    lnS(short) = c - R2*TE(short) - 0.3*zeta*(dw*tau(short)).^2;
    S = exp(lnS);
    x = 2.67502e8*0.0048*(45e-6)*tau/2;               % Gz = 45 uT/m, 4.8 mm slice
    F = abs(sin(x)./x);  F(x==0) = 1;                 % down to ~0.2 at the longest tau
    sigma = [4; 3.5];
    Ym = (S.*F).' + sigma(echo).' .* randn(nMC, numel(S));
    Y = reshape(Ym ./ F.', nMC,1,1,[]);
    w = repmat(reshape(F.^2,1,1,1,[]), nMC,1,1,1);
    yL = Y(:,:,:,long);  yS = Y(:,:,:,short);  wL = w(:,:,:,long);  wS = w(:,:,:,short);
    XL = [ones(nnz(long),1), -TE(long), -tau(long)];
    XS = [ones(nnz(short),1), tau(short).^2];
    TES = TE(short);
    fitOne = @(wl, ws) localFit(yL, yS, XL, XS, TES, wl, ws);
    frW = fitOne(wL, wS);  frU = fitOne([], []);
    sl = struct('sigma', repmat(reshape(sigma,1,1,1,2), nMC,1,1,1), ...
                'echoL', echo(long), 'facL', ones(nnz(long),1), 'echoS', echo(short), 'facS', ones(nnz(short),1), ...
                'wL', wL, 'wS', wS);
    o = vcbvConsistency(yL, yS, XL, XS, TES, frW, sl);
    sdW = std(frW.vCBV(:));  sdU = std(frU.vCBV(:));
    fprintf('weighted: MC SD %.4g, analytic SE %.4g | unweighted MC SD %.4g\n', sdW, median(o.SE_int(:)), sdU);
    tc.verifyEqual(median(o.SE_int(:)), sdW, 'RelTol', 0.05);
    tc.verifyLessThan(sdW, sdU);
    tc.verifyEqual(mean(frW.vCBV(:)), zeta, 'AbsTol', 3*sdW/sqrt(nMC) + 0.002);
end

function fr = localFit(yL, yS, XL, XS, TES, wL, wS)
    [BL,~,~] = fitMultiRegressor(yL, XL, [], [], wL);
    modS = yS .* exp(BL(:,:,:,2) .* reshape(TES,1,1,1,[]));
    [BS,~,~] = fitMultiRegressor(modS, XS, [], [], wS);
    fr.vCBV = BL(:,:,:,1) - BS(:,:,:,1);  fr.R2P = BL(:,:,:,3);  fr.k = BS(:,:,:,2);
    fr.vCBV_quad = -0.3*fr.R2P.^2./fr.k;  fr.vCBV_quad(~(fr.k<0)) = NaN;
end
