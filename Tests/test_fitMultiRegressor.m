function tests = test_fitMultiRegressor
% Weighted path of fitMultiRegressor (z-shim F^2 weights).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
end

function testUnitWeightsBitIdentical(tc)
    rng(1);
    [Y, X] = makeData(500, 40);
    Y(3,1,1,5) = NaN;  Y(7,1,1,9) = -1;
    [B0,r0,n0] = fitMultiRegressor(Y, X);
    [B1,r1,n1] = fitMultiRegressor(Y, X, [], [], ones(size(Y)));
    tc.verifyEqual(B1, B0);  tc.verifyEqual(r1, r0);  tc.verifyEqual(n1, n0);
end

function testMatchesLscov(tc)
    rng(2);
    [Y, X] = makeData(50, 40);
    w = rand(size(Y)).^2 + 0.01;
    Y(4,1,1,[2 3]) = NaN;
    [B,r,n] = fitMultiRegressor(Y, X, [], [], w);
    for v = [1 4 20]
        y = squeeze(Y(v,1,1,:));  ww = squeeze(w(v,1,1,:));  ok = isfinite(y) & y > 0;
        b = lscov(X(ok,:), log(y(ok)), ww(ok));
        tc.verifyEqual(squeeze(B(v,1,1,:)), b, 'RelTol', 1e-8);
        res = log(y(ok)) - X(ok,:)*b;
        tc.verifyEqual(r(v)^2*n(v), sum(ww(ok).*res.^2), 'RelTol', 1e-8);   % weighted RSS
    end
end

function [Y, X] = makeData(nv, nf)
    t = linspace(0.02, 0.08, nf).';
    X = [ones(nf,1), -t];
    S = exp(log(500) - 6*t);
    Y = reshape(S.' .* (1 + 0.02*randn(nv, nf)), nv,1,1,nf);
end
