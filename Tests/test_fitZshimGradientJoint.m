function tests = test_fitZshimGradientJoint
% Joint (all-echo) through-slice gradient fit.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
    M1 = [3819 3055 2291 1527 763 0 -763 -1527].'*1e-9;      % bch046 echo-1 moments (T*s/m)
    tc.TestData.M  = {M1, M1*76/44};                          % echo-2 moments scale with t
    tc.TestData.t  = [0.044; 0.076];
    tc.TestData.dz = [4.8e-3; 4.8e-3];
end

function y = synth(tc, e, Gz, A, snr, n)
    gamma = 2.67502e8;  M = tc.TestData.M{e};
    x = gamma*tc.TestData.dz(e)*(Gz*1e-6*tc.TestData.t(e) - M.')/2;
    f = abs(sin(x)./x); f(x==0) = 1;
    y = reshape(A*f + (A/snr)*randn(n, numel(M)), n,1,1,[]);
end

function testNoiselessExact(tc)
    for G = [-40 0 17.3 60]
        y = {synth(tc,1,G,300,Inf,1), synth(tc,2,G,200,Inf,1)};
        o = fitZshimGradientJoint(y, tc.TestData.M, tc.TestData.t, tc.TestData.dz, 'GzStep', 0.5);
        tc.verifyEqual(o.Gz, G, 'AbsTol', 0.1);
        tc.verifyEqual(squeeze(o.A).', [300 200], 'RelTol', 1e-3);
    end
end

function testSingleEchoMatchesPerEchoFit(tc)
    rng(2);
    y = synth(tc,1,23,400,60,500);
    oj = fitZshimGradientJoint({y}, tc.TestData.M(1), tc.TestData.t(1), tc.TestData.dz(1));
    o1 = fitZshimGradient(y, tc.TestData.M{1}, tc.TestData.t(1), tc.TestData.dz(1));
    tc.verifyEqual(oj.Gz, o1.Gz, 'AbsTol', 1e-10);
    tc.verifyEqual(oj.A, o1.A, 'AbsTol', 1e-8);
end

function testJointMoreRobustThanEcho2Alone(tc)
    % realistic SNR (echo 1 ~80, echo 2 ~60): joint error no worse than echo 1 alone and
    % clearly better than echo 2 alone at gradients where echo 2's profile is poorly sampled
    rng(5);
    for G = [-30 10 40]
        y1 = synth(tc,1,G,500,80,3000);  y2 = synth(tc,2,G,350,60,3000);
        oj = fitZshimGradientJoint({y1,y2}, tc.TestData.M, tc.TestData.t, tc.TestData.dz);
        o1 = fitZshimGradient(y1, tc.TestData.M{1}, tc.TestData.t(1), tc.TestData.dz(1));
        o2 = fitZshimGradient(y2, tc.TestData.M{2}, tc.TestData.t(2), tc.TestData.dz(2));
        ej = prctile(abs(oj.Gz(:)-G),90);  e1 = prctile(abs(o1.Gz(:)-G),90);  e2 = prctile(abs(o2.Gz(:)-G),90);
        fprintf('Gz=%g: p90 error joint %.2f, echo1 %.2f, echo2 %.2f uT/m\n', G, ej, e1, e2);
        tc.verifyLessThan(ej, 1.5);
        tc.verifyLessThanOrEqual(ej, 1.1*min(e1, e2));
    end
end
