function tests = test_fitZshimGradient
% Unit tests for fitZshimGradient / zshimAttenuation. Run: runtests('Tests')
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
    % test-protocol values (bch046 echo 1): moments in T*s/m, 44 ms, 4.8 mm
    tc.TestData.M  = [3819 3055 2291 1527 763 0 -763 -1527].'*1e-9;
    tc.TestData.t  = 0.044;
    tc.TestData.dz = 4.8e-3;
end

function y = synth(tc, Gz, A, snr, n)
    gamma = 2.67502e8;
    x = gamma*tc.TestData.dz*(Gz*1e-6*tc.TestData.t - tc.TestData.M.')/2;
    f = abs(sin(x)./x); f(x==0) = 1;
    y = A*f + (A/snr)*randn(n, numel(tc.TestData.M));     % n x nShim
    y = reshape(y, n,1,1,[]);
end

function testRecoversGzRealisticNoise(tc)
    % SNR 80 at the peak (median tissue SNR of the test subject is ~77-102)
    rng(10);
    for GzTrue = [-30 -10 0 12.3 40 75]
        y = synth(tc, GzTrue, 500, 80, 2000);
        o = fitZshimGradient(y, tc.TestData.M, tc.TestData.t, tc.TestData.dz);
        err = abs(o.Gz(:) - GzTrue);
        tc.verifyLessThan(median(err), 2, sprintf('Gz=%g: median error %.2f', GzTrue, median(err)));
        tc.verifyLessThan(prctile(err,90), 2, sprintf('Gz=%g: p90 error %.2f', GzTrue, prctile(err,90)));
        tc.verifyEqual(median(o.A(:)), 500, 'RelTol', 0.02);
    end
end

function testNoiselessExact(tc)
    y = synth(tc, 23.4, 300, Inf, 1);
    o = fitZshimGradient(y, tc.TestData.M, tc.TestData.t, tc.TestData.dz, 'GzStep', 0.5);
    tc.verifyEqual(o.Gz, 23.4, 'AbsTol', 0.1);
    tc.verifyEqual(o.edge, uint8(0));
    F0 = zshimAttenuation(o.Gz, tc.TestData.t, tc.TestData.dz);
    tc.verifyEqual(o.atten, F0, 'AbsTol', 1e-3);   % atten ratio = model attenuation at M=0
end

function testEdgeFlags(tc)
    % peak outside the sampled moment range (> 3.819e-6/0.044 = 86.8 uT/m)
    y = synth(tc, 110, 300, Inf, 1);
    o = fitZshimGradient(y, tc.TestData.M, tc.TestData.t, tc.TestData.dz);
    tc.verifyTrue(bitand(o.edge, uint8(2)) > 0);
    % true Gz beyond the search grid -> grid-edge bit
    y = synth(tc, 80, 300, Inf, 1);
    o = fitZshimGradient(y, tc.TestData.M, tc.TestData.t, tc.TestData.dz, 'GzRange', [-50 50]);
    tc.verifyTrue(bitand(o.edge, uint8(1)) > 0);
end

function testBackgroundNaN(tc)
    o = fitZshimGradient(zeros(3,3,1,8), tc.TestData.M, tc.TestData.t, tc.TestData.dz);
    tc.verifyTrue(all(isnan(o.Gz(:))));
end

function testAttenuationUsesUnnormalizedSinc(tc)
    gamma = 2.67502e8; dz = 4.8e-3; t = 0.044;
    Gz = pi/(gamma*dz*t/2)*1e6;          % first null of sin(x)/x
    tc.verifyEqual(zshimAttenuation(Gz, t, dz), 0, 'AbsTol', 1e-12);
    tc.verifyEqual(zshimAttenuation(0, t, dz), 1);
end
