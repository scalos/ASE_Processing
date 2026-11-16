function tests = test_estimateRepNoise
% Unit tests for estimateRepNoise. Run: results = runtests('Tests')
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Aux_Mat_Functions')));
end

function testTwoRepsWithDrift(tc)
    rng(1);
    nSweep = 45; nReps = 2; sz = [20 20 5];
    repInds = reshape(1:nSweep*nReps, nSweep, nReps).';
    taus = linspace(0,0.044,nSweep);
    sigTrue = 7;
    clean = 1000*exp(-20*taus);                              % same signal each rep
    data = zeros([sz, nSweep*nReps + 8]);
    drift = [0 25];                                          % rep 2 offset (drift)
    for r = 1:nReps
        data(:,:,:,repInds(r,:)) = reshape(clean + drift(r), 1,1,1,[]) + sigTrue*randn([sz nSweep]);
    end
    data(:,:,:,end-7:end) = 500;                             % z-shim frames: ignored
    [sigma, S0] = estimateRepNoise(data, repInds, 'baseInds', 1);
    tc.verifyEqual(median(sigma(:)), sigTrue, 'RelTol', 0.03);
    tc.verifyEqual(median(S0(:)), 1000 + mean(drift), 'RelTol', 0.01);
end

function testThreeRepsMatchesTruth(tc)
    rng(2);
    nSweep = 30; nReps = 3; sz = [30 30 4];
    repInds = reshape(1:nSweep*nReps, nSweep, nReps).';
    data = zeros([sz, nSweep*nReps]);
    for r = 1:nReps
        data(:,:,:,repInds(r,:)) = 300 + 10*r + 3*randn([sz nSweep]);
    end
    sigma = estimateRepNoise(data, repInds);
    tc.verifyEqual(median(sigma(:)), 3, 'RelTol', 0.03);
end

function testTwoRepFormulaEquivalence(tc)
    % for 2 reps the pooled form equals sum((d-mean d)^2)/(2(n-1))
    rng(3);
    nSweep = 45; repInds = [1:45; 46:90];
    data = reshape(100 + randn(1,90), 1,1,1,[]);
    sigma = estimateRepNoise(data, repInds);
    d = squeeze(data(1,1,1,1:45) - data(1,1,1,46:90));
    ref = sqrt(sum((d-mean(d)).^2)/(2*(nSweep-1)));
    tc.verifyEqual(sigma, ref, 'AbsTol', 1e-10);
end

function testSingleRepWarnsNaN(tc)
    data = rand(4,4,2,10);
    sigma = tc.verifyWarning(@() estimateRepNoise(data, 1:10), 'estimateRepNoise:tooFew');
    tc.verifyTrue(all(isnan(sigma(:))));
end
