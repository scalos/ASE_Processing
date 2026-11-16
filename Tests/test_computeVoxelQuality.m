function tests = test_computeVoxelQuality
% Each criterion sets only its own bit; disabled criteria never set bits.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'Aux_Mat_Functions')));
end

function [raw, qual, bds] = cleanInputs()
    sz = [11 1 1];                                  % voxel i (2..11) breaks criterion i-2
    raw = struct('R2',12*ones(sz),'R2P',5*ones(sz),'vCBV',0.04*ones(sz),'OEF',0.3*ones(sz));
    qual = struct('SNR',100*ones([sz 2]),'csfFrac',zeros(sz),'Gz',zeros([sz 2]), ...
                  'zshimEdge',zeros([sz 2]),'zshimFracLow',zeros(sz),'vCBV_quad',0.04*ones(sz), ...
                  'vCBV_z',zeros(sz),'kZ',-5*ones(sz),'R2P_pval',0.5*ones(sz));
    bds = struct('R2',[0 Inf],'R2P',[0 Inf],'vCBV',[0 1],'OEF',[0 1]);
    raw.OEF(2) = 1.5;              % bit 0
    qual.kZ(3) = 3;                % bit 1 (k significantly >= 0)
    qual.vCBV_z(4) = -4;           % bit 2
    qual.SNR(5,1,1,2) = 5;         % bit 3 (echo 2 only)
    qual.csfFrac(6) = 0.5;         % bit 4
    raw.vCBV(7) = 0.3;             % bit 5
    qual.zshimEdge(8,1,1,1) = 2;   % bit 6
    qual.Gz(9,1,1,2) = 50;         % bit 7
    qual.zshimFracLow(10) = 0.5;   % bit 8
    qual.R2P_pval(11) = 1e-4;      % bit 9
end

function testEachBitIsolated(tc)
    [raw, qual, bds] = cleanInputs();
    q = qualityDefaults(struct('vein_mode','absolute','enabled',true,'useR2PEcho',true));   % bit 9 is off by default
    [keep, flags, info] = computeVoxelQuality(raw, qual, bds, q);
    tc.verifyEqual(flags(1), uint16(0));
    tc.verifyTrue(keep(1));
    for b = 0:9
        tc.verifyEqual(flags(b+2), uint16(2^b), sprintf('bit %d', b));
    end
    tc.verifyEqual(nnz(keep), 1);
    tc.verifyEqual(info.nFlagged, ones(10,1));
end

function testDisabledNeverSets(tc)
    [raw, qual, bds] = cleanInputs();
    names = {'useBounds','useQuadSign','useVcbvZ','useSNR','useCSF','useVein', ...
             'useZshimEdge','useZshimEcho','useZshimFrac','useR2PEcho'};
    for b = 0:9
        q = qualityDefaults(struct('vein_mode','absolute','enabled',true, names{b+1}, false));
        [keep, flags] = computeVoxelQuality(raw, qual, bds, q);
        tc.verifyEqual(nnz(bitand(flags, uint16(2^b))), 0, sprintf('bit %d set while disabled', b));
        tc.verifyTrue(keep(b+2), sprintf('voxel for bit %d should be kept', b));
    end
end

function testMissingInputWarnsAndSkips(tc)
    [raw, qual, bds] = cleanInputs();
    qual = rmfield(qual, 'csfFrac');
    q = qualityDefaults(struct('vein_mode','absolute','enabled',true));
    [~, flags, info] = tc.verifyWarning(@() computeVoxelQuality(raw, qual, bds, q), 'computeVoxelQuality:missingInput');
    tc.verifyFalse(info.Available(5));
    tc.verifyEqual(nnz(bitand(flags, uint16(16))), 0);
end

function testNaNFitIsImplausible(tc)
    [raw, qual, bds] = cleanInputs();
    raw.R2P(1) = NaN;
    [~, flags] = computeVoxelQuality(raw, qual, bds, qualityDefaults(struct('vein_mode','absolute','enabled',true)));
    tc.verifyEqual(flags(1), uint16(1));
end

function testVeinRelativeThreshold(tc)
% threshold = vein_relK x median vCBV over tissue voxels with plausible fits
    rng(3);
    n = 400;  sz = [n 1 1];
    raw = struct('R2',12*ones(sz),'R2P',5*ones(sz),'vCBV',0.05 + 0.02*rand(sz),'OEF',0.3*ones(sz));
    qual = struct('csfFrac', zeros(sz));
    bds = struct('R2',[0 Inf],'R2P',[0 Inf],'vCBV',[0 1],'OEF',[0 1]);
    qual.csfFrac(1:100) = 0.5;  raw.vCBV(1:100) = 0.9;     % CSF voxels: excluded from the reference
    raw.OEF(101:120) = 2;  raw.vCBV(101:120) = 0.8;        % implausible: excluded from the reference
    raw.vCBV(200) = 0.5;                                   % a vein
    q = qualityDefaults(struct('enabled',true,'vein_relK',4));
    [~, flags, ~, info] = computeVoxelQuality(raw, qual, bds, q);
    ref = raw.vCBV(121:end);  ref = ref(:);
    tc.verifyEqual(info.veinRef, median(ref), 'AbsTol', 1e-12);
    tc.verifyEqual(info.veinThreshold, 4*median(ref), 'AbsTol', 1e-12);
    tc.verifyEqual(info.veinRefN, numel(ref));
    isVein = bitand(flags, uint16(32)) > 0;
    tc.verifyEqual(find(isVein & (1:n)' > 120), 200);      % among non-CSF, plausible voxels only the vein
    tc.verifyTrue(all(isVein(1:100)));                     % 0.9 > threshold, flagged regardless
end

function testVeinRelativeNeedsCSF(tc)
    [raw, qual, bds] = cleanInputs();
    qual = rmfield(qual, 'csfFrac');
    q = qualityDefaults(struct('enabled',true,'useCSF',false));
    [~, flags, info] = tc.verifyWarning(@() computeVoxelQuality(raw, qual, bds, q), 'computeVoxelQuality:missingInput');
    tc.verifyFalse(info.Available(6));
    tc.verifyEqual(nnz(bitand(flags, uint16(32))), 0);
end

function testVeinModeValidated(tc)
    tc.verifyError(@() qualityDefaults(struct('vein_mode','percentile')), ?MException);
end
