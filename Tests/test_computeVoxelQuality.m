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
                  'vCBV_z',zeros(sz),'R2P_pval',0.5*ones(sz));
    bds = struct('R2',[0 Inf],'R2P',[0 Inf],'vCBV',[0 1],'OEF',[0 1]);
    raw.OEF(2) = 1.5;              % bit 0
    qual.vCBV_quad(3) = NaN;       % bit 1 (k >= 0)
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
    q = qualityDefaults(struct('enabled',true));
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
        q = qualityDefaults(struct('enabled',true, names{b+1}, false));
        [keep, flags] = computeVoxelQuality(raw, qual, bds, q);
        tc.verifyEqual(nnz(bitand(flags, uint16(2^b))), 0, sprintf('bit %d set while disabled', b));
        tc.verifyTrue(keep(b+2), sprintf('voxel for bit %d should be kept', b));
    end
end

function testMissingInputWarnsAndSkips(tc)
    [raw, qual, bds] = cleanInputs();
    qual = rmfield(qual, 'csfFrac');
    q = qualityDefaults(struct('enabled',true));
    [~, flags, info] = tc.verifyWarning(@() computeVoxelQuality(raw, qual, bds, q), 'computeVoxelQuality:missingInput');
    tc.verifyFalse(info.Available(5));
    tc.verifyEqual(nnz(bitand(flags, uint16(16))), 0);
end

function testNaNFitIsImplausible(tc)
    [raw, qual, bds] = cleanInputs();
    raw.R2P(1) = NaN;
    [~, flags] = computeVoxelQuality(raw, qual, bds, qualityDefaults(struct('enabled',true)));
    tc.verifyEqual(flags(1), uint16(1));
end
