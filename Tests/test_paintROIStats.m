function tests = test_paintROIStats
% ASE_paintROIStats on a small synthetic template: placement order, overlap dropping, NaN fill
% and the combined report.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'Aux_Mat_Functions')));
end

function setup(tc)
    d = tempname;  mkdir(d);  tc.TestData.dir = d;
    tc.addTeardown(@() rmdir(d, 's'));
    % template A: region 10 (x 1-3) and 11 (x 4-6); template B: region 20 (x 3-4, overlaps both)
    A = zeros(8,8,8,'int16');  A(1:3,:,:) = 10;  A(4:6,:,:) = 11;
    B = zeros(8,8,8,'int16');  B(3:4,:,:) = 20;  B(8,:,:) = 21;
    niftiwrite(A, fullfile(d,'tA.nii'));  niftiwrite(B, fullfile(d,'tB.nii'));
    writelines(["10	RegA";"11	RegB";"20	RegOver";"21	RegC";"40	RegAbsent"], fullfile(d,'lut.txt'));
    t = table([0.3;0.4;0.5;NaN;0.9], [3;4;5;6;9], 'VariableNames', {'OEF','R2P'}, ...
        'RowNames', {'RegA','RegB','RegOver','RegC','NotInLut'});
    writetable(t, fullfile(d,'tbl.csv'), 'WriteRowNames', true);
end

function testPaintAndReport(tc)
    d = tc.TestData.dir;
    rep = ASE_paintROIStats({fullfile(d,'tA.nii'), fullfile(d,'tB.nii')}, {fullfile(d,'lut.txt')}, ...
        fullfile(d,'tbl.csv'), fullfile(d,'painted_%STAT%'), 'stats', {'OEF','R2P','nope'}, 'verbose', false);
    tc.verifyEqual(rep.painted, ["OEF";"R2P"]);
    V = niftiread(fullfile(d,'painted_OEF.nii.gz'));
    tc.verifyEqual(double(V(2,1,1)), 0.3, 'AbsTol', 1e-6);
    tc.verifyEqual(double(V(5,1,1)), 0.4, 'AbsTol', 1e-6);
    tc.verifyTrue(isnan(V(7,1,1)));                     % outside every region -> NaN fill
    tc.verifyTrue(isnan(V(8,1,1)));                     % RegC has a NaN value
    tc.verifyFalse(any(V(:) == single(0.5)));           % RegOver overlaps RegA/RegB -> dropped whole
    s = rep.segReport;
    tc.verifyEqual(s.status(s.name == "RegOver"), "overlap");
    tc.verifyEqual(s.status(s.name == "RegAbsent"), "not in volumes");
    tc.verifyEqual(rep.notInLUT, "NotInLut");
    tc.verifyEqual(rep.nanValues.OEF, "RegC");
    txt = fileread(rep.reportFile);
    tc.verifyEqual(string(rep.reportFile), string(fullfile(d,'painted_paint_report.txt')));
    for k = ["RegOver","RegAbsent","NotInLut","RegC","nope"], tc.verifySubstring(txt, char(k)); end
end

function testLutOrderDecides(tc)
    % listing the overlapping region first places it, and the later regions it overlaps are dropped
    d = tc.TestData.dir;
    writelines(["20	RegOver";"10	RegA";"11	RegB";"21	RegC"], fullfile(d,'lut2.txt'));
    rep = ASE_paintROIStats({fullfile(d,'tA.nii'), fullfile(d,'tB.nii')}, {fullfile(d,'lut2.txt')}, ...
        fullfile(d,'tbl.csv'), fullfile(d,'p2_%STAT%'), 'stats', {'OEF'}, 'verbose', false);
    V = niftiread(fullfile(d,'p2_OEF.nii.gz'));
    tc.verifyEqual(double(V(3,1,1)), 0.5, 'AbsTol', 1e-6);
    s = rep.segReport;
    tc.verifyEqual(s.status(s.name == "RegA"), "overlap");
    tc.verifyEqual(s.status(s.name == "RegB"), "overlap");
end

function testFillValue(tc)
    d = tc.TestData.dir;
    ASE_paintROIStats({fullfile(d,'tA.nii')}, {fullfile(d,'lut.txt')}, fullfile(d,'tbl.csv'), ...
        fullfile(d,'p3_%STAT%'), 'stats', {'OEF'}, 'defaultVal', 0, 'verbose', false);
    V = niftiread(fullfile(d,'p3_OEF.nii.gz'));
    tc.verifyEqual(double(V(7,1,1)), 0);
end
