function tests = test_groupROIStats
% ASE_groupROIStats on a synthetic cohort written to a temporary folder: subject selection from
% the covariate table, fit-run choice by tag, settings check, and the per-ROI linear model.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    root = fileparts(fileparts(mfilename('fullpath')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'CoreMatFiles')));
    tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'Aux_Mat_Functions')));
end

function setup(tc)
    tc.TestData.dir = tempname;  mkdir(tc.TestData.dir);
    tc.addTeardown(@() rmdir(tc.TestData.dir, 's'));
end

%% ---- fixture helpers ----
function writeRun(root, id, runName, oef, nVox, fit)
% one fit run with a 3-ROI table: OEF given, R2P = 10*OEF, vCBV fixed
    d = fullfile(root, id, 'session_ase', 'fits', runName);  mkdir(d);
    rois = {'Left-Caudate'; 'Right-Caudate'; 'ctx-lh-insula'};
    t = table(oef(:), 10*oef(:), 0.03*ones(3,1), nVox(:), 'VariableNames', {'OEF','R2P','vCBV','nVox'}, 'RowNames', rois);
    writetable(t, fullfile(d, 'ase_fit_results_ROI_Stats.csv'), 'WriteRowNames', true);
    if nargin < 6, fit = struct('StO2_max',0.6,'tCrit_buffer',[3 12],'verbose',true); end
    writelines(jsonencode(struct('fit',fit,'B0',2.89)), fullfile(d, 'config_resolved.json'));
end

function [root, ages] = makeCohort(tc)
    root = tc.TestData.dir;
    ages = [25 35 45 55 65 75 30 60];
    slope = [0.002; -0.001; 0];          % per ROI, OEF per year
    base  = [0.30; 0.35; 0.40];
    rng(3);
    for i = 1:numel(ages)
        id = sprintf('s%02d', i);
        oef = base + slope*(ages(i) - 50) + 0.0005*randn(3,1);
        writeRun(root, id, '20260101_120000', oef + 0.5, [100 100 100]);      % older untagged run: must be ignored
        writeRun(root, id, '20260201_120000_calib', oef, [100 100 3]);       % tagged run
    end
    % subject with a tagged run but no row in the table (should be warned about)
    writeRun(root, 'x99', '20260201_120000_calib', [0.3 0.3 0.3], [100 100 100]);
    % commented-out subject with a run
    writeRun(root, 'c01', '20260201_120000_calib', [0.3 0.3 0.3], [100 100 100]);
    % listed subject with no tagged run
    writeRun(root, 'n01', '20260101_120000', [0.3 0.3 0.3], [100 100 100]);
    lines = "#subjinit sex age group";
    grp = ["ctrl","mdt","ctrl","mdt","ctrl","mdt","ctrl","AD"];
    for i = 1:numel(ages)
        lines(end+1) = sprintf('s%02d %d %d %s', i, mod(i,2)+1, ages(i), grp(i)); %#ok<AGROW>
    end
    lines(end+1) = "#c01 1 50 ctrl";
    lines(end+1) = "n01 1 50 ctrl";
    lines(end+1) = "s09 2 NA ctrl";            % no folder and missing age
    writelines(lines, fullfile(root, 'cov_table_ase'));
end

%% ---- tests ----
function testSelectionAndModel(tc)
    [root, ages] = makeCohort(tc);
    out = fullfile(root, 'out');
    res = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), out, ...
        'select', "group=ctrl,mdt", 'covariates', "age", 'fitTag', "calib", 'minN', 3, 'verbose', false);
    inc = res.subjects.subject(res.subjects.included);
    tc.verifyEqual(sort(inc), "s0" + string(1:7).');           % s08 is AD
    tc.verifyFalse(any(res.subjects.subject == "s08"));
    r = res.subjects(res.subjects.subject == "n01", :);
    tc.verifyFalse(r.included);  tc.verifySubstring(char(r.reason), 'no run tagged');
    r = res.subjects(res.subjects.subject == "s09", :);
    tc.verifyFalse(r.included);  tc.verifySubstring(char(r.reason), 'missing covariate');
    tc.verifyEqual(res.unlisted, "x99");
    % tagged runs were used (the untagged runs are offset by +0.5)
    tc.verifyTrue(all(endsWith(res.subjects.runDir(res.subjects.included), "_calib")));
    T = res.table;
    tc.verifyEqual(T{'Left-Caudate','OEF_age_beta'},  0.002, 'AbsTol', 2e-4);
    tc.verifyEqual(T{'Right-Caudate','OEF_age_beta'}, -0.001, 'AbsTol', 2e-4);
    tc.verifyLessThan(T{'Left-Caudate','OEF_age_q'}, 1e-3);
    % adjusted mean = value at the cohort's mean age
    a = mean(ages(1:7));
    tc.verifyEqual(T{'Left-Caudate','OEF_adjMean'}, 0.30 + 0.002*(a-50), 'AbsTol', 5e-4);
    tc.verifyEqual(res.centers, a, 'AbsTol', 1e-12);
    tc.verifyEqual(T{'Left-Caudate','R2P_age_beta'}, 0.02, 'AbsTol', 2e-3);
    tc.verifyEqual(T{'Left-Caudate','OEF_n'}, 7);
    tc.verifyTrue(all(isfile(res.files)));
    rep = fileread(res.files(4));
    tc.verifySubstring(rep, 'x99');  tc.verifySubstring(rep, 'c01');
end

function testMinVoxAndMinN(tc)
    root = makeCohort(tc);
    res = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'out'), ...
        'select', "group=ctrl,mdt", 'covariates', "age", 'fitTag', "calib", 'minVox', 10, 'minN', 3, 'verbose', false);
    % ctx-lh-insula has nVox = 3 in every tagged run -> no usable values
    tc.verifyEqual(res.table{'ctx-lh-insula','OEF_n'}, 0);
    tc.verifyTrue(isnan(res.table{'ctx-lh-insula','OEF_age_beta'}));
    res2 = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'out2'), ...
        'select', "group=ctrl,mdt", 'covariates', "age", 'fitTag', "calib", 'minN', 10, 'verbose', false);
    tc.verifyTrue(all(isnan(res2.table.OEF_age_beta)));
    tc.verifyEqual(res2.table{'Left-Caudate','OEF_n'}, 7);
end

function testCategoricalCovariate(tc)
    root = makeCohort(tc);
    res = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'out'), ...
        'select', ["group=ctrl,mdt","age>=30"], 'covariates', ["age","sex"], 'categorical', "sex", ...
        'fitTag', "calib", 'minN', 3, 'verbose', false);
    tc.verifyEqual(nnz(res.subjects.included), 6);              % s02-s07 (s01 is 25; n01 has no tagged run)
    tc.verifyTrue(ismember('OEF_sex_2_beta', res.table.Properties.VariableNames));
    tc.verifyEqual(res.table{'Left-Caudate','OEF_sex_2_beta'}, 0, 'AbsTol', 3e-3);
end

function testSettingsMismatch(tc)
    root = makeCohort(tc);
    f = struct('StO2_max',0.7,'tCrit_buffer',[3 12],'verbose',false);
    writeRun(root, 's03', '20260301_120000_calib', [0.3 0.35 0.4], [100 100 100], f);
    args = {'select', "group=ctrl,mdt", 'covariates', "age", 'fitTag', "calib", 'minN', 3, 'verbose', false};
    tc.verifyError(@() ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o1'), args{:}), ?MException);
    try
        ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o1b'), args{:});
    catch ME
        tc.verifySubstring(ME.message, 'fit.StO2_max differs');
        tc.verifyFalse(contains(ME.message, 'verbose'));         % output-only key is ignored
    end
    res = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o2'), args{:}, 'settingsCheck', "warn");
    tc.verifyNotEmpty(res.settingsMsg);
end

function testNoOverwrite(tc)
    root = makeCohort(tc);
    args = {'select', "group=ctrl", 'fitTag', "calib", 'minN', 3, 'verbose', false};
    ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o'), args{:});
    tc.verifyError(@() ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o'), args{:}), ?MException);
    res = ASE_groupROIStats(root, fullfile(root,'cov_table_ase'), fullfile(root,'o'), args{:}, 'overwrite', true);
    tc.verifyFalse(any(contains(res.table.Properties.VariableNames, '_beta')));   % no covariates
end

function testFindFitRun(tc)
    root = tc.TestData.dir;
    f = fullfile(root, 'fits');
    mkdir(fullfile(f, '20260101_120000'));  mkdir(fullfile(f, '20260301_120000'));
    mkdir(fullfile(f, '20260201_120000_a'));  mkdir(fullfile(f, '20260101_090000_a'));
    mkdir(fullfile(f, '20260401_120000_ab'));  mkdir(fullfile(f, 'notes'));
    tc.verifyEqual(findFitRun(f, "a"), string(fullfile(f, '20260201_120000_a')));
    tc.verifyEqual(findFitRun(f, ""), string(fullfile(f, '20260401_120000_ab')));
    [r, info] = findFitRun(f, "zz");
    tc.verifyEqual(r, "");  tc.verifySubstring(char(info.reason), 'no run tagged');
end

function testFilters(tc)
    T = table(["a";"b";"c";"d"], [1;2;NaN;4], ["x";"y";"x";missing], 'VariableNames', {'subject','age','grp'});
    tc.verifyEqual(applyCovFilters(T, "grp=x"), [true;false;true;false]);
    tc.verifyEqual(applyCovFilters(T, "grp!=x"), [false;true;false;false]);
    tc.verifyEqual(applyCovFilters(T, ["age>1","age<=4"]), [false;true;false;true]);
    tc.verifyEqual(applyCovFilters(T, "age=1,4"), [true;false;false;true]);
    tc.verifyError(@() applyCovFilters(T, "grp>1"), ?MException);
    tc.verifyError(@() applyCovFilters(T, "nope=1"), ?MException);
    tc.verifyEqual(bhFDR([0.01 0.04 NaN 0.03]), [0.03 0.04 NaN 0.04], 'AbsTol', 1e-12);
end
