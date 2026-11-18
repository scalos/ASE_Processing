function out = calibrateQuality(rootDir, opts)
% calibrateQuality  Cohort calibration of the voxel quality thresholds
% (docs/calibration-plan.md, Phase 3).
%
%   out = calibrateQuality(rootDir)
%   out = calibrateQuality(rootDir, calibSubjects=["mdh186_240824" ...], outDir="...")
%
% For every subject folder in rootDir, reads the latest fit run with each tag in runTags
% (aseFit.sh --tag; default calDiag-dB and calDiag-noDB, written with every quality diagnostic
% and writeRawMaps=true) and the subject's segmentation on the ASE side (segName). Nothing is
% refitted: thresholds are applied to the saved diagnostic maps with computeVoxelQuality itself,
% so the flags are exactly what the pipeline would produce.
%
% Tissues, on the ASE grid (occupancy from segToAseMembership, aparc+aseg labels):
%   cleanDeepWM : cerebral WM (2/41) occupancy >= 0.9, eroded by a 3x3x3 cube, CSF fraction 0
%   WM, cortex, deepGM : occupancy >= 0.5;  brain : total brain occupancy >= 0.5
%
% Steps:
%   1. per subject and run: distribution of each criterion's statistic, per tissue
%   2. per run tag: threshold = p<targetPct> (p<100-targetPct> for SNR) of clean deep WM pooled
%      over the calibration subjects, with the per-subject spread and the leave-one-subject-out
%      range. The fields in calibrateFields form the "calibrated" threshold set.
%   3. per subject, run and threshold set: fraction flagged by each criterion and combined,
%      per tissue
%
% Subjects and roles come from the covariate table (default <rootDir>/cov_table_ase):
%   healthy = group in healthyGroups;
%   calib   = calibSubjects (default: every healthy subject, i.e. no hold-out; the summary says so);
%   valid   = the other healthy subjects;  other = everyone else (reported, never used for thresholds).
%
% Threshold sets (step 3): "runDefaults" (each run's own quality block), "calibrated" (run block +
% the calibrated fields of that run tag), plus any in opts.sets: struct array with fields name and
% quality (a struct of quality-field overrides applied to the run block).
%
% Output folder (default <rootDir>/calibration/<timestamp>): subjects.csv, distributions.csv,
% thresholds.csv, flagRates.csv, summary.txt. The tables are also returned in out.
    arguments
        rootDir {mustBeTextScalar}
        opts.runTags string = ["calDiag-dB", "calDiag-noDB"]
        opts.covTable string = ""
        opts.healthyGroups string = ["ctrl", "mdt"]
        opts.subjects string = strings(0,1)
        opts.calibSubjects string = strings(0,1)
        opts.targetPct (1,1) double {mustBeInRange(opts.targetPct, 50, 100)} = 99
        opts.calibrateFields string = ["vcbvZ_max", "zshimEcho_maxDiff"]
        opts.sets struct = struct('name', {}, 'quality', {})
        opts.combinedTarget (1,1) double = 0.10
        opts.segName string = "anat_fs_aparc_aseg_al_ase.nii.gz"
        opts.outDir string = ""
        opts.verbose (1,1) logical = true
    end
    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(repo, 'CoreMatFiles'), fullfile(repo, 'Aux_Mat_Functions'));
    rootDir = string(rootDir);
    v = opts.verbose;

    %% subjects and roles
    if opts.covTable == "", opts.covTable = fullfile(rootDir, "cov_table_ase"); end
    cov = table();
    if isfile(opts.covTable), cov = readCovTable(opts.covTable);
    else, warning('calibrateQuality:noCov', 'No covariate table (%s): every subject is "other".', opts.covTable);
    end
    if isempty(opts.subjects)
        d = dir(rootDir);
        d = d([d.isdir] & ~startsWith({d.name}, '.'));
        names = string({d.name}).';
        keep = arrayfun(@(s) isfolder(fullfile(rootDir, s, "session_ase", "fits")), names);
        opts.subjects = names(keep);
    end
    subjects = opts.subjects(:);
    nS = numel(subjects);
    assert(nS > 0, 'No subject folders with session_ase/fits in %s', rootDir);

    group = strings(nS, 1);  group(:) = "NA";
    if ~isempty(cov) && ismember('group', cov.Properties.VariableNames)
        [tf, loc] = ismember(subjects, cov.subject);
        group(tf) = string(cov.group(loc(tf)));
    end
    healthy = ismember(group, opts.healthyGroups);
    noHoldout = isempty(opts.calibSubjects);
    if noHoldout
        isCalib = healthy;
    else
        isCalib = ismember(subjects, opts.calibSubjects);
        bad = opts.calibSubjects(~ismember(opts.calibSubjects, subjects(healthy)));
        if ~isempty(bad)
            warning('calibrateQuality:calib', 'Calibration subjects not found or not healthy, ignored: %s', strjoin(bad, ', '));
        end
        isCalib = isCalib & healthy;
    end
    role = repmat("other", nS, 1);
    role(healthy) = "valid";
    role(isCalib) = "calib";

    %% criteria: statistic, side of the threshold, quality field it calibrates
    crit = struct( ...
        'name',  {'vcbvZ', 'zshimEcho', 'quadSign', 'lowSNR', 'vein', 'zshimFrac', 'r2pEchoEffect'}, ...
        'stat',  {'vcbvZ', 'zshimEcho', 'kZ', 'snrMin', 'veinRatio', 'zshimFrac', 'r2pRelDiff'}, ...
        'side',  {'upper', 'upper', 'upper', 'lower', 'upper', 'upper', 'upper'}, ...
        'field', {'vcbvZ_max', 'zshimEcho_maxDiff', 'quadSign_kZmin', 'snr_min', 'vein_relK', 'zshimFrac_max', ''});
    tissues = ["cleanDeepWM", "WM", "cortex", "deepGM", "brain"];
    pcts = [1 5 50 90 95 99];
    nT = numel(opts.runTags);

    %% pass 1: distributions, clean-WM values for the thresholds
    distRows = {};
    cleanVals = cell(nS, nT);         % struct of clean-deep-WM vectors per statistic
    runDirs = strings(nS, nT);
    masks = cell(nS, 1);
    qByTag = cell(1, nT);
    subjNote = strings(nS, 1);
    veinRef = NaN(nS, nT);
    for i = 1:nS
        s = subjects(i);
        segFile = fullfile(rootDir, s, "session_ase", opts.segName);
        for t = 1:nT
            [runDir, info] = findFitRun(fullfile(rootDir, s, "session_ase", "fits"), opts.runTags(t));
            if runDir == ""
                subjNote(i) = addNote(subjNote(i), opts.runTags(t) + ": " + info.reason);
                continue
            end
            if v, fprintf('[%d/%d] %s  %s\n', i, nS, s, runDir); end
            R = loadRun(runDir);
            if ~isfield(R.qual, 'csfFrac')
                subjNote(i) = addNote(subjNote(i), opts.runTags(t) + ": no csfFrac map (quality masks missing at fit time)");
                continue
            end
            if ~R.rawUsed
                subjNote(i) = addNote(subjNote(i), opts.runTags(t) + ": no _raw maps (writeRawMaps off); bit 0 uses clipped maps");
            end
            if isempty(masks{i})
                if ~isfile(segFile)
                    subjNote(i) = addNote(subjNote(i), "segmentation not found: " + segFile);
                    break
                end
                masks{i} = tissueMasks(segFile, R.aseInfo, R.qual.csfFrac);
            end
            runDirs(i, t) = runDir;
            if isempty(qByTag{t}), qByTag{t} = R.q; end
            [~, ~, ~, qi] = computeVoxelQuality(R.raw, R.qual, R.physBds, R.q);
            veinRef(i, t) = qi.veinRef;
            S = voxelStats(R, qi.veinRef);
            statNames = fieldnames(S);
            for k = 1:numel(tissues)
                m = masks{i}.(tissues(k));
                for j = 1:numel(statNames)
                    x = S.(statNames{j})(m);
                    distRows(end+1, :) = [{s, group(i), role(i), opts.runTags(t), tissues(k), string(statNames{j}), ...
                        nnz(m), nnz(isfinite(x)), mean(x(isfinite(x)))}, num2cell(pctl(x, pcts))]; %#ok<AGROW>
                end
            end
            cleanVals{i, t} = structfun(@(x) x(masks{i}.cleanDeepWM), S, 'UniformOutput', false);
        end
    end
    distT = cell2table(distRows, 'VariableNames', ...
        [{'subject', 'group', 'role', 'run', 'tissue', 'stat', 'nVox', 'nFinite', 'mean'}, ...
         arrayfun(@(p) sprintf('p%g', p), pcts, 'UniformOutput', false)]);

    %% pass 2: thresholds per run tag (calibration subjects, clean deep WM)
    pHi = opts.targetPct;  pLo = 100 - opts.targetPct;
    thrRows = {};
    calibQ = cell(1, nT);
    for t = 1:nT
        calibQ{t} = struct();
        use = find(isCalib & ~cellfun(@isempty, cleanVals(:, t)));
        for c = 1:numel(crit)
            p = pHi;  if strcmp(crit(c).side, 'lower'), p = pLo; end
            vals = cellfun(@(z) z.(crit(c).stat), cleanVals(use, t), 'UniformOutput', false);
            perSubj = cellfun(@(x) pctl(x, p), vals);
            pooled = NaN;  nVox = 0;  loso = NaN(numel(use), 1);
            if ~isempty(use)
                allv = vertcat(vals{:});
                pooled = pctl(allv, p);  nVox = nnz(isfinite(allv));
                for k = 1:numel(use)
                    loso(k) = pctl(vertcat(vals{[1:k-1, k+1:end]}), p);
                end
            end
            cur = NaN;
            if ~isempty(crit(c).field) && ~isempty(qByTag{t}), cur = qByTag{t}.(crit(c).field); end
            applied = ~isempty(crit(c).field) && ismember(crit(c).field, opts.calibrateFields) && isfinite(pooled);
            if applied, calibQ{t}.(crit(c).field) = round(pooled, 3, 'significant'); end
            thrRows(end+1, :) = {opts.runTags(t), string(crit(c).name), string(crit(c).stat), string(crit(c).side), ...
                p, pooled, numel(use), nVox, min(perSubj), median(perSubj), max(perSubj), ...
                std(perSubj) / mean(perSubj), min(loso), max(loso), string(crit(c).field), cur, applied}; %#ok<AGROW>
        end
    end
    thrT = cell2table(thrRows, 'VariableNames', {'run', 'criterion', 'stat', 'side', 'pct', 'pooled', ...
        'nSubjects', 'nVox', 'perSubj_min', 'perSubj_median', 'perSubj_max', 'perSubj_cv', ...
        'loso_min', 'loso_max', 'field', 'current', 'inCalibratedSet'});

    %% pass 3: flag rates per subject, run and threshold set
    flagRows = {};
    for t = 1:nT
        sets = [struct('name', "runDefaults", 'quality', struct()), ...
                struct('name', "calibrated", 'quality', calibQ{t})];
        for k = 1:numel(opts.sets)
            sets(end+1) = struct('name', string(opts.sets(k).name), 'quality', opts.sets(k).quality); %#ok<AGROW>
        end
        for i = 1:nS
            if runDirs(i, t) == "", continue; end
            R = loadRun(runDirs(i, t));
            for k = 1:numel(sets)
                qs = R.q;
                f = fieldnames(sets(k).quality);
                for j = 1:numel(f), qs.(f{j}) = sets(k).quality.(f{j}); end
                qs = qualityDefaults(qs);
                [keep, flags, bitInfo] = computeVoxelQuality(R.raw, R.qual, R.physBds, qs);
                for m = 1:numel(tissues)
                    msk = masks{i}.(tissues(m));
                    n = nnz(msk);
                    for b = 1:height(bitInfo)
                        if ~bitInfo.Enabled(b), continue; end
                        fr = mean(bitand(flags(msk), uint16(2^bitInfo.Bit(b))) > 0);
                        flagRows(end+1, :) = {subjects(i), group(i), role(i), opts.runTags(t), sets(k).name, ...
                            tissues(m), n, string(bitInfo.Name{b}), fr}; %#ok<AGROW>
                    end
                    flagRows(end+1, :) = {subjects(i), group(i), role(i), opts.runTags(t), sets(k).name, ...
                        tissues(m), n, "combined", mean(~keep(msk))}; %#ok<AGROW>
                end
            end
        end
    end
    flagT = cell2table(flagRows, 'VariableNames', {'subject', 'group', 'role', 'run', 'set', 'tissue', ...
        'nVox', 'criterion', 'fracFlagged'});

    subjT = table(subjects, group, role, veinRef, subjNote, 'VariableNames', ...
        {'subject', 'group', 'role', 'veinRef', 'note'});
    for t = 1:nT
        subjT.("run_" + matlab.lang.makeValidName(opts.runTags(t))) = runDirs(:, t);
    end

    %% write
    if opts.outDir == ""
        opts.outDir = fullfile(rootDir, "calibration", string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    end
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end
    writetable(subjT, fullfile(opts.outDir, 'subjects.csv'));
    writetable(distT, fullfile(opts.outDir, 'distributions.csv'));
    writetable(thrT, fullfile(opts.outDir, 'thresholds.csv'));
    writetable(flagT, fullfile(opts.outDir, 'flagRates.csv'));
    txt = summaryText(subjT, thrT, flagT, opts, noHoldout);
    fid = fopen(fullfile(opts.outDir, 'summary.txt'), 'w');  fprintf(fid, '%s', txt);  fclose(fid);
    if v, fprintf('%s\nWritten to %s\n', txt, opts.outDir); end

    out = struct('subjects', subjT, 'distributions', distT, 'thresholds', thrT, ...
                 'flagRates', flagT, 'calibratedQuality', {calibQ}, 'outDir', opts.outDir);
end

% ---------------------------------------------------------------------------------------------
function R = loadRun(runDir)
% Raw maps, quality diagnostics, bounds and quality options of one fit run (files looked up by
% name in the run folder: the manifest's absolute paths may be from another machine).
    cfg = jsondecode(fileread(fullfile(runDir, 'config_resolved.json')));
    fit = cfg.fit;
    R.q = qualityDefaults(fit.quality);
    for nm = ["R2", "R2P", "vCBV", "OEF"]
        R.physBds.(nm) = toNum(fit.("physBd_" + nm));
    end
    R.rawUsed = true;
    for nm = ["OEF", "R2P", "vCBV", "R2"]
        f = mapFile(runDir, nm + "_raw");
        if f == "", f = mapFile(runDir, nm); R.rawUsed = false; end
        assert(f ~= "", 'Map %s not found in %s', nm, runDir);
        R.raw.(nm) = double(niftiread(f));
        if nm == "OEF", R.aseInfo = niftiinfo(f); end
    end
    R.qual = struct();
    for nm = ["SNR", "Gz", "zshimEdge"]          % per echo -> 4th dimension
        e = 1;  x = [];
        while true
            f = mapFile(runDir, sprintf('%s_e%d', nm, e));
            if f == "", break; end
            x = cat(4, x, double(niftiread(f)));  e = e + 1;
        end
        if ~isempty(x), R.qual.(nm) = x; end
    end
    for nm = ["csfFrac", "zshimFracLow", "vCBV_quad", "vCBV_z", "kZ", "R2P_e1", "R2P_e2", "R2P_Fstat"]
        f = mapFile(runDir, nm);
        if f ~= "", R.qual.(nm) = double(niftiread(f)); end
    end
end

function f = mapFile(runDir, name)
    d = dir(fullfile(runDir, "*_" + name + ".nii*"));
    d = d(~startsWith({d.name}, '._'));
    f = "";
    if isscalar(d), f = string(fullfile(d.folder, d.name)); end
end

function x = toNum(x)
% JSON bounds may hold "Inf" / "-Inf" strings (mixed arrays decode to cell)
    if iscell(x)
        x = cellfun(@(z) ifelse(ischar(z) || isstring(z), str2double(z), z), x);
    end
    x = double(x(:).');
end

function y = ifelse(c, a, b)
    if c, y = a; else, y = b; end
end

function M = tissueMasks(segFile, aseInfo, csf)
% Tissue masks on the ASE grid from the aparc+aseg labels (high-res grid in ASE space).
    segInfo = niftiinfo(segFile);
    seg = round(double(niftiread(segInfo)));
    wm     = [2 41];
    deepGM = [10 11 12 13 17 18 26 28 49 50 51 52 53 54 58 60];
    csfLab = [4 5 14 15 24 31 43 44 63 72];
    cls = zeros(size(seg));
    cls(ismember(seg, wm)) = 1;
    cls((seg >= 1000 & seg < 3000) | ismember(seg, [3 42])) = 2;
    cls(ismember(seg, deepGM)) = 3;
    cls(seg > 0 & cls == 0 & ~ismember(seg, csfLab)) = 4;     % other brain (cerebellum, stem, ...)
    LUT = table((1:4).', ["WM"; "cortex"; "deepGM"; "otherBrain"], 'VariableNames', {'code', 'name'});
    Mem = segToAseMembership(cls, LUT, segInfo, aseInfo, "minFrac", 0);
    sz = aseInfo.ImageSize(1:3);
    occ = reshape(full(Mem).', [sz 4]);
    pureWM = double(occ(:,:,:,1) >= 0.9);
    eroded = convn(pureWM, ones(3,3,3), 'same') >= 27 - 0.5;
    M.cleanDeepWM = eroded & isfinite(csf) & csf <= 1e-6;
    M.WM     = occ(:,:,:,1) >= 0.5;
    M.cortex = occ(:,:,:,2) >= 0.5;
    M.deepGM = occ(:,:,:,3) >= 0.5;
    M.brain  = sum(occ, 4) >= 0.5;
end

function S = voxelStats(R, veinRef)
% The statistic each criterion thresholds (NaN where the diagnostic is missing).
    sz = size(R.raw.vCBV);  nanv = NaN(sz);
    g = @(f) getOr(R.qual, f, nanv);
    S.vcbvZ = abs(g('vCBV_z'));
    Gz = g('Gz');
    if size(Gz, 4) >= 2, S.zshimEcho = abs(Gz(:,:,:,1) - Gz(:,:,:,2)); else, S.zshimEcho = nanv; end
    S.kZ = g('kZ');
    S.snrMin = min(g('SNR'), [], 4);
    S.veinRatio = R.raw.vCBV ./ veinRef;
    S.csf = g('csfFrac');
    S.zshimFrac = g('zshimFracLow');
    S.r2pF = g('R2P_Fstat');
    S.r2pRelDiff = abs(g('R2P_e1') - g('R2P_e2')) ./ abs(R.raw.R2P);
end

function x = getOr(s, f, d)
    if isfield(s, f), x = s.(f); else, x = d; end
end

function v = pctl(x, p)
% Percentiles of the finite values (prctile's convention: x(i) sits at 100*(i-0.5)/n).
    x = sort(x(isfinite(x)));  n = numel(x);
    if n == 0, v = NaN(size(p)); return; end
    if n == 1, v = repmat(x, size(p)); return; end
    pos = min(max(p / 100 * n + 0.5, 1), n);
    v = interp1(1:n, x, pos);
end

function s = addNote(s, msg)
    if s == "", s = msg; else, s = s + "; " + msg; end
end

function txt = summaryText(subjT, thrT, flagT, opts, noHoldout)
    L = strings(0,1);
    L(end+1) = "calibrateQuality summary  " + string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
    L(end+1) = sprintf("Subjects: %d (calib %d, valid %d, other %d)", height(subjT), ...
        nnz(subjT.role == "calib"), nnz(subjT.role == "valid"), nnz(subjT.role == "other"));
    evalRole = "valid";
    if noHoldout || ~any(subjT.role == "valid")
        evalRole = "calib";
        L(end+1) = "NO HOLD-OUT: thresholds are calibrated and evaluated on the same subjects (optimistic).";
    end
    bad = subjT(subjT.note ~= "", :);
    for k = 1:height(bad), L(end+1) = "  note " + bad.subject(k) + ": " + bad.note(k); end %#ok<AGROW>
    for r = unique(thrT.run, 'stable').'
        L(end+1) = "";  L(end+1) = "== run " + r + " ==";
        L(end+1) = sprintf("%-14s %6s %10s %10s %24s %20s", "criterion", "pct", "pooled", "current", ...
            "per-subject min/med/max", "leave-one-out range");
        T = thrT(thrT.run == r, :);
        for k = 1:height(T)
            mark = "";  if T.inCalibratedSet(k), mark = "  <- calibrated set"; end
            L(end+1) = sprintf("%-14s %6g %10.4g %10.4g %8.3g/%7.3g/%7.3g %9.4g - %-9.4g%s", T.criterion(k), ...
                T.pct(k), T.pooled(k), T.current(k), T.perSubj_min(k), T.perSubj_median(k), ...
                T.perSubj_max(k), T.loso_min(k), T.loso_max(k), mark); %#ok<AGROW>
        end
        F = flagT(flagT.run == r & flagT.role == evalRole & flagT.tissue == "cleanDeepWM" ...
                  & flagT.criterion == "combined", :);
        for st = unique(F.set, 'stable').'
            x = F.fracFlagged(F.set == st);
            if isempty(x), continue; end
            L(end+1) = sprintf("set %-12s combined false exclusion in clean deep WM (%s, n=%d): median %.1f%%, max %.1f%%, %d over %.0f%%", ...
                st, evalRole, numel(x), 100*median(x), 100*max(x), nnz(x > opts.combinedTarget), 100*opts.combinedTarget); %#ok<AGROW>
        end
    end
    txt = strjoin(L, newline) + newline;
end
