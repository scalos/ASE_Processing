function res = ASE_groupROIStats(allSubjDir, covTableFile, outDir, opts)
% Group statistics on per-subject ROI tables (ase_fit_results_ROI_Stats.csv), one linear model
% per ROI and stat:   value ~ 1 + covariates.
%
% Subjects: the covariate table defines who is in the analysis. Every uncommented row that
% passes opts.select is used if its fit can be found. Problems are reported, not fatal:
%   - a selected subject with no folder, no matching fit run, no ROI table, or a missing
%     covariate value is left out (reason in group_subjects.csv and the report);
%   - a subject folder in allSubjDir with a matching fit run but no row in the table is
%     listed as "not in covariate table" (warning), since it may have been left out by mistake.
%
% Fit run per subject: findFitRun(<subj>/<fitsSubdir>, opts.fitTag), i.e. the latest run
% tagged opts.fitTag (aseFit.sh --tag), or the latest run of any kind when no tag is given.
% The "fit" sections of the chosen runs' config_resolved.json must match (opts.settingsCheck),
% ignoring output-only keys (opts.ignoreFitKeys).
%
% Model: numeric covariates are centred on the mean of the subjects in the analysis, so
% adjMean is the ROI mean at the cohort's mean covariate values (averaged over the categorical
% covariates' observed mix). Text columns, and numeric columns named in opts.categorical, are
% categorical (reference level = first in sorted order). An ROI value is used when it is finite
% and its nVox >= opts.minVox. ROIs with fewer than max(opts.minN, nCoefficients+2) usable
% subjects get NaN results. q = Benjamini-Hochberg FDR across ROIs, per stat and term.
%
% Outputs in outDir:
%   group_ROI_Stats.csv  one row per ROI (row names = ROI names, like the per-subject tables);
%                        columns <stat>_n, _mean, _sd, _adjMean, and <stat>_<term>_beta/_se/_t/_p/_q
%                        for each model term. Can be painted with ASE_paintROIStats.
%   group_ROI_long.csv   every value used: subject x ROI x stat, with nVox and the covariates
%   group_subjects.csv   every selected subject: included or not (reason), fit run, covariates
%   group_report.txt     inputs, cohort, exclusions, settings check, centring values, ROI counts
    arguments
        allSubjDir {mustBeTextScalar}
        covTableFile {mustBeTextScalar}
        outDir {mustBeTextScalar}
        opts.select = strings(0,1)            % filters, see applyCovFilters (ANDed)
        opts.covariates = strings(0,1)        % covariate table columns used as covariates
        opts.categorical = strings(0,1)       % numeric columns to treat as categorical (e.g. sex)
        opts.stats = ["OEF","R2P","vCBV"]
        opts.fitTag = ""
        opts.fitsSubdir = "session_ase/fits"  % relative to each subject folder
        opts.tableName = "ase_fit_results_ROI_Stats.csv"
        opts.minVox (1,1) double = 0
        opts.minN (1,1) double = 5
        opts.settingsCheck {mustBeMember(opts.settingsCheck,["error","warn","off"])} = "error"
        opts.ignoreFitKeys = ["verbose","writeOpts","overwriteOpts","writeRawMaps"]
        opts.overwrite (1,1) logical = false
        opts.verbose (1,1) logical = true
    end
    allSubjDir = string(allSubjDir);  outDir = string(outDir);
    covs = toStrings(opts.covariates);  catCols = toStrings(opts.categorical);
    stats = toStrings(opts.stats);
    select = string(opts.select);  select = select(:);   % filters contain commas: never split
    v = opts.verbose;
    assert(isfolder(allSubjDir), 'allSubjDir not found: %s', allSubjDir);
    assert(~isempty(stats), 'No stats requested.');
    outFiles = fullfile(outDir, ["group_ROI_Stats.csv","group_ROI_long.csv","group_subjects.csv","group_report.txt"]);
    if ~opts.overwrite && any(isfile(outFiles))
        error('Output already exists in %s (set overwrite to replace it).', outDir);
    end
    report = strings(0,1);
    addR = @(varargin) sprintf(varargin{:});

    %% ---- covariate table and selection ----
    [C, cinfo] = readCovTable(covTableFile);
    [keep, fdesc] = applyCovFilters(C, select);
    vn = string(C.Properties.VariableNames);
    covVar = strings(size(covs));
    for k = 1:numel(covs)
        hit = find(strcmpi(vn, covs(k)) | strcmpi(vn, matlab.lang.makeValidName(covs(k))), 1);
        assert(~isempty(hit), 'Covariate "%s" is not a column of %s (columns: %s).', ...
            covs(k), covTableFile, strjoin(vn(2:end), ', '));
        assert(hit > 1, 'The subject ID column cannot be a covariate.');
        covVar(k) = vn(hit);
    end
    for k = 1:numel(catCols)
        hit = find(strcmpi(vn, catCols(k)) | strcmpi(vn, matlab.lang.makeValidName(catCols(k))), 1);
        assert(~isempty(hit), 'Categorical column "%s" is not in %s.', catCols(k), covTableFile);
        catCols(k) = vn(hit);
    end
    isCat = ismember(covVar, catCols) | arrayfun(@(c) ~isnumeric(C.(c)), covVar);

    S = C(keep, :);
    nSel = height(S);
    subj = struct('id',S.subject, 'included',true(nSel,1), 'reason',strings(nSel,1), 'runDir',strings(nSel,1));
    if v, fprintf('%d of %d subject(s) in %s pass the selection.\n', nSel, height(C), covTableFile); end

    % missing covariate values
    for k = 1:numel(covVar)
        miss = ismissing(S.(covVar(k)));
        for i = find(miss & subj.included).'
            subj.included(i) = false;  subj.reason(i) = "missing covariate " + covVar(k);
        end
    end

    %% ---- find fit runs and read tables ----
    tbls = cell(nSel,1);  fitCfg = cell(nSel,1);
    for i = 1:nSel
        if ~subj.included(i), continue; end
        sdir = fullfile(allSubjDir, subj.id(i));
        if ~isfolder(sdir)
            subj.included(i) = false;  subj.reason(i) = "no subject folder";  continue
        end
        [runDir, rinfo] = findFitRun(fullfile(sdir, opts.fitsSubdir), opts.fitTag);
        if strlength(runDir) == 0
            subj.included(i) = false;  subj.reason(i) = rinfo.reason;  continue
        end
        subj.runDir(i) = runDir;
        tf = fullfile(runDir, opts.tableName);
        if ~isfile(tf)
            subj.included(i) = false;  subj.reason(i) = "no " + opts.tableName + " in " + runDir;  continue
        end
        t = readtable(tf, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
        missStats = setdiff(stats, string(t.Properties.VariableNames), 'stable');
        if numel(missStats) == numel(stats)
            subj.included(i) = false;  subj.reason(i) = "ROI table has none of the stats";  continue
        end
        if opts.minVox > 0 && ~ismember("nVox", string(t.Properties.VariableNames))
            error('minVox is set but %s has no nVox column.', tf);
        end
        if ~isempty(missStats)
            report(end+1,1) = addR('Warning: %s: ROI table lacks %s (treated as missing).', ...
                subj.id(i), strjoin(missStats, ', ')); %#ok<AGROW>
        end
        tbls{i} = t;
        cf = fullfile(runDir, 'config_resolved.json');
        if isfile(cf)
            c = jsondecode(fileread(cf));
            if isfield(c, 'fit'), fitCfg{i} = c.fit; end
        end
    end

    %% ---- subjects with fits that are not in the table ----
    listed = [C.subject; cinfo.commented];
    d = dir(allSubjDir);
    d = d([d.isdir] & ~startsWith({d.name}, '.'));
    unlisted = strings(0,1);
    for k = 1:numel(d)
        id = string(d(k).name);
        if ismember(id, listed), continue; end
        if strlength(findFitRun(fullfile(allSubjDir, id, opts.fitsSubdir), opts.fitTag)) > 0
            unlisted(end+1,1) = id; %#ok<AGROW>
        end
    end
    commentedWithFit = strings(0,1);
    for id = cinfo.commented.'
        if isfolder(fullfile(allSubjDir, id)) && ...
                strlength(findFitRun(fullfile(allSubjDir, id, opts.fitsSubdir), opts.fitTag)) > 0
            commentedWithFit(end+1,1) = id; %#ok<AGROW>
        end
    end

    %% ---- settings check ----
    inc = find(subj.included);
    settingsMsg = strings(0,1);
    if opts.settingsCheck ~= "off" && ~isempty(inc)
        noCfg = inc(cellfun(@isempty, fitCfg(inc)));
        if ~isempty(noCfg)
            settingsMsg(end+1,1) = "no config_resolved.json (fit section) for: " + strjoin(subj.id(noCfg), ', ');
        end
        have = setdiff(inc, noCfg, 'stable');
        if numel(have) > 1
            keys = strings(0,1);
            for i = have.', keys = union(keys, string(fieldnames(fitCfg{i}))); end
            keys = setdiff(keys, toStrings(opts.ignoreFitKeys));
            for key = keys.'
                enc = strings(numel(have),1);
                for j = 1:numel(have)
                    f = fitCfg{have(j)};
                    if isfield(f, key), enc(j) = string(jsonencode(f.(key))); else, enc(j) = "<absent>"; end
                end
                [u,~,g] = unique(enc);
                if numel(u) > 1
                    parts = strings(numel(u),1);
                    for m = 1:numel(u)
                        parts(m) = u(m) + " (" + strjoin(subj.id(have(g==m)), ', ') + ")";
                    end
                    settingsMsg(end+1,1) = "fit." + key + " differs: " + strjoin(parts, '; '); %#ok<AGROW>
                end
            end
        end
        if ~isempty(settingsMsg)
            msg = sprintf('Fit settings differ between subjects:\n  %s', strjoin(settingsMsg, newline + "  "));
            if opts.settingsCheck == "error"
                error('%s\n(Refit with the same settings, choose runs with --fitTag, or set settingsCheck to "warn".)', msg);
            else
                warning('%s', msg);
            end
        end
    end

    %% ---- long table ----
    inc = find(subj.included);
    assert(~isempty(inc), 'No subjects left to analyse. Reasons:\n  %s', ...
        strjoin(subj.id + ": " + subj.reason, newline + "  "));
    roiNames = strings(0,1);
    for i = inc.', roiNames = union(roiNames, string(tbls{i}.Properties.RowNames), 'stable'); end
    L = cell(numel(inc),1);
    for j = 1:numel(inc)
        i = inc(j);  t = tbls{i};
        rn = string(t.Properties.RowNames);
        nv = NaN(numel(rn),1);
        if ismember("nVox", string(t.Properties.VariableNames)), nv = t.nVox; end
        parts = cell(numel(stats),1);
        for s = 1:numel(stats)
            if ismember(stats(s), string(t.Properties.VariableNames)), val = t.(stats(s)); else, val = NaN(numel(rn),1); end
            parts{s} = table(repmat(subj.id(i),numel(rn),1), rn, repmat(stats(s),numel(rn),1), val, nv, ...
                'VariableNames', {'subject','roi','stat','value','nVox'});
        end
        L{j} = vertcat(parts{:});
    end
    L = vertcat(L{:});
    Ccov = S(inc, ["subject", covVar(:).']);
    L = join(L, Ccov, 'Keys', 'subject');

    %% ---- covariates: centring and categorical coding ----
    centers = NaN(numel(covVar),1);
    for k = 1:numel(covVar)
        if isCat(k)
            L.(covVar(k)) = removecats(categorical(string(L.(covVar(k)))));
        else
            centers(k) = mean(Ccov.(covVar(k)));
            L.(covVar(k)) = L.(covVar(k)) - centers(k);
        end
    end

    %% ---- per-ROI models ----
    % expected coefficient names (fitlm naming): numeric -> name; categorical -> name_level for
    % each non-reference level present in the cohort
    expTerms = strings(0,1);
    for k = 1:numel(covVar)
        if isCat(k)
            lv = string(categories(L.(covVar(k))));
            expTerms = [expTerms; matlab.lang.makeValidName(covVar(k) + "_" + lv(2:end))]; %#ok<AGROW>
        else
            expTerms(end+1,1) = covVar(k); %#ok<AGROW>
        end
    end
    formula = "value ~ 1";
    if ~isempty(covVar), formula = formula + " + " + strjoin(covVar, " + "); end
    nR = numel(roiNames);
    out = table('RowNames', cellstr(roiNames));
    termCols = struct();       % per stat: list of term names seen
    roiNotes = strings(0,1);
    ws = warning('off', 'stats:LinearModel:RankDefDesignMat');
    cleanupWarn = onCleanup(@() warning(ws));
    for s = 1:numel(stats)
        st = stats(s);
        n = zeros(nR,1);  mu = NaN(nR,1);  sd = NaN(nR,1);  adj = NaN(nR,1);
        coef = struct();       % term -> [beta se t p] per ROI; every expected term gets a column
        for tn = expTerms.', coef.(tn) = NaN(nR,4); end
        Ls = L(L.stat == st, :);
        for r = 1:nR
            d = Ls(Ls.roi == roiNames(r), :);
            use = isfinite(d.value);
            if opts.minVox > 0, use = use & d.nVox >= opts.minVox; end
            d = d(use, :);
            n(r) = height(d);
            if n(r) > 0, mu(r) = mean(d.value); end
            if n(r) > 1, sd(r) = std(d.value); end
            nCoef = 1;
            for k = 1:numel(covVar)
                if isCat(k), nCoef = nCoef + max(numel(categories(removecats(d.(covVar(k))))) - 1, 0);
                else, nCoef = nCoef + 1; end
            end
            if n(r) < max(opts.minN, nCoef + 2)
                if n(r) > 0
                    roiNotes(end+1,1) = sprintf('%s / %s: n = %d, below the minimum %d (no model).', ...
                        st, roiNames(r), n(r), max(opts.minN, nCoef + 2)); %#ok<AGROW>
                end
                continue
            end
            for k = find(isCat).', d.(covVar(k)) = removecats(d.(covVar(k))); end
            mdl = fitlm(d, char(formula));
            if mdl.NumEstimatedCoefficients < mdl.NumCoefficients
                roiNotes(end+1,1) = sprintf('%s / %s: rank-deficient design (no model).', st, roiNames(r)); %#ok<AGROW>
                continue
            end
            dc = d;
            for k = find(~isCat).', dc.(covVar(k)) = zeros(height(dc),1); end
            adj(r) = mean(predict(mdl, dc));
            ct = mdl.Coefficients;
            names = string(ct.Properties.RowNames);
            for m = 1:numel(names)
                if names(m) == "(Intercept)", continue; end
                key = matlab.lang.makeValidName(names(m));
                if ~isfield(coef, key), coef.(key) = NaN(nR,4); end
                coef.(key)(r,:) = [ct.Estimate(m), ct.SE(m), ct.tStat(m), ct.pValue(m)];
            end
        end
        out.(st + "_n") = n;  out.(st + "_mean") = mu;  out.(st + "_sd") = sd;  out.(st + "_adjMean") = adj;
        terms = string(fieldnames(coef));
        termCols.(st) = terms;
        for m = 1:numel(terms)
            c = coef.(terms(m));
            out.(st + "_" + terms(m) + "_beta") = c(:,1);
            out.(st + "_" + terms(m) + "_se")   = c(:,2);
            out.(st + "_" + terms(m) + "_t")    = c(:,3);
            out.(st + "_" + terms(m) + "_p")    = c(:,4);
            out.(st + "_" + terms(m) + "_q")    = bhFDR(c(:,4));
        end
    end

    %% ---- write ----
    if ~isfolder(outDir), mkdir(outDir); end
    writetable(out, outFiles(1), 'WriteRowNames', true);
    Lw = L;
    for k = find(~isCat).', Lw.(covVar(k)) = Lw.(covVar(k)) + centers(k); end   % write raw values
    writetable(Lw, outFiles(2));
    ST = table(subj.id, subj.included, subj.reason, subj.runDir, 'VariableNames', {'subject','included','reason','runDir'});
    ST = [ST, S(:, 2:end)];
    writetable(ST, outFiles(3));

    R = strings(0,1);
    R(end+1,1) = "ASE group ROI statistics";
    R(end+1,1) = "created:        " + string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    R(end+1,1) = "allSubjDir:     " + allSubjDir;
    R(end+1,1) = "covariate table: " + string(covTableFile);
    R(end+1,1) = "fits:           <subj>/" + opts.fitsSubdir + ", tag: " + tern(strlength(string(opts.fitTag)) > 0, string(opts.fitTag), "(none: latest run)");
    R(end+1,1) = "ROI table:      " + opts.tableName;
    R(end+1,1) = "stats:          " + strjoin(stats, ', ');
    R(end+1,1) = "model:          " + formula;
    R(end+1,1) = "minVox: " + opts.minVox + ", minN: " + opts.minN + ", settings check: " + opts.settingsCheck;
    R(end+1,1) = "";
    R(end+1,1) = "Selection (" + height(C) + " uncommented rows):";
    if isempty(fdesc), R(end+1,1) = "  (no filters: all rows)"; else, R = [R; "  " + fdesc]; end
    R(end+1,1) = "  selected: " + nSel + ", included: " + nnz(subj.included);
    if ~isempty(cinfo.commented)
        R(end+1,1) = "  commented-out rows (skipped): " + strjoin(cinfo.commented, ', ');
    end
    R(end+1,1) = "";
    R(end+1,1) = "Covariates:";
    if isempty(covVar), R(end+1,1) = "  (none)"; end
    for k = 1:numel(covVar)
        if isCat(k)
            lv = categories(L.(covVar(k)));
            R(end+1,1) = "  " + covVar(k) + ": categorical, levels " + strjoin(string(lv), ', ') + ...
                " (reference " + string(lv{1}) + ")"; %#ok<AGROW>
        else
            x = Ccov.(covVar(k));
            R(end+1,1) = sprintf('  %s: numeric, centred on %.4g (range %.4g to %.4g)', covVar(k), centers(k), min(x), max(x)); %#ok<AGROW>
        end
    end
    R(end+1,1) = "";
    R(end+1,1) = "Subjects:";
    for i = 1:nSel
        if subj.included(i)
            R(end+1,1) = "  + " + subj.id(i) + "  " + subj.runDir(i); %#ok<AGROW>
        else
            R(end+1,1) = "  - " + subj.id(i) + "  EXCLUDED: " + subj.reason(i); %#ok<AGROW>
        end
    end
    if ~isempty(unlisted)
        R(end+1,1) = "";
        R(end+1,1) = "Warning: subject folder(s) with a matching fit run but no row in the covariate table:";
        R(end+1,1) = "  " + strjoin(unlisted, ', ');
    end
    if ~isempty(commentedWithFit)
        R(end+1,1) = "";
        R(end+1,1) = "Note: commented-out subject(s) that do have a matching fit run: " + strjoin(commentedWithFit, ', ');
    end
    R(end+1,1) = "";
    R(end+1,1) = "Settings check:";
    if opts.settingsCheck == "off", R(end+1,1) = "  off";
    elseif isempty(settingsMsg), R(end+1,1) = "  all included runs have the same fit settings";
    else, R = [R; "  WARNING: " + settingsMsg];
    end
    if ~isempty(report), R(end+1,1) = ""; R = [R; report]; end
    R(end+1,1) = "";
    for s = 1:numel(stats)
        st = stats(s);
        nm = out.(st + "_n");
        R(end+1,1) = sprintf('%s: %d ROIs, %d with a model; subjects per ROI %d to %d', st, nR, ...
            nnz(isfinite(out.(st + "_adjMean"))), min(nm), max(nm)); %#ok<AGROW>
        for term = termCols.(st).'
            q = out.(st + "_" + term + "_q");
            R(end+1,1) = sprintf('  %s: %d ROI(s) with q < 0.05', term, nnz(q < 0.05)); %#ok<AGROW>
        end
    end
    if ~isempty(roiNotes)
        R(end+1,1) = "";
        R(end+1,1) = "ROIs without a model:";
        R = [R; "  " + roiNotes];
    end
    writelines(R, outFiles(4));
    if v, fprintf('%s\n', strjoin(R, newline)); fprintf('Wrote %s\n', strjoin(outFiles, ', ')); end

    res = struct('table', out, 'long', Lw, 'subjects', ST, 'unlisted', unlisted, ...
        'settingsMsg', settingsMsg, 'centers', centers, 'covariates', covVar, 'files', outFiles);
end

function s = toStrings(x)
    if isempty(x), s = strings(0,1); return; end
    s = string(x);  s = s(:);
    if isscalar(s) && contains(s, ","), s = strtrim(split(s, ",")); end
    s = s(strlength(s) > 0);
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end
