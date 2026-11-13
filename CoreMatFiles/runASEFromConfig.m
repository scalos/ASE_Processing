function runASEFromConfig(B0, hdr_File, sorted_ASE_Files, runOutDir, opts)
% Reads the "fit" section of a JSON config, applies command-line overrides,
% and calls ASE_fitting. Per-subject inputs are passed directly.
%
%   B0, hdr_File, sorted_ASE_Files : passed straight to ASE_fitting
%   runOutDir     : run folder; config_resolved.json is written here
%   opts.cfgFile  : path to config json (optional; empty = MATLAB defaults)
%   opts.overrides: JSON object string of fit overrides, e.g. '{"Hct":0.4}'
%   opts.segFiles, opts.segLUTs : cell arrays of paths (may be empty)
%   opts.qualityMaskFiles : cell array of quality masks on the ASE grid (may be empty)
%   opts.outPref  : full output prefix (path included)
%
% Config values: "Inf", "-Inf", "NaN" (strings) become numbers; null becomes [].
% Overrides of nested blocks (e.g. --fit 'quality={"enabled":true}') are merged into the
% config's block key by key rather than replacing it.
    arguments
        B0 (1,1) double
        hdr_File
        sorted_ASE_Files (:,:) cell
        runOutDir
        opts.cfgFile = ''
        opts.overrides = ''
        opts.segFiles = {}
        opts.segLUTs = {}
        opts.qualityMaskFiles = {}
        opts.outPref = ''
    end
    assert(~isempty(opts.outPref), 'outPref must be provided.');
    reserved = {'segFiles','segLUTs','qualityMaskFiles','outPref'};

    % ---- read config ----
    fit = struct();
    if ~isempty(opts.cfgFile)
        assert(isfile(opts.cfgFile), 'Config file not found: %s', opts.cfgFile);
        cfg = jsondecode(fileread(opts.cfgFile));
        if isfield(cfg,'fit') && isstruct(cfg.fit)
            fit = cfg.fit;
        end
    end

    % ---- apply overrides ----
    if ~isempty(strtrim(opts.overrides))
        try
            ov = jsondecode(opts.overrides);
        catch ME
            error('Could not parse --fit overrides as JSON (%s): %s', opts.overrides, ME.message);
        end
        for f = fieldnames(ov)'
            if isfield(fit, f{1}) && isstruct(fit.(f{1})) && isstruct(ov.(f{1}))
                % nested blocks (e.g. "quality") are merged key by key, so an override
                % changes only the keys it names
                for g = fieldnames(ov.(f{1}))'
                    fit.(f{1}).(g{1}) = ov.(f{1}).(g{1});
                end
            else
                fit.(f{1}) = ov.(f{1});
            end
        end
    end

    fit = convertSpecial(fit);
    if isfield(fit, 'quality')
        fit.quality = qualityDefaults(fit.quality);   % validate keys; record the full block
    end

    bad = intersect(fieldnames(fit), reserved);
    if ~isempty(bad)
        error('Config/override keys not allowed in "fit" (set per subject by the batch script): %s', strjoin(bad, ', '));
    end

    % ---- build name-value list ----
    nv = {'segFiles', opts.segFiles, 'segLUTs', opts.segLUTs, 'qualityMaskFiles', opts.qualityMaskFiles, ...
          'outPref', opts.outPref};
    for f = fieldnames(fit)'
        nv(end+1:end+2) = {f{1}, fit.(f{1})}; %#ok<AGROW>
    end

    % ---- run; write the resolved config on success, or on failure if the folder has content ----
    try
        ASE_fitting(B0, hdr_File, sorted_ASE_Files, nv{:});
    catch ME
        d = dir(runOutDir);
        if any(~ismember({d.name}, {'.','..'}))
            writeResolved(runOutDir, B0, hdr_File, sorted_ASE_Files, opts, fit);
        end
        rethrow(ME);
    end
    writeResolved(runOutDir, B0, hdr_File, sorted_ASE_Files, opts, fit);
end

% ---------------------------------------------------------------------------
function writeResolved(runOutDir, B0, hdr_File, aseFiles, opts, fit)
    s = struct('created', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
               'cfgFile', char(opts.cfgFile), ...
               'overrides', char(opts.overrides), ...
               'B0', B0, ...
               'hdr_File', char(hdr_File), ...
               'aseFiles', {num2cell(aseFiles,2)}, ...
               'segFiles', {opts.segFiles}, ...
               'segLUTs', {opts.segLUTs}, ...
               'qualityMaskFiles', {opts.qualityMaskFiles}, ...
               'outPref', char(opts.outPref), ...
               'fit', encodeSpecial(fit));
    fid = fopen(fullfile(runOutDir,'config_resolved.json'), 'w');
    if fid == -1
        warning('Could not write config_resolved.json in %s', runOutDir);
        return
    end
    fprintf(fid, '%s', jsonencode(s, 'PrettyPrint', true));
    fclose(fid);
end

function v = convertSpecial(v)
% JSON-decoded value -> value ASE_fitting expects.
    if isstruct(v)
        for f = fieldnames(v)'
            v.(f{1}) = convertSpecial(v.(f{1}));
        end
    elseif iscell(v)
        v = cellfun(@convertSpecial, v, 'UniformOutput', false);
        isNum = cellfun(@(c) (isnumeric(c) || islogical(c)) && isscalar(c), v);
        if ~isempty(v) && all(isNum(:))
            v = [v{:}];                      % e.g. {0;Inf} -> [0 Inf]
        else
            v = v(:).';                      % cell row (e.g. string lists)
        end
    elseif ischar(v) || isstring(v)
        switch lower(strtrim(char(v)))
            case {'inf','+inf','infinity'}, v = Inf;
            case {'-inf','-infinity'},      v = -Inf;
            case 'nan',                     v = NaN;
        end
    elseif (isnumeric(v) || islogical(v)) && isvector(v)
        v = v(:).';                          % row vectors
    end
end

function v = encodeSpecial(v)
% Inverse of convertSpecial for writing valid JSON.
    if isstruct(v)
        for f = fieldnames(v)'
            v.(f{1}) = encodeSpecial(v.(f{1}));
        end
    elseif iscell(v)
        v = cellfun(@encodeSpecial, v, 'UniformOutput', false);
    elseif isnumeric(v) && any(~isfinite(v(:)))
        c = num2cell(v);
        for i = 1:numel(c)
            if isnan(c{i}),        c{i} = 'NaN';
            elseif c{i} == Inf,    c{i} = 'Inf';
            elseif c{i} == -Inf,   c{i} = '-Inf';
            end
        end
        if isscalar(c), v = c{1}; else, v = c; end
    end
end