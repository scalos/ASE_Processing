function rep = ASE_paintROIStats(templateFile, segLUT, statsFile, outPath, opts, segOpts)
% Paints ROI values from a table onto template segmentation(s): every voxel of a region gets
% that region's value. Works for per-subject tables (ase_fit_results_ROI_Stats.csv) and group
% tables (group_ROI_Stats.csv; e.g. stats "OEF_adjMean","OEF_age_t").
%
% Regions are placed with combineSegVols(..., keepOverlaps=false): in LUT order, and a region
% that overlaps one already placed is dropped whole (no trimming). To paint a different set of
% regions, pass a different (painting) LUT. Voxels with no painted region get opts.defaultVal
% (NaN by default, so "no value" is never confused with 0).
%
% One report is written next to the outputs (outPath with %STAT% -> "paint_report", .txt), listing:
% regions painted; template regions dropped (overlap, not in template, ...); table rows whose
% name is not in the LUT; LUT regions with no table row; and, per stat, regions whose value is NaN.
    arguments
        templateFile        %template segmentation file(s) (cell or char)
        segLUT              %LUT file(s) (code<TAB>name) for the template segmentation(s)
        statsFile           %ROI table (csv, row names = ROI names)
        outPath             %output path without extension. %STAT% is replaced by the stat name.
        opts.compressOutputs = true;
        opts.stats = {'OEF','R2P','vCBV'}
        opts.verbose = true
        opts.overwrite = true
        opts.defaultVal = NaN;      %value of voxels outside every painted region
        segOpts.omitCodes
        segOpts.omitNames
        segOpts.warnUnusedCodes
        segOpts.bkdCode
    end
    v = opts.verbose;
    if ~iscell(templateFile), templateFile = {templateFile}; end
    if ~iscell(segLUT), segLUT = {segLUT}; end
    stats = string(opts.stats);  stats = stats(:).';
    segArgs = namedargs2cell(segOpts);
    [segTemp, LUT, segRep] = combineSegVols(templateFile, segLUT, ...
        "verbose", false, "keepOverlaps", false, "warnSkipped", false, segArgs{:});
    assert(~isempty(segTemp), 'No LUT region could be placed on the template(s).');
    tmpHdr = niftiinfo(templateFile{1});
    volSize = size(segTemp);

    statsTbl = readtable(statsFile, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
    roiNames = string(statsTbl.Properties.RowNames);
    lutNames = string(LUT.name);  lutCodes = double(LUT.code);

    % match table rows to placed regions (case-insensitive, like before)
    [inLUT, loc] = ismember(lower(roiNames), lower(lutNames));
    allLutNames = string(segRep.name);
    notInLUT   = roiNames(~ismember(lower(roiNames), lower(allLutNames)));
    droppedRow = roiNames(~inLUT & ismember(lower(roiNames), lower(allLutNames)));
    noValue    = lutNames(~ismember(lower(lutNames), lower(roiNames)));

    % voxel -> placed-region index, computed once
    [~, regIdx] = ismember(double(segTemp), lutCodes);       % 0 = background / unplaced
    nVoxRegion = accumarray(regIdx(regIdx > 0), 1, [numel(lutCodes) 1]);

    [outDir, outBase] = fileparts(string(outPath));
    rep = struct('painted', strings(0,1), 'files', strings(0,1), 'segReport', segRep, ...
        'notInLUT', notInLUT, 'noTableRow', noValue, 'nanValues', struct());
    R = strings(0,1);
    R(end+1,1) = "ROI painting report";
    R(end+1,1) = "created:   " + string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    R(end+1,1) = "table:     " + string(statsFile);
    R(end+1,1) = "templates: " + strjoin(string(templateFile), ', ');
    R(end+1,1) = "LUTs:      " + strjoin(string(segLUT), ', ');
    R(end+1,1) = "fill value outside painted regions: " + string(opts.defaultVal);
    R(end+1,1) = "";
    R(end+1,1) = sprintf('Template regions: %d in LUT(s), %d placed on the template.', height(segRep), numel(lutCodes));
    skipped = segRep(segRep.status ~= "placed", :);
    for st = ["overlap","volumes disagree","not in volumes","omitted"]
        rows = skipped(skipped.status == st, :);
        if isempty(rows), continue; end
        R(end+1,1) = sprintf('  Dropped (%s): %d', st, height(rows)); %#ok<AGROW>
        for k = 1:height(rows)
            line = sprintf('    %d %s', rows.code(k), rows.name(k));
            if strlength(rows.detail(k)) > 0, line = line + "  [" + rows.detail(k) + "]"; end
            R(end+1,1) = line; %#ok<AGROW>
        end
    end
    R(end+1,1) = sprintf('Table rows: %d; painted: %d.', numel(roiNames), nnz(inLUT));
    if ~isempty(notInLUT)
        R(end+1,1) = "  Not in the LUT(s) (not painted): " + strjoin(notInLUT, ', ');
    end
    if ~isempty(droppedRow)
        R(end+1,1) = "  In the LUT but dropped from the template (not painted): " + strjoin(droppedRow, ', ');
    end
    if ~isempty(noValue)
        R(end+1,1) = "  Placed regions with no table row (left at the fill value): " + strjoin(noValue, ', ');
    end

    colNames = string(statsTbl.Properties.VariableNames);
    for s = 1:numel(stats)
        statName = stats(s);
        if ~ismember(statName, colNames)
            R(end+1,1) = "Stat " + statName + ": not a column of the table (skipped)."; %#ok<AGROW>
            warning('Stat "%s" not found in %s; skipped.', statName, statsFile);
            continue
        end
        vals = NaN(numel(lutCodes),1);
        vals(loc(inLUT)) = statsTbl.(statName)(inLUT);
        lut = [opts.defaultVal; vals];                          % index 1 = outside
        outVol = reshape(lut(regIdx + 1), volSize);
        nanRegions = lutNames(isnan(vals) & ~ismember(lutNames, noValue));
        rep.nanValues.(matlab.lang.makeValidName(statName)) = nanRegions;
        line = sprintf('Stat %s: %d regions painted, %d voxels', statName, nnz(~isnan(vals)), sum(nVoxRegion(~isnan(vals))));
        if ~isempty(nanRegions), line = line + "; NaN value in: " + strjoin(nanRegions, ', '); end
        R(end+1,1) = line; %#ok<AGROW>

        outName = replace(outBase, "%STAT%", statName);
        outFileFinal = fullfile(outDir, outName + ".nii");
        if opts.compressOutputs, outFileFinal = outFileFinal + ".gz"; end
        outFileFinal = findGoodPath(outFileFinal, opts.overwrite);
        hdr = tmpHdr;
        hdr.Datatype = 'single';  hdr.BitsPerPixel = 32;
        hdr.ImageSize = volSize;  hdr.PixelDimensions = hdr.PixelDimensions(1:3);
        writePath = regexprep(outFileFinal, "\.gz$", "");     % niftiwrite adds .gz when compressing
        niftiwrite(single(outVol), writePath, hdr, 'Compressed', opts.compressOutputs);
        rep.painted(end+1,1) = statName;  rep.files(end+1,1) = outFileFinal;
        if v, fprintf('Wrote %s -> %s\n', statName, outFileFinal); end
    end
    if contains(outBase, "%STAT%")
        repFile = fullfile(outDir, replace(outBase, "%STAT%", "paint_report") + ".txt");
    else
        repFile = fullfile(outDir, outBase + "_paint_report.txt");
    end
    writelines(R, repFile);
    rep.reportFile = repFile;
    if v
        fprintf('%s\n', strjoin(R, newline));
        fprintf('Report: %s\n', repFile);
    end
end
