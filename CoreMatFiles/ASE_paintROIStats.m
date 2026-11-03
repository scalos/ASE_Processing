function ASE_paintROIStats(templateFiles, segLUTs, statsFile, outPath, opts, segOpts)
    arguments
        templateFiles
        segLUTs       %path to combined GM/WM/BS lookup table (code<TAB>name)
        statsFile    %path to per-subject ROI stats csv (tab-delimited, row names = ROI names)
        outPath      %output path. Use glob %STAT% in name to include stat name.
        opts.compressOutputs = true;
        opts.stats = {'OEF','R2P','vCBV'}
        opts.verbose = true
        opts.overwrite = true
        opts.defaultVal = 0; %initialization value used for painted volume
        opts.comb_verbose = false;
        segOpts.omitCodes
        segOpts.omitNames
        segOpts.warnUnusedCodes
        segOpts.bkdCode
        
    end

    v = opts.verbose;

    % combine LUT and template files:
    if v; fprintf('Combining template and LUT files...\n');end
    fwd = namedargs2cell(segOpts);
    [segTemp,LUT] = combineSegVols(templateFiles,segLUTs,'verbose',opts.comb_verbose,fwd{:});

    tmpHdr = niftiinfo(templateFiles{1});
    % --- load ROI stats table ---
    statsTbl = readtable(statsFile,'ReadRowNames',true);
    roiNames = statsTbl.Properties.RowNames;
    maxNameLen = max(strlength(LUT.name),[],'all');
    maxCodeLen = max(strlength(string(num2str(LUT.code(:)))),[],'all');
    maxStatLen = max(strlength(opts.stats),[],'all');
    maxVoxLen = strlength(string(numel(segTemp)));

    [~,statsFileName] = fileparts(statsFile); 

    volSize = size(segTemp);
    lutNames = {};
    for idx = 1:numel(segLUTs)
        [~,fileName] = fileparts(segLUTs{idx}); 
        lutNames(end+1) = {char(fileName)}; %#ok
    end
    lutNames = strjoin(lutNames,', ');
    if v; fprintf('Painting Stats...\n'); end

    badNames = {};

    for s = 1:numel(opts.stats)
        statName = opts.stats{s};
        colNames = statsTbl.Properties.VariableNames;
        if ~ismember(statName, colNames)
            warning(['Stat "%s" not found in %s with columns: %s\n' ...
                     '=> Skipping this stat.'], statName, statsFileName,strjoin(colNames,', '));
            continue;
        end

        outVol = ones(volSize)*opts.defaultVal;

        for r = 1:numel(roiNames)
            name = roiNames{r};
            if ismember(name,badNames)
                if v; fprintf('Skipping region "%s" (previously determined invalid)',name); end
                continue;
            end
            val  = statsTbl.(statName)(r);

            code = LUT.code(strcmpi(LUT.name, name));
            if isempty(code)
                warning(['Region with name "%s" not found in LUTs ("%s")\n' ...
                         '=> Skipping this region.'], name, lutNames);
                badNames{end+1} = name; %#ok
                continue;
            end
            code = code(1);
            
            if any(segTemp(:)==code)
                mask = (segTemp == code);
                outVol(mask) = val;
            else
                warning(['Code %d (%s) not found in any template\n' ...
                         '=> Skipping this region.'], code, name);
                badNames{end+1} = name; %#ok
                continue;
            end

            if v
                fprintf('%*s (code %*d): painted %*d voxels with %*s = %g\n',...
                    maxNameLen,name,maxCodeLen,code,maxVoxLen,nnz(mask),maxStatLen,statName,val);
            end
        end
        
        [outPath_,outName] = fileparts(outPath);
        outName = replace(outName,"%STAT%",statName);
        outFileFinal = sprintf("%s.nii.gz",fullfile(outPath_,outName));
        outFileFinal = findGoodPath(outFileFinal, opts.overwrite);

        hdr = tmpHdr;
        hdr.Datatype = 'single';
        hdr.BitsPerPixel = 32;
        hdr.ImageSize = volSize;
        hdr.PixelDimensions = hdr.PixelDimensions(1:3);
        
        niftiwrite(single(outVol), outFileFinal, hdr, 'Compressed', true);

        if v
            fprintf('Wrote %s stat volume to %s\n', statName, outFileFinal);
        end
    end
end
