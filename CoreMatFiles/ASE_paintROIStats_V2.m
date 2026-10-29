function ASE_paintROIStats_V2(templateFiles, segLUTs, statsFile, outPref, opts)
    arguments
        templateFiles
        segLUTs       %path to combined GM/WM/BS lookup table (code<TAB>name)
        statsFile    %path to per-subject ROI stats csv (tab-delimited, row names = ROI names)
        outPref      %output prefix; files written as <outPref>_template_<STAT>.nii.gz
        opts.stats = {'OEF','R2P','vCBV'}
        opts.verbose = true
        opts.overwrite = true
    end

    v = opts.verbose;

    % combine LUT and template files:
    [segTemp,LUT] = combineSegVols(templateFiles,segLUTs);
    tmpHdr = niftiinfo(segTemp);
    
    % --- load ROI stats table ---
    statsTbl = readtable(statsFile,'FileType','text','Delimiter',',','ReadRowNames',true);
    roiNames = statsTbl.Properties.RowNames;

    for s = 1:numel(opts.stats)
        statName = opts.stats{s};
        if ~ismember(statName, statsTbl.Properties.VariableNames)
            warning('Stat "%s" not found in %s -- skipping.', statName, statsFile);
            continue;
        end

        outVol = nan(volSize);

        for r = 1:numel(roiNames)
            name = roiNames{r};
            val  = statsTbl.(statName)(r);

            code = LUT.code(strcmp(LUT.name, name));
            if isempty(code)
                warning('ROI "%s" not found in LUT %s -- skipping.', name, segLUTs);
                continue;
            end
            code = code(1);
            
            if any(segTemp(:)==code)
                mask = (segTemp == code);
            else
                warning('Code %d (%s) not found in any template -- skipping.', code, name);
                continue;
            end

            outVol(mask) = val;

            if v
                fprintf('  %s (code %d): painted %d voxels with %s = %g\n', name, code, nnz(mask), statName, val);
            end
        end

        outFileFinal = sprintf('%s_template_%s.nii', outPref, statName);
        outFileFinal = findGoodPath(outFileFinal, opts.overwrite);

        hdr = tmpHdr;
        hdr.Datatype = 'single';
        hdr.BitsPerPixel = 32;
        hdr.ImageSize = volSize;
        hdr.PixelDimensions = hdr.PixelDimensions(1:3);
        
        niftiwrite(single(outVol), writeName, hdr, 'Compressed', true);

        if v
            fprintf('Wrote %s stat volume to %s\n', statName, outFileFinal);
        end
    end
end
