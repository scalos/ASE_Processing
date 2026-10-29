function outPath = writeVec(templateHdr,vecData,outFile,opts)
    arguments
        templateHdr %nifti header for map
        vecData (:,:,:,:)     %vector data. assumed that 4th dimension is vectors
        outFile     %out file name (includes full path)
        opts.verbose = true
        opts.overwrite = false
        opts.vecName = '';
        opts.vecDesc = '';
        opts.intent = 'FitParams'
    end

    outPath = findGoodPath(outFile,opts.overwrite);
    vecLen = size(vecData,4);
    templateHdr.Datatype = 'single';
    templateHdr.BitsPerPixel = 32;
    templateHdr.ImageSize = size(vecData);
    templateHdr.PixelDimensions = [templateHdr.PixelDimensions(1:3),vecLen];
    templateHdr.raw.intent_code = 1007;
    templateHdr.raw.intent_name = opts.intent;
    if ~isempty(opts.vecName)
        templateHdr.Description = opts.vecDesc;
    end
    niftiwrite(single(vecData),outPath,templateHdr);
    if opts.verbose
        vecName = opts.vecName;
        if isempty(vecName)
            vecName = 'a'; 
        end
        fprintf('Wrote %s vec file to %s\n',vecName,outPath);
    end
end