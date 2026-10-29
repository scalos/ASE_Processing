function outPath = writeMap(templateHdr,mapData,outFile,opts)
    arguments
        templateHdr %nifti header for map
        mapData     %map data
        outFile     %out file name (includes full path)
        opts.verbose = true
        opts.overwrite = false
        opts.mapName = ''
    end

    outPath = findGoodPath(outFile,opts.overwrite);
    
    templateHdr.Datatype = 'single';
    templateHdr.BitsPerPixel = 32;
    templateHdr.ImageSize = size(mapData);
    templateHdr.PixelDimensions = templateHdr.PixelDimensions(1:3);
    niftiwrite(single(mapData),outPath,templateHdr);
    if opts.verbose
        mapName = opts.mapName;
        if isempty(mapName)
            mapName = 'a'; 
        end
        fprintf('Wrote %s map to %s\n',mapName,outPath);
    end
end