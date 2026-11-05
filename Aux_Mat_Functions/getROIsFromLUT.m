function [ROIs,ROI_Names] = getROIsFromLUT(segFiles,segLUTs,opts)
    arguments
        segFiles
        segLUTs
        opts.verbose = true;
    end
    v = opts.verbose;
    if v; fprintf('Deriving ROI''s from segmented volume(s)\n'); end

    [segVol,LUT] = combineSegVols(segFiles,segLUTs,"keepOverlaps",true,"verbose",v);
    if isempty(segVol) || isempty(LUT)
        error('Unable to derive ROI''s');
    end

    nRegions  = height(LUT);
    nLayers   = size(segVol,4);          % assumes 3D seg volumes
    volSize   = size(segVol,1:3);
    lutCodes  = double(LUT.code);
    ROIs      = cell(nRegions,1);
    ROI_Names = cell(nRegions,1);
    found     = false(nRegions,1);

    for lIdx = 1:nLayers
        layerVol = segVol(:,:,:,lIdx);               % sliced once per layer
        [sortedVals, order] = sort(layerVol(:));     % one sort per layer
        clear layerVol

        % boundaries of each run of identical labels in the sorted list
        starts = find([true; sortedVals(2:end) ~= sortedVals(1:end-1)]);
        ends   = [starts(2:end)-1; numel(sortedVals)];
        uVals  = double(sortedVals(starts));
        clear sortedVals

        % which LUT regions live in this layer, and where in uVals
        [inLayer, loc] = ismember(lutCodes, uVals);

        for idx = find(inLayer & ~found)'
            k = loc(idx);
            roi = false(volSize);
            roi(order(starts(k):ends(k))) = true;
            ROIs{idx}      = roi;
            ROI_Names{idx} = char(LUT.name(idx));
            found(idx)     = true;
        end
    end

    for idx = find(~found)'
        warning('ROI %s is empty!', char(LUT.name(idx)));
    end
end