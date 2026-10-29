function [ROIs,ROI_Names] = getROIsFromLUT(segFiles,segLUTs,opts)
    arguments
        segFiles
        segLUTs
        opts.verbose = true;
    end
    v = opts.verbose;
    if v; fprintf('Deriving ROI''s from segmented volume(s)\n');end
    [segVol,LUT] = combineSegVols(segFiles,segLUTs);
    if isempty(segVol)||isempty(LUT)
        error('Unable to derive ROI''s');
    end
    nRegions = height(LUT);
    ROIs = cell(nRegions,1);
    ROI_Names = cell(nRegions,1);
    for idx = 1:nRegions
        code = LUT.code(idx);
        name = LUT.name{idx};
        roi = segVol;
        roi(roi~=code) = 0;
        roi(roi>0) = 1;
        if sum(roi==1,'all')==0
            warning('ROI %s is empty!',name);
        end
        ROIs{idx} = roi;
        ROI_Names{idx} = name;
    end
end