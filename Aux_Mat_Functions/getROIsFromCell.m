function [ROIs,ROI_Names] = getROIsFromCell(ROI_Files,opts)
    arguments
        ROI_Files
        opts.verbose = true;
    end
    v = opts.verbose;
    if v; fprintf('Deriving ROI''s from cell array\n');end
    
    ROIs = cell(numel(ROI_Files),1);
    ROI_Names = cell(numel(ROI_Files),1);
    for idx = 1:numel(ROI_Files)
        [~,name,~] = niftiFileparts(ROI_Files{idx});
        roiData = loadNifti(ROI_Files{idx},'niiName',sprintf('ROI %s',name),'warn');
        if ~isempty(roiData)
            ROIs{idx} = roiData;
            ROI_Names{idx} = name;
        end
    end
    if all(cellfun(@(x) isempty(x),ROIs))
        %if all failed, set to empty singular
        ROIs = {};
        ROI_Names = {};
    end
end