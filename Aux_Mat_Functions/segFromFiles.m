function [segVol,LUT] = segFromFiles(regionFiles,opts)
    arguments 
        regionFiles %cell list of paths to nifti roi's
        opts.LUT_tmp = [] %optional template (code | name)
    end
    %check that regionFiles can be compared:
    isequalNiftiGeometry(regionFiles,"onMismatch","error");
    LUT_cell = cell(numel(regionFiles),2);
    LUT_tmp = opts.LUT_tmp;
    is1to1 = true;
    if ~isempty(LUT_tmp)
        if height(LUT_tmp)~=numel(regionFiles)
            is1to1 = false;
            warning(['Provided LUT cannot be 1:1 with %d entries for %d regions. ' ...
                'Attempting to use as partial LUT'],height(LUT_tmp),numel(regionFiles));
        end
    end
    
    for idx = 1:numel(regionFiles)
        %If is1to1, treat LUT_tmp as a direct mapping:
        if is1to1
            LUT_cell(idx,:) = LUT_tmp{idx,:};
        end

    end
        
end