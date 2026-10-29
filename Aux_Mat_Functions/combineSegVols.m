function [segVol,LUT_found] = combineSegVols(segFiles,segLUTs)
    %Combines segmented volumes using codes found in LUTs.
    arguments
        segFiles
        segLUTs
    end
    if ~iscell(segFiles); segFiles = {segFiles}; end
    if ~iscell(segLUTs); segLUTs = {segLUTs}; end
    
    %check that segFiles have matching geometry
    isequalNiftiGeometry(segFiles, 'onMismatch', 'error');
    

    LUT = readtable(segLUTs{1}, 'FileType','text', 'Delimiter','\t', ...
             'ReadVariableNames', false);
    LUT.Properties.VariableNames = {'code','name'};
    for idx = 2:numel(segLUTs)
        thisLUT = readtable(segLUTs{idx}, 'FileType','text', 'Delimiter','\t', ...
                              'ReadVariableNames', false);
        thisLUT.Properties.VariableNames = {'code','name'};
        LUT = [LUT; thisLUT];
    end
    nSegs = numel(segFiles);
    segVols = cell(nSegs,1);
    for idx = 1:nSegs
        segVols{idx} = niftiread(segFiles{idx});
    end
    

    segVol = zeros(size(segVols{1}))-1;
    LUT_found = cell(height(LUT),2);
    for idx = 1:height(LUT)
        code = LUT{idx,1};
        name = LUT{idx,2};
        name = name{1};
        foundSeg = zeros(nSegs,1);
        for sIdx = 1:nSegs
            seg = segVols{sIdx};
            if any(seg(:) == code)
                foundSeg(sIdx) = 1;
            end
        end

        matchingVolInds = [];
        if sum(foundSeg) == 0
            warning('Could not find code %d (%s) in segmented volumes. Skipping this region...',code,name);
        elseif sum(foundSeg) == 1
            matchingVolInds = segVols{foundSeg==1} == code;
        elseif sum(foundSeg) > 1
            matchingSegs = segVols(foundSeg==1);
            firstVolInds = matchingSegs{1} == code;
            segsMatch = true;
            for mIdx = 2:numel(matchingSegs)
                segsMatch = segsMatch && isequal(matchingSegs{mIdx} == code, firstVolInds);
            end
            if segsMatch
                matchingVolInds = firstVolInds;
            else
                warning('Code %d (%s) appears in multiple segmented volumes and covers different regions. Skipping this region...',code,name);
            end
        end
        if ~isempty(matchingVolInds)
            segVol(matchingVolInds) = code;
            LUT_found(idx,:) = {code,name};
        end
    end

    if all(cellfun(@(x) isempty(x),LUT_found))
        %if all failed, set to empty singular
        segVol = [];
        LUT_found = table;
    else
        segVol = int16(segVol);
        %get only non-empty rows:
        LUT_found = LUT_found(any(~cellfun('isempty', LUT_found), 2), :);
        LUT_found = cell2table(LUT_found,"VariableNames",{'code','name'});
    end
end