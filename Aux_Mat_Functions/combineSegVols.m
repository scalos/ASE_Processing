function [segVol,LUT_found] = combineSegVols(segFiles,segLUTs,opts)
    %Combines segmented volumes using codes found in LUTs. 
    % Duplicate entries in LUT's must match completely otherwise mismatch 
    % errors will be triggered. 
    % Duplicate codes across segmented volumes will be allowed provided
    % that they cover exactly the same region in each volume
    arguments
        segFiles
        segLUTs
        opts.bkdCode int32 = 0;
        opts.verbose = false;
        opts.omitCodes = [];
        opts.omitNames cell = {};
        opts.warnUnusedCodes = false; %warn if some codes in volumes are not found in LUT's
    end
    v = opts.verbose;

    function tbl = readLUT(lutFile)
        tbl = readtable(lutFile, 'ReadVariableNames', false, 'TextType', 'string');
        if width(tbl) ~= 2
            error('LUT:badFormat', 'LUT "%s" has %d columns, expected 2 (code,name).', ...
                lutFile, width(tbl));
        end
        tbl.Properties.VariableNames = {'code','name'};
    
        if ~isnumeric(tbl.code)
            % Text in the code column means the first row was a header: drop it
            tbl(1,:) = [];
            tblCodes = str2double(tbl.code);
            if any(isnan(tblCodes))
                error('LUT:badCodes', 'LUT "%s" has non-numeric codes: %s', ...
                    lutFile, strjoin("'" + tbl.code(isnan(tblCodes)) + "'", ', '));
            end
            tbl.code = tblCodes;
        end
        tbl.name = string(tbl.name);
    end


    if ~iscell(segFiles); segFiles = {segFiles}; end
    if ~iscell(segLUTs); segLUTs = {segLUTs}; end
    
    %check that segFiles have matching geometry
    isequalNiftiGeometry(segFiles, 'onMismatch', 'error');
    
    LUT = readLUT(segLUTs{1});
    for idx = 2:numel(segLUTs)
        LUT = [LUT;readLUT(segLUTs{idx})]; %#ok<AGROW>
    end
    
    % Trim exact duplicate rows (same code AND same name), keeping first occurrence
    LUT = unique(LUT, 'stable');
    
    % Whatever duplicates remain are, by definition, mismatches
    msgs = strings(0,1);
    
    % Same code, different names
    [uCodes, ~, codeGrp] = unique(LUT.code);
    for k = 1:numel(uCodes)
        rows = (codeGrp == k);
        if nnz(rows) > 1
            names = LUT.name(rows);
            msgs(end+1,1) = sprintf('Code %g has multiple names: %s', ...
                uCodes(k), strjoin("'" + names + "'", ', ')); %#ok<AGROW>
        end
    end
    
    % Same name, different codes
    [uNames, ~, nameGrp] = unique(LUT.name);
    for k = 1:numel(uNames)
        rows = (nameGrp == k);
        if nnz(rows) > 1
            codes = LUT.code(rows);
            msgs(end+1,1) = sprintf('Name ''%s'' has multiple codes: %s', ...
                uNames(k), strjoin(string(codes'), ', ')); %#ok<AGROW>
        end
    end
    
    % 3) Report everything at once and error out
    if ~isempty(msgs)
        error('LUT:conflictingEntries', ...
            'Found %d conflicting lookup table entries:\n  %s', ...
            numel(msgs), strjoin(msgs, newline + "  "));
    end

    nSegs = numel(segFiles);
    segVols = cell(nSegs,1);
    for idx = 1:nSegs
        segVols{idx} = niftiread(segFiles{idx});
    end
    uniqueCodes = cellfun(@(s) double(unique(s(:))), segVols, 'UniformOutput', false);
    missing = setdiff(unique(vertcat(uniqueCodes{:})), [LUT.code; opts.bkdCode]);
    if ~isempty(missing)&&opts.warnUnusedCodes
        missingStr = compressIntArr(missing,"joinChar",',');
        warning('Volume codes not in LUT(s): [%s]', missingStr);
    end
    

    segVol = ones(size(segVols{1}),'int32').*opts.bkdCode;
    LUT_found = cell(height(LUT),2);
    maxCodeLen = max(strlength(string(num2str(LUT.code(:)))),[],'all');
    for idx = 1:height(LUT)
        code = LUT{idx,1};
        name = LUT{idx,2};
        if ismember(code,opts.omitCodes)
            if v; fprintf('Skipping code %d as requested\n',code);end
            continue;
        end
        if ~isempty(opts.omitNames)
            if ismember(name,cell2mat(opts.omitNames))
                if v; fprintf('Skipping name %s as requested\n',name);end
                continue;
            end
        end
        foundSeg = zeros(nSegs,1);
        for sIdx = 1:nSegs
            segCodes = uniqueCodes{sIdx};
            if ismember(code,segCodes)
                foundSeg(sIdx) = 1;
            end
        end
        fileNames = {};
        for sIdx = 1:nSegs
            if foundSeg(sIdx)
                [~,fileName] = niftiFileparts(segFiles{sIdx}); 
                fileNames(end+1) = {char(fileName)}; %#ok
            end
        end
        fileNames = strjoin(fileNames,', ');

        matchingVolInds = [];
        if sum(foundSeg) == 0
            %code was not found in any volumes
            warning(['Issue with code %d (%s):\n' ...
                     '-  Could not find code in segmented volumes.\n' ...
                     '=> Skipping this region...'],code,name);
        elseif sum(foundSeg) == 1
            %code was found exactly once in all volumes
            matchingVolInds = segVols{foundSeg==1} == code;
        elseif sum(foundSeg) > 1
            %code was found in multiple volumes
            matchingSegs = segVols(foundSeg==1);
            firstVolInds = matchingSegs{1} == code;
            segsMatch = all(cellfun(@(s) isequal(s == code, firstVolInds), matchingSegs(2:end)));
            if segsMatch
                %code covers the same region in all volumes
                matchingVolInds = firstVolInds;
            else
                %volumes do not agree with code
                
                warning(['Issue with code %d (%s) from files "%s":\n' ...
                         '-  Code appears in multiple segmented volumes and covers different regions.\n' ...
                         '=> Skipping this region...'],code,name,fileNames);
            end
        end
        if ~isempty(matchingVolInds)
            prev = segVol(matchingVolInds);
            clash = prev ~= opts.bkdCode & prev ~= code;
            if any(clash)
                clashCodes = unique(prev(clash));
                clashNames = strings(size(clashCodes));
                for k = 1:numel(clashCodes)
                    nm = LUT.name(LUT.code == clashCodes(k));
                    if isempty(nm); nm = "(not in LUT)"; end   % defensive fallback
                    clashNames(k) = nm(1);
                end
                clashTbl = strjoin(compose('   %6d   %s',clashCodes(:), clashNames(:)), newline);
            
                warning(['Issue with code %d (%s) from file(s) "%s":\n' ...
                         '-  Overlaps previously assigned voxels:\n' ...
                         '     code   name\n' ...
                         '%s\n' ...
                         '=> Skipping this region...'], code, name, fileNames, clashTbl);
            else
                if v; fprintf('Wrote code %d (%s) from file(s) "%s"\n',code,name,fileNames);end
                segVol(matchingVolInds) = code;
                LUT_found(idx,:) = {code,name};
            end
        end
    end

    if all(cellfun(@(x) isempty(x),LUT_found))
        %if all failed, set to empty singular
        segVol = [];
        LUT_found = table;
    else
        maxCode = max(segVol(:));
        if maxCode > intmax('int32')
            error('Codes exceed int32 range.');
        elseif maxCode > intmax('int16')
            warning('Some codes exceed int16 range, using int32...');
            segVol = int32(segVol);
        else
            segVol = int16(segVol);
        end
        %get only non-empty rows:
        LUT_found = LUT_found(any(~cellfun('isempty', LUT_found), 2), :);
        LUT_found = cell2table(LUT_found,"VariableNames",{'code','name'});
    end
end