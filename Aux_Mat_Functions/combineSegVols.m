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
        opts.keepOverlaps = true; % stack overlapping regions along an extra (last) dimension
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

    volSize = size(segVols{1});
    nVox    = numel(segVols{1});
    if opts.keepOverlaps
        nLayers = nSegs+1; %extra layer to safeguard against greedy assignment     
    else
        nLayers = 1;
    end
    segVol = zeros(nVox, nLayers,'int32') + int32(opts.bkdCode);   % one column per layer
    LUT_found = cell(height(LUT),2);

    maxNameLen = max(strlength(string(LUT.name)));
    maxCodeLen = max(strlength(string(LUT.code)));
    maxLayerDigits = strlength(string(nLayers));
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
            mask = matchingVolInds(:);
            placed = false;
            allClashCodes = [];
            for L = 1:nLayers
                prev  = segVol(mask, L);
                clash = prev ~= opts.bkdCode & prev ~= code;
                if ~any(clash)
                    segVol(mask, L) = code;
                    LUT_found(idx,:) = {code,name};
                    placed = true;
                    if v; fprintf('Wrote code %*d (%*s) to layer %*d from file(s) "%s"\n', ...
                            maxCodeLen,code,maxNameLen,name,maxLayerDigits,L,fileNames); end
                    break
                end
                allClashCodes = [allClashCodes; unique(prev(clash))]; %#ok<AGROW>
            end
        
            if ~placed
                clashCodes = unique(allClashCodes);
                clashNames = strings(size(clashCodes));
                for k = 1:numel(clashCodes)
                    nm = LUT.name(LUT.code == clashCodes(k));
                    if isempty(nm); nm = "(not in LUT)"; end
                    clashNames(k) = nm(1);
                end
                codeW    = max(4, max(strlength(string(clashCodes))));
                header   = "     " + pad("code", codeW, 'left') + "   name";
                clashTbl = strjoin("     " + pad(string(clashCodes(:)), codeW, 'left') + "   " + clashNames(:), newline);
        
                if opts.keepOverlaps
                    outcome = sprintf('No free layer (of %d) for this region.', nLayers);
                else
                    outcome = 'Overlaps are not being kept (opts.keepOverlaps = false).';
                end
                warning(['Issue with code %d (%s) from file(s) "%s":\n' ...
                         '-  Overlaps previously assigned voxels in all %d layer(s):\n' ...
                         '%s\n%s\n' ...
                         '-  %s\n' ...
                         '=> Skipping this region...'], ...
                         code, name, fileNames, nLayers, header, clashTbl, outcome);
            end
        end
    end

    if all(cellfun(@(x) isempty(x),LUT_found))
        %if all failed, set to empty singular
        segVol = [];
        LUT_found = table;
    else
        % drop unused trailing layers
        lastUsed = find(any(segVol ~= opts.bkdCode, 1), 1, 'last');
        segVol = segVol(:, 1:lastUsed);
        
        minCode = min(segVol(:)); maxCode = max(segVol(:));
        if maxCode > intmax('int32') || minCode < intmin('int32')
            error('Codes exceed int32 range.');
        elseif maxCode > intmax('int16') || minCode < intmin('int16')
            warning('Some codes exceed int16 range, using int32...');
            segVol = int32(segVol);
        else
            segVol = int16(segVol);
        end
        segVol = reshape(segVol, [volSize, lastUsed]);
        emptyIdx = all(cellfun(@isempty, LUT_found), 2);
        LUT_found(emptyIdx,:) = [];
        LUT_found = cell2table(LUT_found,"VariableNames",{'code','name'});
    end
end