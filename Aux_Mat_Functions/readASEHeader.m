function dTEs = readASEHeader(hdr_File)
    hdrTable = readtable(hdr_File);
    dTE_col = hdrTable.Var4;

    dTEs = nan(size(dTE_col));
    for i = 1:numel(dTE_col)
        tok = regexp(dTE_col{i}, '(?<![A-Za-z])DTE=(\d+\.?\d*)', 'tokens');
        if ~isempty(tok)
            dTEs(i) = str2double(tok{1}{1});
        else
            warning('No DTE found in string %d: %s', i, dTE_col{i});
        end
    end
end