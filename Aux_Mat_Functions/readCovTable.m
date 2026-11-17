function [T, info] = readCovTable(covFile)
% Reads a whitespace-delimited covariate table (e.g. cov_table_ase).
%
%   First line: header, optionally starting with '#' (e.g. "#subjinit sex age bmi group").
%   Other lines: one subject per line; the first column is the subject ID.
%   Lines starting with '#' after the header are commented-out subjects: they are skipped and
%   listed in info.commented. Blank lines are ignored. "NA", "NaN" and "" are missing values.
%
%   T    : table; T.subject (string) plus one variable per remaining header column. A column whose
%          non-missing entries are all numeric becomes double (missing = NaN); any other column
%          stays string (missing = <missing>).
%   info : struct with header (original names), commented (IDs of commented-out rows), file
    arguments
        covFile {mustBeTextScalar}
    end
    covFile = char(covFile);
    assert(isfile(covFile), 'Covariate table not found: %s', covFile);
    lines = splitlines(string(fileread(covFile)));
    lines = strtrim(erase(lines, char(13)));
    lines = lines(strlength(lines) > 0);
    assert(~isempty(lines), 'Covariate table is empty: %s', covFile);

    header = split(strtrim(regexprep(lines(1), '^#+', ''))).';
    nCol = numel(header);
    assert(nCol >= 1 && strlength(header(1)) > 0, 'Covariate table %s has no header.', covFile);

    body = lines(2:end);
    isComment = startsWith(body, '#');
    commented = strings(0,1);
    for k = find(isComment).'
        f = split(strtrim(regexprep(body(k), '^#+', '')));
        if ~isempty(f) && strlength(f(1)) > 0, commented(end+1,1) = f(1); end %#ok<AGROW>
    end
    body = body(~isComment);

    cells = strings(numel(body), nCol);
    for k = 1:numel(body)
        f = split(body(k)).';
        assert(numel(f) == nCol, ...
            'Covariate table %s: line for "%s" has %d columns, header has %d.', ...
            covFile, f(1), numel(f), nCol);
        cells(k,:) = f;
    end

    ids = cells(:,1);
    [u,~,j] = unique(ids);
    dup = u(accumarray(j,1) > 1);
    assert(isempty(dup), 'Covariate table %s lists subject(s) more than once: %s', ...
        covFile, strjoin(dup, ', '));

    varNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(header));
    varNames(1) = "subject";
    T = table(ids, 'VariableNames', "subject");
    for c = 2:nCol
        col = cells(:,c);
        miss = ismember(upper(col), ["NA","NAN",""]);
        num = str2double(col);
        if all(~isnan(num) | miss)
            num(miss) = NaN;
            T.(varNames(c)) = num;
        else
            col(miss) = missing;
            T.(varNames(c)) = col;
        end
    end
    info = struct('header', header, 'varNames', varNames, 'commented', commented, 'file', string(covFile));
end
