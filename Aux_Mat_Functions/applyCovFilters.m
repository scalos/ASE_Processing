function [keep, desc] = applyCovFilters(T, filters)
% Selects rows of a covariate table (from readCovTable) with simple column filters.
%
%   filters : string array / cellstr; every filter must hold (AND). Forms:
%       "col=a,b,c"    value is one of a, b, c   (text or numeric columns)
%       "col!=a,b"     value is none of a, b
%       "col>x", "col>=x", "col<x", "col<=x"     numeric columns only
%   Rows with a missing value in a filtered column fail that filter.
%
%   keep : logical, one per row of T
%   desc : one line per filter: "<filter>: <n> of <N> rows pass"
    arguments
        T table
        filters = strings(0,1)
    end
    filters = string(filters);  filters = filters(strlength(strtrim(filters)) > 0);
    keep = true(height(T),1);
    desc = strings(0,1);
    for k = 1:numel(filters)
        f = strtrim(filters(k));
        tok = regexp(f, '^\s*([^!<>=\s]+)\s*(!=|>=|<=|=|>|<)\s*(.*)$', 'tokens', 'once');
        assert(~isempty(tok), ['Cannot parse filter "%s". Use col=a,b  col!=a  col>x  col>=x  ' ...
            'col<x  col<=x.'], f);
        col = string(tok{1});  op = string(tok{2});  rhs = strtrim(string(tok{3}));
        vn = string(T.Properties.VariableNames);
        hit = find(strcmpi(vn, col) | strcmpi(vn, matlab.lang.makeValidName(col)), 1);
        assert(~isempty(hit), 'Filter "%s": column "%s" not in the covariate table (columns: %s).', ...
            f, col, strjoin(vn, ', '));
        x = T.(vn(hit));
        switch op
            case {"=","!="}
                vals = strtrim(split(rhs, ","));
                if isnumeric(x)
                    v = str2double(vals);
                    assert(all(~isnan(v)), 'Filter "%s": column "%s" is numeric.', f, col);
                    in = ismember(x, v);
                else
                    in = ismember(x, vals);
                end
                ok = ~ismissing(x);
                if op == "=", pass = in & ok; else, pass = ~in & ok; end
            otherwise
                assert(isnumeric(x), 'Filter "%s": %s needs a numeric column.', f, op);
                v = str2double(rhs);
                assert(~isnan(v), 'Filter "%s": "%s" is not a number.', f, rhs);
                switch op
                    case ">",  pass = x >  v;
                    case ">=", pass = x >= v;
                    case "<",  pass = x <  v;
                    case "<=", pass = x <= v;
                end
                pass = pass & ~isnan(x);
        end
        keep = keep & pass;
        desc(end+1,1) = sprintf('%s: %d of %d rows pass', f, nnz(pass), height(T)); %#ok<AGROW>
    end
end
