function q = bhFDR(p)
% Benjamini-Hochberg adjusted p-values (q-values). NaN entries are ignored and stay NaN;
% the number of tests is the number of finite p-values.
    q = NaN(size(p));
    ok = isfinite(p);
    pv = p(ok);  m = numel(pv);
    if m == 0, return; end
    [ps, order] = sort(pv(:));
    qs = ps .* m ./ (1:m).';
    qs = min(1, flipud(cummin(flipud(qs))));     % enforce monotonicity from the top
    out = zeros(m,1);  out(order) = qs;
    q(ok) = out;
end
