function s = compressIntArr(v,opts)
    arguments
        v
        opts.joinChar = ',';
    end
    v = v(:).';                          % force row vector
    breaks = [0, find(diff(v) ~= 1), numel(v)];  % run boundaries
    parts = {};
    for k = 1:numel(breaks)-1
        run = v(breaks(k)+1 : breaks(k+1));
        if numel(run) >= 3
            parts{end+1} = sprintf('%d:%d', run(1), run(end)); %#ok
        else
            parts{end+1} = strjoin(compose('%d', run), ','); %#ok
        end
    end
    s = strjoin(parts, opts.joinChar);
end