function goodPath = findGoodPath(path,overwrite)
    arguments
        path
        overwrite = false
    end
    goodPath = path;
    if isfile(path)
        if overwrite
            delete(path);
            return;
        end
        goodPath = [];
        for idx = 1:100
            [path_,name,ext] = niftiFileparts(path);
            goodPath = fullfile(path_,sprintf("%s_%d%s",name,idx,ext));
            if ~isfile(goodPath)
                return
            end
        end
    end
end