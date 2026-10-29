function [path,name,ext] = niftiFileparts(path)
    [path,name,ext] = fileparts(path);
    [~,name2,ext2] = fileparts(name);
    if ~isempty(ext2)
        name = name2;
        ext = sprintf("%s%s",ext2,ext);
    end
end