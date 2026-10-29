function [data,info] = loadNifti(niiPath,opts)
    arguments
        niiPath
        opts.niiName = "nifti file";
        opts.failType {mustBeMember(opts.failType,{'warn','error','silent'})} = 'warn';
    end
    
    try
        data = niftiread(niiPath);
    catch ME
        data = [];
        mes = sprintf('Unable to load %s data.\n@: %s\nME: (%s;%s)', ...
                    opts.niiName,niiPath,ME.identifier,ME.message);
        switch opts.failType
            case 'warn'
                warning(mes);
            case 'error'
                error(mes);
        end
    end

    try
        info = niftiinfo(niiPath);
    catch ME
        info = [];
        mes = sprintf('Unable to load %s info.\n@: %s\nME: (%s;%s)', ...
                    opts.niiName,niiPath,ME.identifier,ME.message);
        switch opts.failType
            case 'warn'
                warning(mes);
            case 'error'
                error(mes);
        end
    end

end