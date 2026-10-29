function aseData = combineAseReps(aseFiles,opts)
    arguments
        aseFiles
        opts.method = 'avg';
    end
    %check if aseFiles have matching geometry:
    isequalNiftiGeometry(aseFiles);
    aseData = niftiread(aseFiles{1});
    % if sorted_ASE_Data is cell array, assume members are averages to be
    % combined:
    switch opts.method
        case 'avg'
            for rep = 2:numel(aseFiles)
                aseData = aseData+niftiread(aseFiles{rep});
            end
            aseData = double(aseData./numel(aseFiles));
        otherwise
            error('Unrecognized combination method');
    end

end