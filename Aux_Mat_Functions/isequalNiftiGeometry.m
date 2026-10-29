function [ok, report] = isequalNiftiGeometry(vols, opts)
% isequalNiftiGeometry  Check that NIfTI volumes share the same spatial geometry.
%   vols : cell array of file paths and/or niftiinfo structs
%   Compares (first 3 dims only, so 3D/4D volumes are comparable):
%     - ImageSize
%     - PixelDimensions (voxel size)
%     - Transform.T (full affine: orientation, origin, and voxel-to-world mapping)
%     - SpaceUnits
%   Returns ok (logical) and report (cell array of mismatch descriptions).
    arguments
        vols
        opts.tol = 1e-4                 % tolerance for voxel size / affine
        opts.onMismatch {mustBeMember(opts.onMismatch,{'warn','error','silent'})} = 'warn'
    end
    if ~iscell(vols); vols = {vols}; end
    n = numel(vols);
    infos = cell(n,1);
    names = cell(n,1);
    for i = 1:n
        if isstruct(vols{i})
            infos{i} = vols{i};
            names{i} = infos{i}.Filename;
        else
            infos{i} = niftiinfo(vols{i});
            names{i} = vols{i};
        end
    end

    report = {};
    ref = infos{1};
    refSize = pad3(ref.ImageSize);
    refPix  = pad3(ref.PixelDimensions);
    refT    = ref.Transform.T;

    for i = 2:n
        cur = infos{i};
        tag = sprintf('"%s" vs reference "%s"', names{i}, names{1});

        if ~isequal(pad3(cur.ImageSize), refSize)
            report{end+1} = sprintf('%s: matrix size differs [%s] vs [%s]', tag, ...
                num2str(pad3(cur.ImageSize)), num2str(refSize)); %#ok<AGROW>
        end
        if any(abs(pad3(cur.PixelDimensions) - refPix) > opts.tol)
            report{end+1} = sprintf('%s: voxel size differs [%s] vs [%s]', tag, ...
                num2str(pad3(cur.PixelDimensions)), num2str(refPix)); %#ok<AGROW>
        end
        if max(abs(cur.Transform.T - refT),[],'all') > opts.tol
            report{end+1} = sprintf('%s: affine transform (orientation/origin) differs', tag); %#ok<AGROW>
        end
        if ~strcmp(cur.SpaceUnits, ref.SpaceUnits)
            report{end+1} = sprintf('%s: space units differ (%s vs %s)', tag, ...
                cur.SpaceUnits, ref.SpaceUnits); %#ok<AGROW>
        end
    end

    ok = isempty(report);
    if ~ok
        msg = strjoin(report, newline);
        switch opts.onMismatch
            case 'error'; error('isequalNiftiGeometry:mismatch', ...
                            'Nifti geometry mismatch:\n%s', msg);
            case 'warn';  warning('isequalNiftiGeometry:mismatch', ...
                            'Nifti geometry mismatch:\n%s', msg);
        end
    end
end

function v = pad3(v)
    v = v(:)';
    v(end+1:3) = 1;
    v = v(1:3);
end