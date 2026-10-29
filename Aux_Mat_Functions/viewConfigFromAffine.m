function views = viewConfigFromAffine(A, display)
% Returns a struct with fields axial/coronal/sagittal, each containing
%   h, v, f   : raw dim on horizontal axis / vertical axis / held fixed
%   xdir,ydir : axis directions to apply AFTER imagesc
%   vox       : voxel size along raw dims (for daspect)
% display: 'radiological' (default; patient L on screen right) or 'neurological'
    arguments
        A
        display {mustBeMember(display,{'radio','neuro'})} = 'radio'
    end
    [~, perm, sgn, vox] = aff2axcodes(A);

    dimOfWorld = zeros(1,3);  signOfWorld = zeros(1,3);
    dimOfWorld(perm)  = 1:3;              % which raw dim follows world axis w
    signOfWorld(perm) = sgn;              % and in which sense

    if strcmp(display, 'radio'), sx = -1; else, sx = +1; end
    % world axes: 1=x (R+), 2=y (A+), 3=z (S+)
    %           name        horiz  screen-right  vert  screen-up  fixed
    spec = {'axial',        1,     sx,            2,    +1,        3;
            'coronal',      1,     sx,            3,    +1,        2;
            'sagittal',     2,     +1,            3,    +1,        1};

    views = struct();
    for n = 1:size(spec,1)
        [name, hw, sr, vw, su, fw] = spec{n,:};
        cfg.h = dimOfWorld(hw);
        cfg.v = dimOfWorld(vw);
        cfg.f = dimOfWorld(fw);
        % index grows toward screen-right (up) when sense * desired direction > 0
        cfg.xdir = tern(signOfWorld(hw) * sr > 0, 'normal', 'reverse');
        cfg.ydir = tern(signOfWorld(vw) * su > 0, 'normal', 'reverse');
        cfg.vox  = vox;
        views.(name) = cfg;
    end
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end