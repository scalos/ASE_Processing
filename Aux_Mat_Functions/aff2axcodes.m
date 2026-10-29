function [codes, perm, sgn, vox] = aff2axcodes(A)
% codes : 1x3 char, the world direction each raw array axis points TOWARD
%         (e.g. 'LIA' means i->Left, j->Inferior, k->Anterior)
% perm  : perm(d) = world axis (1=x, 2=y, 3=z) that raw dim d follows
% sgn   : sgn(d)  = +1 if raw dim d increases along that world axis, else -1
% vox   : voxel size along each raw dim (mm)
    R   = A(1:3,1:3);
    vox = sqrt(sum(R.^2, 1));
    R   = R ./ vox;                       % unit direction of each voxel axis

    M    = abs(R);
    perm = zeros(1,3);
    sgn  = zeros(1,3);
    dom  = zeros(1,3);
    for n = 1:3                           % greedy assignment -> valid permutation
        [~, lin] = max(M(:));
        [w, d]   = ind2sub([3 3], lin);
        perm(d)  = w;
        sgn(d)   = sign(R(w,d)) + (R(w,d) == 0);
        dom(d)   = abs(R(w,d));
        M(w,:)   = -1;
        M(:,d)   = -1;
    end
    if any(dom < 0.99)
        warning('aff2axcodes:oblique', ...
            ['Volume is oblique (min alignment %.3f). Axis codes are only ' ...
             'approximate; slice display will be slightly off-angle.'], min(dom));
    end

    labels = ['L' 'R'; 'P' 'A'; 'I' 'S'];    % row = world axis, col = (-,+)
    codes  = char(zeros(1,3));
    for d = 1:3
        codes(d) = labels(perm(d), (sgn(d) > 0) + 1);
    end
end