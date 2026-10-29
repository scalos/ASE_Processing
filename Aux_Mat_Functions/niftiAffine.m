function A = niftiAffine(raw)
% 4x4 affine mapping 0-based voxel indices [i j k 1]' to world RAS+ mm.
% raw = info.raw from niftiinfo. Prefers sform, falls back to qform.
    if raw.sform_code > 0
        A = [raw.srow_x; raw.srow_y; raw.srow_z; 0 0 0 1];
    elseif raw.qform_code > 0
        b = raw.quatern_b; c = raw.quatern_c; d = raw.quatern_d;
        a = sqrt(max(0, 1 - b^2 - c^2 - d^2));
        R = [a^2+b^2-c^2-d^2, 2*(b*c-a*d),       2*(b*d+a*c);
             2*(b*c+a*d),     a^2+c^2-b^2-d^2,   2*(c*d-a*b);
             2*(b*d-a*c),     2*(c*d+a*b),       a^2+d^2-c^2-b^2];
        qfac = raw.pixdim(1);
        if qfac < 0, R(:,3) = -R(:,3); end
        A = eye(4);
        A(1:3,1:3) = R * diag(raw.pixdim(2:4));
        A(1:3,4)   = [raw.qoffset_x; raw.qoffset_y; raw.qoffset_z];
    else
        error('No sform or qform. Orientation is unknown (method 1 legacy).');
    end
end