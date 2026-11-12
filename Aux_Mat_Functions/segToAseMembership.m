function M = segToAseMembership(segVol, LUT, segInfo, aseInfo, opts)
% Fractional-occupancy membership matrix on the ASE grid.
%   segVol   : Xs x Ys x Zs x nLayers label volume (high-res grid)
%   LUT      : table with .code
%   segInfo, aseInfo : niftiinfo structs (same world space assumed)
%   opts.minFrac  : drop weights below this occupancy (e.g. 0.25), 0 to keep all
%   opts.maxSampling : max number of voxels in each spatial domain to look
%   through when computing occupancy.
% Returns sparse M (nReg x nVox), entries in [0,1] = fraction of ASE voxel covered.
    arguments
        segVol
        LUT
        segInfo
        aseInfo
        opts.minFrac (1,1) = 0;
        opts.maxSampling (1,3) = [5,5,5];
    end
    minFrac = opts.minFrac;
    assert(0<=minFrac & minFrac <=1,'ERROR: minFrac must be between 0 and 1');
    aseSz = aseInfo.ImageSize(1:3);
    nReg = height(LUT);  nVox = prod(aseSz);
    
    %MODIFY Tc here to include tForm
    % ASE index -> seg index (row-vector convention: idx = [i j k 1]*T)
    Tc = aseInfo.Transform.T/segInfo.Transform.T;

    % number of sub-samples per axis from the voxel-size ratio (capped)
    vA = sqrt(sum(aseInfo.Transform.T(1:3,1:3).^2, 2))';
    vS = sqrt(sum(segInfo.Transform.T(1:3,1:3).^2, 2))';
    nSub = max(1, min(opts.maxSampling, ceil(vA ./ vS)));

    [I,J,K] = ndgrid(1:aseSz(1), 1:aseSz(2), 1:aseSz(3));
    I = I(:); J = J(:); K = K(:);
    voxIdx = (1:nVox)';
    segSz = size(segVol,1:3);
    nPts  = prod(nSub);

    M = sparse(nReg, nVox);
    oi = ((1:nSub(1)) - 0.5)/nSub(1) - 0.5;
    oj = ((1:nSub(2)) - 0.5)/nSub(2) - 0.5;
    ok = ((1:nSub(3)) - 0.5)/nSub(3) - 0.5;
    for a = oi
      for b = oj
        for c = ok
            P = [I+a, J+b, K+c, ones(nVox,1)] * Tc;           % seg-grid coordinates
            p = round(P(:,1:3));
            inb = all(p >= 1 & p <= segSz, 2);
            lin = sub2ind(segSz, p(inb,1), p(inb,2), p(inb,3));
            vox = voxIdx(inb);
            for L = 1:size(segVol,4)
                lab = segVol(:,:,:,L);
                [tf,loc] = ismember(lab(lin), LUT.code);
                M = M + sparse(loc(tf), vox(tf), 1/nPts, nReg, nVox);
            end
        end
      end
    end
    if minFrac > 0
        M(M < minFrac) = 0;     % drops low-occupancy voxels from a region
    end
end