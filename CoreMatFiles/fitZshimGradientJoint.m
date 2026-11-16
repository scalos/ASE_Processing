function out = fitZshimGradientJoint(yShim, moments, tDeph, sliceThk, opts)
% fitZshimGradientJoint  One through-slice gradient Gz per voxel from the z-shim frames of
% several echoes (the same physical gradient), with a separate amplitude per echo.
%
%   yShim    : 1 x nE cell, each X x Y x Z x nShim_e  (z-shim frames of echo e)
%   moments  : 1 x nE cell, each nShim_e x 1 applied moments (T*s/m), same frame order
%   tDeph    : nE x 1 dephasing time of each echo's shim block (s)
%   sliceThk : nE x 1 slice thickness per echo (m)
%   opts.GzRange, opts.GzStep : search grid (uT/m), as in fitZshimGradient
%
% Model per echo (rect slice profile, unnormalized sinc):
%   y_e(M) = A_e * |sinc( gamma*dz_e*(Gz*t_e - M)/2 )|
% For each grid Gz_j, A_e is solved in closed form per echo and the residual is
% R_j = sum_e [ sum(y_e^2) - (sum(y_e f_ej))^2 / sum(f_ej^2) ] (equal noise per echo assumed),
% minimized over the grid and refined with a parabola, as in fitZshimGradient.
%
% out fields: Gz (X x Y x Z), A / atten / edge (X x Y x Z x nE), R (X x Y x Z, total)
%   edge bit 0 (1) = minimum on the grid edge; bit 1 (2) = peak moment Gz*t_e outside echo e's
%   sampled moment range.
    arguments
        yShim cell
        moments cell
        tDeph (:,1) double {mustBePositive}
        sliceThk (:,1) double {mustBePositive}
        opts.GzRange (1,2) double = [-150 150]
        opts.GzStep (1,1) double {mustBePositive} = 1
    end
    gamma = 2.67502e8;
    nE = numel(yShim);
    assert(numel(moments) == nE && numel(tDeph) == nE && numel(sliceThk) == nE, ...
        'yShim, moments, tDeph and sliceThk must have one entry per echo.');
    sz = size(yShim{1}, 1:3);  nVox = prod(sz);
    grid = (opts.GzRange(1):opts.GzStep:opts.GzRange(2)).';  nG = numel(grid);

    Y = cell(nE,1);  ok = true(nVox,1);  R = zeros(nVox, nG);
    for e = 1:nE
        M = moments{e}(:);
        Y{e} = reshape(double(yShim{e}), nVox, []);
        assert(size(Y{e},2) == numel(M), 'Echo %d: %d moments but %d frames.', e, numel(M), size(Y{e},2));
        ok = ok & all(isfinite(Y{e}), 2) & any(Y{e} > 0, 2);
        c = gamma*sliceThk(e)/2;
        F = abssinc(c*(grid.'*1e-6*tDeph(e) - M));          % nShim_e x nG
        sFF = sum(F.^2, 1);  sFF(sFF < eps) = Inf;
        R = R + sum(Y{e}.^2, 2) - (Y{e}*F).^2 ./ sFF;
    end
    [Rmin, j] = min(R, [], 2);

    off = zeros(nVox,1);
    in = j > 1 & j < nG;
    li = sub2ind([nVox nG], find(in), j(in)-1);
    ci = sub2ind([nVox nG], find(in), j(in));
    ri = sub2ind([nVox nG], find(in), j(in)+1);
    den = R(li) - 2*R(ci) + R(ri);
    o = 0.5*(R(li) - R(ri)) ./ den;
    o(~isfinite(o) | den <= 0) = 0;
    off(in) = max(-0.5, min(0.5, o));
    Gz = grid(j) + off*opts.GzStep;

    A = NaN(nVox, nE);  atten = NaN(nVox, nE);  edge = zeros(nVox, nE, 'uint8');
    onEdge = j == 1 | j == nG;
    for e = 1:nE
        M = moments{e}(:);
        f = abssinc(gamma*sliceThk(e)/2*(Gz*1e-6*tDeph(e) - M.'));     % nVox x nShim_e
        A(:,e) = sum(Y{e}.*f, 2) ./ sum(f.^2, 2);
        z0 = find(M == 0, 1);
        if ~isempty(z0), atten(:,e) = Y{e}(:,z0) ./ A(:,e); end
        pk = Gz*1e-6*tDeph(e);
        out_e = pk < min(M) | pk > max(M);
        edge(onEdge,e) = bitor(edge(onEdge,e), uint8(1));
        edge(out_e,e)  = bitor(edge(out_e,e),  uint8(2));
    end

    bad = ~ok | any(~isfinite(A) | A <= 0, 2);
    Gz(bad) = NaN;  Rmin(bad) = NaN;  A(bad,:) = NaN;  atten(bad,:) = NaN;  edge(bad,:) = 0;
    out = struct('Gz', reshape(Gz, sz), 'R', reshape(Rmin, sz), ...
                 'A', reshape(A, [sz nE]), 'atten', reshape(atten, [sz nE]), 'edge', reshape(edge, [sz nE]));
end

function s = abssinc(x)
    s = abs(sin(x)./x);
    s(x == 0) = 1;
end
