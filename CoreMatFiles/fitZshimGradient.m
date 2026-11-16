function out = fitZshimGradient(yShim, moments, tDeph, sliceThk, opts)
% fitZshimGradient  Voxelwise through-slice background gradient (Gz) from z-shim frames.
%
%   yShim    : X x Y x Z x nShim  signal of the z-shim frames of ONE echo
%   moments  : nShim x 1          applied z-shim moments (T*s/m), same order as yShim frames
%   tDeph    : scalar             dephasing time of the shim block (s), i.e. its |tau| entry
%   sliceThk : scalar             slice thickness (m)
%   opts.GzRange : [lo hi] search range (uT/m), default [-150 150]
%   opts.GzStep  : grid step (uT/m), default 1
%
% Model (rect slice profile; UNNORMALIZED sinc, sinc(x) = sin(x)/x):
%   y(M) = A * |sinc( gamma*dz*(Gz*tDeph - M)/2 )|
% For each grid value Gz_j the template f_j is fixed, so A_j = sum(y f_j)/sum(f_j^2) and
% R_j = sum((y - A_j f_j)^2) = sum(y^2) - (sum(y f_j))^2/sum(f_j^2), evaluated for all voxels
% at once. The best j is refined with a parabola through R(j-1), R(j), R(j+1).
%
% out fields (X x Y x Z):
%   Gz    : estimated gradient (uT/m), NaN where not estimable
%   A     : fitted amplitude at Gz
%   R     : residual sum of squares at the grid minimum
%   atten : model-free attenuation ratio y(M=0)/A (NaN if no zero-moment frame)
%   edge  : uint8 bitmask; bit 0 (1) = minimum on the grid edge,
%           bit 1 (2) = peak moment Gz*tDeph outside the sampled moment range
    arguments
        yShim {mustBeNumeric}
        moments (:,1) double
        tDeph (1,1) double {mustBePositive}
        sliceThk (1,1) double {mustBePositive}
        opts.GzRange (1,2) double = [-150 150]
        opts.GzStep (1,1) double {mustBePositive} = 1
    end
    gamma = 2.67502e8;                         % rad/s/T
    sz = size(yShim); sz(end+1:4) = 1;
    nShim = sz(4);
    assert(numel(moments) == nShim, 'moments has %d entries but yShim has %d frames.', numel(moments), nShim);
    nVox = prod(sz(1:3));
    Y = reshape(double(yShim), nVox, nShim);
    ok = all(isfinite(Y), 2) & any(Y > 0, 2);

    grid = (opts.GzRange(1):opts.GzStep:opts.GzRange(2)).';   % uT/m
    nG = numel(grid);
    c = gamma*sliceThk/2;
    % templates: nShim x nG
    F = abssinc(c*(grid.'*1e-6*tDeph - moments));
    sFF = sum(F.^2, 1);                                     % 1 x nG
    sFF(sFF < eps) = Inf;                                   % degenerate template -> no fit
    sYY = sum(Y.^2, 2);
    YF  = Y * F;                                            % nVox x nG
    R   = sYY - YF.^2 ./ sFF;                               % nVox x nG
    [Rmin, j] = min(R, [], 2);

    % parabolic refinement (interior points only)
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

    % amplitude at the refined Gz
    f = abssinc(c*(Gz*1e-6*tDeph - moments.'));            % nVox x nShim
    A = sum(Y.*f, 2) ./ sum(f.^2, 2);

    % model-free attenuation: zero-moment frame / fitted amplitude
    z0 = find(moments == 0, 1);
    if isempty(z0)
        atten = NaN(nVox,1);
    else
        atten = Y(:,z0) ./ A;
    end

    edge = zeros(nVox,1,'uint8');
    edge(j == 1 | j == nG) = bitor(edge(j == 1 | j == nG), uint8(1));
    pk = Gz*1e-6*tDeph;                                     % moment at the profile peak
    outside = pk < min(moments) | pk > max(moments);
    edge(outside) = bitor(edge(outside), uint8(2));

    bad = ~ok | ~isfinite(A) | A <= 0;
    Gz(bad) = NaN;  A(bad) = NaN;  Rmin(bad) = NaN;  atten(bad) = NaN;  edge(bad) = 0;

    out = struct('Gz',   reshape(Gz,   sz(1:3)), ...
                 'A',    reshape(A,    sz(1:3)), ...
                 'R',    reshape(Rmin, sz(1:3)), ...
                 'atten',reshape(atten,sz(1:3)), ...
                 'edge', reshape(edge, sz(1:3)));
end

function s = abssinc(x)
% |sin(x)/x| (unnormalized sinc; NOT MATLAB's sinc, which is sin(pi x)/(pi x))
    s = abs(sin(x)./x);
    s(x == 0) = 1;
end
