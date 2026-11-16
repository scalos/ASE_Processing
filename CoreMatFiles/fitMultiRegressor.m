function [B, rms, nValid] = fitMultiRegressor(ydata, Xdes, maskArray, chunkSize, weights)
% Voxelwise (optionally weighted) least squares of log(ydata) against the columns of Xdes, with per-voxel
% exclusion of invalid frames (NaN, <=0, Inf).
%
% ydata     : Nx x Ny x Nz x Nd x Nr   (NaN / <=0 = excluded)
% Xdes      : (Nd*Nr) x p  or  Nd x p  design matrix (include ones column
%             for an intercept). Row order: dTE fastest, then repetition.
% maskArray : optional Nx x Ny x Nz logical; only these voxels are fitted
% chunkSize : optional, voxels per chunk (default 1e5)
% weights   : optional, same size as ydata: relative inverse log-space variance per
%             voxel and frame (e.g. F^2 after z-shim correction). Empty = unweighted OLS.
%             NaN weights count as 1.
%
% B      : Nx x Ny x Nz x p   coefficients (NaN where fit not possible)
% rms    : Nx x Ny x Nz       sqrt(sum(w.*resid.^2)/n) in log space (divides by n, not n-p);
%                             plain RMS when unweighted, so rms.^2.*n is the (weighted) RSS
% nValid : Nx x Ny x Nz       number of frames used

    if nargin < 3, maskArray = []; end
    if nargin < 4 || isempty(chunkSize), chunkSize = 1e5; end
    if nargin < 5, weights = []; end

    sz = size(ydata); sz(end+1:5) = 1;
    Nx = sz(1); Ny = sz(2); Nz = sz(3); Nd = sz(4); Nr = sz(5);
    Nv = Nx*Ny*Nz;  Nfr = Nd*Nr;  p = size(Xdes,2);

    if size(Xdes,1) == Nfr
        % use as is
    elseif size(Xdes,1) == Nd
        Xdes = repmat(Xdes, Nr, 1);
    else
        error('Xdes must have Nd (%d) or Nd*Nr (%d) rows.', Nd, Nfr);
    end

    % Column scaling improves conditioning of the normal equations
    colScale = sqrt(mean(Xdes.^2, 1));
    colScale(colScale == 0) = 1;
    Xs = bsxfun(@rdivide, Xdes, colScale);

    Y2 = reshape(ydata, Nv, Nfr);
    useW = ~isempty(weights);
    if useW
        assert(numel(weights) == numel(ydata), 'weights must have the same size as ydata.');
        Wt2 = reshape(double(weights), Nv, Nfr);
        Wt2(~isfinite(Wt2)) = 1;
    end
    if isempty(maskArray)
        vox = (1:Nv)';
    else
        vox = find(maskArray(:));
    end

    B2   = NaN(Nv, p);
    rms2 = NaN(Nv, 1);
    n2   = zeros(Nv, 1);
    tol  = 1e-10;                         % relative pivot tolerance

    for s = 1:chunkSize:numel(vox)
        idx = vox(s : min(s+chunkSize-1, numel(vox)));
        Y   = Y2(idx, :);
        W   = (Y > 0) & isfinite(Y);
        L   = zeros(size(Y));
        L(W) = log(Y(W));
        Wd  = double(W);
        n   = sum(Wd, 2);
        if useW
            Wd = Wd .* Wt2(idx,:);        % weight 0 on invalid frames
            L  = L .* Wd;                 % b = X'*W*log(y)
        end
        nv  = numel(idx);

        % ---- normal equations A*beta = b (one system per voxel) ----
        A = cell(p,p);  b = cell(p,1);
        for i = 1:p
            b{i} = L * Xs(:,i);
            for j = i:p
                A{i,j} = Wd * (Xs(:,i) .* Xs(:,j));
            end
        end

        % ---- Cholesky A = Lc*Lc' (elementwise across voxels) ----
        Lc  = cell(p,p);
        bad = (n < p);
        for j = 1:p
            d = A{j,j};
            for k = 1:j-1
                d = d - Lc{j,k}.^2;
            end
            ok = d > tol * A{j,j};
            bad = bad | ~ok;
            d(~ok) = 1;                   % placeholder, voxel flagged bad
            Lc{j,j} = sqrt(d);
            for i = j+1:p
                t = A{j,i};
                for k = 1:j-1
                    t = t - Lc{i,k} .* Lc{j,k};
                end
                Lc{i,j} = t ./ Lc{j,j};
            end
        end

        % ---- forward then back substitution ----
        z = cell(p,1);
        for i = 1:p
            t = b{i};
            for k = 1:i-1
                t = t - Lc{i,k} .* z{k};
            end
            z{i} = t ./ Lc{i,i};
        end
        beta = zeros(nv, p);
        for i = p:-1:1
            t = z{i};
            for k = i+1:p
                t = t - Lc{k,i} .* beta(:,k);
            end
            beta(:,i) = t ./ Lc{i,i};
        end
        beta(bad,:) = NaN;

        % ---- residual RMS in log space ----
        if useW
            Lr = zeros(size(Y));  Lr(W) = log(Y(W));
            resid = (beta * Xs.' - Lr) .* double(W);
            r = sqrt(sum(Wd .* resid.^2, 2) ./ n);
        else
            resid = (beta * Xs.' - L) .* Wd;
            r = sqrt(sum(resid.^2, 2) ./ n);
        end
        r(bad) = NaN;

        B2(idx,:) = bsxfun(@rdivide, beta, colScale);   % undo column scaling
        rms2(idx) = r;
        n2(idx)   = n;
    end

    B      = reshape(B2,   Nx, Ny, Nz, p);
    rms    = reshape(rms2, Nx, Ny, Nz);
    nValid = reshape(n2,   Nx, Ny, Nz);
end