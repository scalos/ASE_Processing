function T = anLin2003DesignFit(sessDir, gzFile, asegFile, opts)
% anLin2003DesignFit  Refit one subject with the An & Lin 2003 (MRM 50:708) ASE analysis.
%
% Their design: one echo; images with identical tau averaged; optional 3x3 in-plane low-pass;
% linear fit of ln S over dephasing times 20-40 ms (slope = R2'); vCBV = extrapolated line at
% 0 minus ln(mean of the spin-echo images); OEF = R2'/(vCBV * (4/3)*pi*gamma*dChi0*Hct*B0).
% This lets the dataset be compared with the paper (GM ~4%, WM ~3-3.5% vCBV with flow
% suppression; ~5% / ~3.6% without) independently of the pipeline's two-echo regime design.
%
%   sessDir  : session_ase folder (ase_hdr.json, ase_e1_proc.nii.gz, xvent mask on ASE grid)
%   gzFile   : echo-1 Gz map (uT/m) from a pipeline run with quality.zshimEstimate ([] = none)
%   asegFile : aparc+aseg resampled to the ASE grid (NN)
%   opts.longRange : dephasing-time window for the linear fit (ms), default [20 40]
%   opts.Fcut      : z-shim frames with F below this are dropped (default 0.1)
%
% Variants (rows of T): raw, smoothed, z-shim corrected (F^2-weighted, as the pipeline's
% dBCorr), and both. Columns: voxel medians of vCBV / R2' / OEF in clean deep WM, cerebral WM
% and cortex (tissue fraction > 0.9).
% Needs CoreMatFiles and Aux_Mat_Functions on the path.
    arguments
        sessDir (1,:) char
        gzFile = []
        asegFile (1,:) char = ''
        opts.longRange (1,2) double = [20 40]
        opts.Fcut (1,1) double = 0.1
        opts.Hct (1,1) double = 0.42*0.85
        opts.dChi0 (1,1) double = 0.18e-6
        opts.verbose (1,1) logical = true
    end
    [dTEs, ~, ~, info] = readASEHeader(fullfile(sessDir, 'ase_hdr.json'));
    Y = double(niftiread(fullfile(sessDir, 'ase_e1_proc.nii.gz')));
    sweep = info.repInds;                                  % nReps x nSweep raw-frame indices
    tau = dTEs(sweep(1,:));  tau = tau(:);                 % ms, identical across reps
    Ym = zeros([size(Y,1:3) numel(tau)]);
    for r = 1:size(sweep,1), Ym = Ym + Y(:,:,:,sweep(r,:)); end
    Ym = Ym / size(sweep,1);                               % average identical-tau images
    oefScale = (4/3)*pi*2.67502e8*opts.dChi0*opts.Hct*info.B0;
    dz = info.SliceThickness_mm(1)/1000;

    Gz = [];
    if ~isempty(gzFile), Gz = double(niftiread(gzFile)); end
    aseg = double(niftiread(asegFile));
    tis  = double(niftiread(fullfile(sessDir, 'anat_fs_xvent_mask_al_ase.nii.gz')));
    wm   = ismember(aseg,[2 41]) & tis > 0.9;
    clean = imerode(ismember(aseg,[2 41]), ones(3,3,3)) & tis > 0.999;
    ctx  = aseg >= 1000 & tis > 0.9;

    lg = tau >= opts.longRange(1) & tau <= opts.longRange(2);
    se = tau == 0;
    variants = {'raw',false,false; 'smoothed 3x3',true,false; 'z-shim corr',false,true; 'z-shim + smoothed',true,true};
    rows = {};
    for v = 1:size(variants,1)
        if variants{v,3} && isempty(Gz), continue; end
        S = Ym;  W = [];
        if variants{v,3}
            W = ones(size(S));
            for f = 1:numel(tau)
                F = zshimAttenuation(Gz, tau(f)/1000, dz);  F(isnan(F)) = 1;
                fr = S(:,:,:,f) ./ F;  fr(F < opts.Fcut) = NaN;
                S(:,:,:,f) = fr;  W(:,:,:,f) = F.^2;
            end
        end
        if variants{v,2}                                     % 3x3 in-plane low-pass
            k = ones(3,3)/9;
            for f = 1:numel(tau)
                for z = 1:size(S,3), S(:,:,z,f) = conv2(S(:,:,z,f), k, 'same'); end
            end
        end
        wl = [];  if ~isempty(W), wl = W(:,:,:,lg); end
        B = fitMultiRegressor(S(:,:,:,lg), [ones(nnz(lg),1), -tau(lg)/1000], [], [], wl);
        lnSE = log(mean(S(:,:,:,se), 4, 'omitnan'));
        vC = B(:,:,:,1) - lnSE;  R2P = B(:,:,:,2);
        O = R2P ./ (vC*oefScale);  O(vC <= 0 | R2P <= 0) = NaN;
        m = @(x,k) median(x(k), 'omitnan');
        rows(end+1,:) = {variants{v,1}, m(vC,clean), m(vC,wm), m(vC,ctx), m(R2P,clean), m(R2P,wm), m(R2P,ctx), ...
                         m(O,clean), m(O,wm), m(O,ctx)}; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'variant','vCBV_cleanWM','vCBV_WM','vCBV_ctx', ...
        'R2P_cleanWM','R2P_WM','R2P_ctx','OEF_cleanWM','OEF_WM','OEF_ctx'});
    if opts.verbose
        fprintf('An & Lin 2003 design: echo 1, long %g-%g ms (%d taus), %d SE images averaged over %d reps\n', ...
            opts.longRange, nnz(lg), nnz(se), size(sweep,1));
        disp(T);
    end
end
