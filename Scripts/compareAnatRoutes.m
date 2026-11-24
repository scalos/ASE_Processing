function S = compareAnatRoutes(sessDir, afniDir, fsDir, lutFile, outDir, opts)
% compareAnatRoutes  Compare the two ways of bringing FreeSurfer anatomy onto the ASE grid.
%
%   Route A (current): FreeSurfer -> session_afni_nortc ("_al" = aligned to the o2co2 BOLD base by
%     clusterFiles/batch_afni_o2co2_kwong_240114_rs_nortc_ForSam: zero-pad, centre shift, 12-parameter
%     affine to the deobliqued BOLD frame 0) -> asePreProcess stage 6 registers ase_base to
%     anat_fs_brain_al (12 parameters) -> labels resampled (NN) through every step.
%   Route B (proposed): FreeSurfer only. Target = T1 x step(aseg), unifized (the recipe of the old
%     target minus the BOLD step); ase_base registered to it with 6, 9 or 12 parameters; labels
%     resampled once.
%
%   S = compareAnatRoutes(sessDir, afniDir, fsDir, lutFile, outDir)
%     sessDir : session_ase of the subject (ase_base, ase_e*_proc, ase_hdr.json, ase_base_al_*mat,
%               anat_fs_*_al[_ase] files from the current pipeline)
%     afniDir : session_afni_nortc (anat_fs_shft.1D, anat_fs_brain_al_mat.aff12.1D, anat_fs_aparc_aseg.nii,
%               anat_fs_aparc_aseg_al+orig, anat_fs_orig.nii) -- used to rebuild route A end to end
%     fsDir   : FreeSurfer volumes converted with mri_convert: T1, aseg, aparc+aseg, wmparc,
%               brainstemSsLabels.v13.FSvoxelSpace (.nii.gz)
%     lutFile : ROI LUT (GM_WM_long.txt)
%     outDir  : figures, summary.txt, summary.mat
%   opts.dofs    : parameter counts for route B (default [6 9 12])
%   opts.mainDof : route B variant compared in the ROI figures (default 12)
%   opts.doFits  : run the stats-only ROI fits (default true)
%
% Every transform is applied by AFNI itself (3dAllineate -1Dmatrix_apply); the route-A composite
% (ASE <- FreeSurfer) is cat_matvec of the stored matrices and is checked against the pipeline's own
% ASE-space segmentation before it is used.
% Needs AFNI (~/abin or on PATH), CoreMatFiles and Aux_Mat_Functions on the MATLAB path.
    arguments
        sessDir (1,:) char
        afniDir (1,:) char
        fsDir (1,:) char
        lutFile (1,:) char
        outDir (1,:) char
        opts.workDir (1,:) char = tempname
        opts.dofs (1,:) double = [6 9 12]
        opts.mainDof (1,1) double = 12
        opts.doFits (1,1) logical = true
        opts.allineateOpts (1,:) char = '-twopass -cost lpa+ -autoweight -fineblur 1 -cmass'
        opts.reuse (1,1) logical = false
    end
    W = opts.workDir;  afni('reuse', opts.reuse);  afni('dir', W);
    if ~isfolder(outDir), mkdir(outDir); end
    if ~isfolder(fullfile(W,'in')), mkdir(fullfile(W,'in')); end
    logLines = {};
    function note(varargin)
        msg = sprintf(varargin{:}); logLines{end+1} = msg; fprintf('%s\n', msg);
    end
    P_ = @(f) fullfile(W, f);
    C.A = [235 104 52]/255;  C.B = [42 120 214]/255;  C.ref = [150 149 145]/255;
    C.ink = [11 11 11]/255;  C.ink2 = [82 81 78]/255;  C.B12 = [98 80 214]/255;
    seqMap = interp1([0 0.5 1], [255 255 255; 42 120 214; 16 42 92]/255, linspace(0,1,256));

    % ---- inputs, linked into a folder without spaces (AFNI does not take them) ----
    L = struct( ...
        'ase_base', {fullfile(sessDir,'ase_base.nii.gz')}, 'e1', {fullfile(sessDir,'ase_e1_proc.nii.gz')}, ...
        'e2', {fullfile(sessDir,'ase_e2_proc.nii.gz')}, 'hdr', {fullfile(sessDir,'ase_hdr.json')}, ...
        'Mase', {fullfile(sessDir,'ase_base_al_mat.aff12.1D')}, 'MaseInv', {fullfile(sessDir,'ase_base_al_inv_mat.aff12.1D')}, ...
        'A_aparc_ase', {fullfile(sessDir,'anat_fs_aparc_aseg_al_ase.nii.gz')}, 'A_wmparc_ase', {fullfile(sessDir,'anat_fs_wmparc_al_ase.nii.gz')}, ...
        'A_bstem', {fullfile(sessDir,'anat_fs_bstem_al.nii.gz')}, 'A_xvent', {fullfile(sessDir,'anat_fs_xvent_mask_al.nii.gz')}, ...
        'shft', {fullfile(afniDir,'anat_fs_shft.1D')}, 'Mal', {fullfile(afniDir,'anat_fs_brain_al_mat.aff12.1D')}, ...
        'A_aparc_sh', {fullfile(afniDir,'anat_fs_aparc_aseg.nii')}, 'A_aparc_alH', {fullfile(afniDir,'anat_fs_aparc_aseg_al+orig.HEAD')}, ...
        'A_aparc_alB', {fullfile(afniDir,'anat_fs_aparc_aseg_al+orig.BRIK')}, 'A_T1', {fullfile(afniDir,'anat_fs_orig.nii')}, ...
        'fs_T1', {fullfile(fsDir,'T1.nii.gz')}, 'fs_aseg', {fullfile(fsDir,'aseg.nii.gz')}, ...
        'fs_aparc', {fullfile(fsDir,'aparc+aseg.nii.gz')}, 'fs_wmparc', {fullfile(fsDir,'wmparc.nii.gz')}, ...
        'fs_bstem', {fullfile(fsDir,'brainstemSsLabels.v13.FSvoxelSpace.nii.gz')});
    fn = fieldnames(L);
    for k = 1:numel(fn)
        src = L.(fn{k});
        [~, b, e] = fileparts(src);  ext = e;  if strcmp(e, '.gz'), [~, ~, e0] = fileparts(b); ext = [e0 e]; end
        dst = ['in/' fn{k} ext];
        if strcmp(fn{k}, 'A_aparc_alB'), dst = 'in/A_aparc_al+orig.BRIK'; end
        if strcmp(fn{k}, 'A_aparc_alH'), dst = 'in/A_aparc_al+orig.HEAD'; end
        if isfile(src), system(sprintf('ln -sfn "%s" "%s"', src, P_(dst))); I.(fn{k}) = dst; else, I.(fn{k}) = ''; end
    end
    req = {'ase_base','e1','e2','hdr','Mase','MaseInv','A_aparc_ase','A_wmparc_ase','A_xvent','shft','Mal', ...
           'fs_T1','fs_aseg','fs_aparc','fs_wmparc','fs_bstem'};
    miss = req(cellfun(@(f) isempty(I.(f)), req));
    if ~isempty(miss), error('Missing inputs: %s', strjoin(miss, ', ')); end

    [T1, iT1] = rd(P_(I.fs_T1));  aseg = rd(P_(I.fs_aseg));  aparc = rd(P_(I.fs_aparc));  bst = rd(P_(I.fs_bstem));
    AT1 = niftiAffine(iT1.raw);  brain = aseg > 0;
    [aseB, iAse] = rd(P_(I.ase_base));  Aase = niftiAffine(iAse.raw);  voxAse = sqrt(sum(Aase(1:3,1:3).^2,1));
    LUT = readLUT(lutFile);

    %% ---------- 0. provenance: route A rebuilt from FreeSurfer ----------
    note('== 0. Provenance ==');
    if ~isempty(I.A_T1)
        S.prov.T1identical = isequal(T1, rd(P_(I.A_T1)));
        note('FreeSurfer T1 (mri_convert) = session_afni_nortc/anat_fs_orig.nii: %d', S.prov.T1identical);
    end
    afni('3dZeropad -S 20 -I -20 -prefix zp_aparc.nii.gz %s', I.fs_aparc);
    afni('3dZeropad -S 20 -I -20 -prefix zp_bstem.nii.gz %s', I.fs_bstem);
    if ~isempty(I.A_aparc_sh)
        for t = {'aparc','bstem'}
            afni('3dAllineate -1Dmatrix_apply %s -interp NN -final NN -master %s -prefix sh_%s.nii.gz -input zp_%s.nii.gz', I.shft, I.A_aparc_sh, t{1}, t{1});
            afni('3dAllineate -1Dmatrix_apply %s -interp NN -final NN -master %s -prefix al_%s.nii.gz -input sh_%s.nii.gz', I.Mal, I.A_aparc_sh, t{1}, t{1});
        end
        S.prov.shiftDiff = nnz(rd(P_('sh_aparc.nii.gz')) ~= rd(P_(I.A_aparc_sh)));
        note('replayed zero-pad + centre shift vs anat_fs_aparc_aseg.nii: %d voxels differ', S.prov.shiftDiff);
        if ~isempty(I.A_aparc_alH) && ~isempty(I.A_aparc_alB)
            afni('3dcalc -a in/A_aparc_al+orig -expr a -prefix A_aparc_al_stored.nii.gz');
            st = rd(P_('A_aparc_al_stored.nii.gz'));  rp = rd(P_('al_aparc.nii.gz'));
            S.prov.alDiff = nnz(st ~= rp);  S.prov.alN = nnz(st);
            note('replayed full chain vs stored anat_fs_aparc_aseg_al: %d of %d labelled voxels differ', S.prov.alDiff, S.prov.alN);
        end
        S.prov.bstemFS = nnz(bst);  S.prov.bstemAl = nnz(rd(P_('al_bstem.nii.gz')));
        note('brainstem labels: %d voxels in FreeSurfer space, %d after the route-A chain (%+.1f%%)', S.prov.bstemFS, S.prov.bstemAl, ...
            100*(S.prov.bstemAl/S.prov.bstemFS - 1));
        S.prov.aparcCropLoss = nnz(aparc) - nnz(rd(P_('zp_aparc.nii.gz')));
        note('labelled voxels lost to the 20 mm inferior crop: %d', S.prov.aparcCropLoss);
    end
    Mal = readAff(P_(I.Mal));
    S.prov.MalDet = det(Mal(:,1:3));  S.prov.MalSV = svd(Mal(:,1:3))';
    note('anatomy -> BOLD affine (anat_fs_brain_al_mat): det %.4f, singular values [%.4f %.4f %.4f]', S.prov.MalDet, S.prov.MalSV);

    %% ---------- 1. registrations ----------
    note('== 1. Registrations ==');
    afni('3dcalc -a %s -b %s -expr "a*step(b)" -datum float -prefix fs_brain.nii.gz', I.fs_T1, I.fs_aseg);
    afni('3dUnifize -prefix fs_brain_uni.nii.gz fs_brain.nii.gz');
    warps = containers.Map({6, 9, 12}, {'shift_rotate', 'shift_rotate_scale', 'affine_general'});
    routes = {'A'};
    % AFNI composes "apply S, then M" as the single matrix M*S (checked: one-step Mal*shft reproduces the
    % sequential replay of section 0). FreeSurfer -> _al is Mal*shft, _al -> ASE is MaseInv, so the
    % matrix that brings ASE data onto the FreeSurfer grid is inv(MaseInv*Mal*shft) = inv(shft)*inv(Mal)*Mase.
    afni('cat_matvec -ONELINE %s -I %s -I %s > A.aff12.1D', I.shft, I.Mal, I.Mase);
    for d = opts.dofs
        r = sprintf('B%d', d);  routes{end+1} = r; %#ok<AGROW>
        afni(['3dAllineate -base fs_brain_uni.nii.gz -master fs_brain_uni.nii.gz -source %s -prefix %s_reg.nii.gz ' ...
              '%s -warp %s -1Dmatrix_save %s.aff12.1D -1Dparam_save %s.param.1D'], I.ase_base, r, opts.allineateOpts, warps(d), r, r);
    end
    mainB = sprintf('B%d', opts.mainDof);
    for k = 1:numel(routes)
        r = routes{k};
        afni('cat_matvec -ONELINE %s.aff12.1D -I > %s_inv.aff12.1D', r, r);
        afni('3dAllineate -1Dmatrix_apply %s.aff12.1D -source %s -master fs_brain_uni.nii.gz -final wsinc5 -prefix base_%s_FS.nii.gz', r, I.ase_base, r);
        M = readAff(P_([r '.aff12.1D']));  Lm = M(:,1:3);
        par = struct('det', det(Lm), 'sv', svd(Lm)');
        if startsWith(r, 'B'), p = readParam(P_([r '.param.1D'])); par.shift = p(1:3); par.angle = p(4:6); par.scale = p(7:9); par.shear = p(10:12); end
        cst = afniOut(W, sprintf('3dAllineate -base fs_brain_uni.nii.gz -source base_%s_FS.nii.gz -autoweight -nocmass -allcostX', r));
        par.cost = parseCosts(cst);
        S.reg.(r) = par;
        note('%s: det %.4f, singular values [%.4f %.4f %.4f]; cost vs FreeSurfer target: lpa+ %.4f, nmi %.4f, lpc %.4f', r, par.det, par.sv, ...
            getc(par.cost,'lpa+'), getc(par.cost,'nmi'), getc(par.cost,'lpc'));
        if isfield(par, 'scale')
            note('    scales [%.4f %.4f %.4f], shears [%.4f %.4f %.4f], angles [%.2f %.2f %.2f] deg', par.scale, par.shear, par.angle);
        end
    end

    %% ---------- 2. displacement between routes ----------
    note('== 2. Where the routes place each brain voxel in the raw ASE data ==');
    for e = {'i','j','k','1'}
        afni('3dcalc -a %s -expr %s -datum float -prefix idx_%s.nii.gz', I.ase_base, e{1}, e{1});
    end
    afni('3dTcat -prefix idx.nii.gz idx_i.nii.gz idx_j.nii.gz idx_k.nii.gz idx_1.nii.gz');
    X = struct();  cover = true(size(brain));
    for k = 1:numel(routes)
        r = routes{k};
        afni('3dAllineate -1Dmatrix_apply %s.aff12.1D -source idx.nii.gz -master fs_brain_uni.nii.gz -final linear -prefix idx_%s.nii.gz', r, r);
        x = single(rd(P_(['idx_' r '.nii.gz'])));  X.(r) = x(:,:,:,1:3);  cover = cover & x(:,:,:,4) > 0.999;
        S.cover.(r) = x(:,:,:,4) > 0.999;
    end
    valid = brain & cover;
    S.disp.brainCoveredFrac = nnz(valid)/nnz(brain);
    note('brain voxels inside the ASE slab for every route: %d (%.1f%% of the FreeSurfer brain)', nnz(valid), 100*S.disp.brainCoveredFrac);
    others = setdiff(arrayfun(@(d) sprintf('B%d', d), opts.dofs, 'UniformOutput', false), mainB, 'stable');
    pairs = {'A', mainB};
    for k = 1:numel(others), pairs(end+1,:) = {'A', others{k}}; end %#ok<AGROW>
    for k = 1:numel(others), pairs(end+1,:) = {mainB, others{k}}; end %#ok<AGROW>
    D = struct('pair', {}, 'median', {}, 'p95', {}, 'max', {}, 'meanRAS', {}, 'medAbsRAS', {});
    dmapMain = [];  dmAll = struct();
    for p = 1:size(pairs,1)
        a = X.(pairs{p,1});  b = X.(pairs{p,2});
        di = reshape(a - b, [], 3);  di = double(di(valid(:),:));
        dW = di * Aase(1:3,1:3).';                                   % world (RAS) mm
        dm = sqrt(sum(dW.^2, 2));
        D(end+1) = struct('pair', sprintf('%s vs %s', pairs{p,:}), 'median', median(dm), 'p95', prctile(dm,95), ...
            'max', max(dm), 'meanRAS', mean(dW,1), 'medAbsRAS', median(abs(dW),1)); %#ok<AGROW>
        note('%-9s: |d| median %.2f mm, 95th %.2f, max %.2f; mean signed (R,A,S) [%.2f %.2f %.2f]; median |d| per axis [%.2f %.2f %.2f]', ...
            D(end).pair, D(end).median, D(end).p95, D(end).max, D(end).meanRAS, D(end).medAbsRAS);
        dmAll.(sprintf('%s_%s', pairs{p,:})) = zeros(size(brain), 'single');  dmAll.(sprintf('%s_%s', pairs{p,:}))(valid) = dm;
        if p == 1, dmapMain = zeros(size(brain), 'single'); dmapMain(valid) = dm; dRASmain = zeros([size(brain) 3], 'single'); for c = 1:3, t = zeros(size(brain),'single'); t(valid) = dW(:,c); dRASmain(:,:,:,c) = t; end, end
    end
    S.disp.pairs = struct2table(D);
    % per region (A vs main B)
    codes = intersect(unique(aparc(valid)), LUT.code);
    R = table('Size', [numel(codes) 5], 'VariableTypes', {'double','string','double','double','double'}, ...
        'VariableNames', {'code','name','nVox1mm','medianDisp','p95Disp'});
    for k = 1:numel(codes)
        m = valid & aparc == codes(k);
        R(k,:) = {codes(k), LUT.name(LUT.code == codes(k)), nnz(m), median(dmapMain(m)), prctile(dmapMain(m), 95)};
    end
    S.disp.byRegion = sortrows(R, 'medianDisp', 'descend');
    note('regions with the largest A vs %s displacement (median mm): %s', mainB, strjoin(compose('%s %.1f', ...
        S.disp.byRegion.name(1:min(8,end)), S.disp.byRegion.medianDisp(1:min(8,end))), ', '));
    figDisp(dmapMain, dRASmain, valid, T1, AT1, S.disp.pairs, mainB, C, seqMap, fullfile(outDir, 'anat1_displacement.png'));
    % overlay slices: where A and the least similar B variant shown (the last column) differ most
    lastB = sprintf('B%d', min(opts.dofs));  dOv = dmapMain;
    if isfield(dmAll, ['A_' lastB]), dOv = dmAll.(['A_' lastB]); end
    figOverlay(routes, mainB, opts.dofs, T1, aseg, dOv, valid, AT1, W, C, fullfile(outDir, 'anat2_overlay.png'), lastB);

    %% ---------- 3. ROI membership on the ASE grid ----------
    note('== 3. ROI membership on the ASE grid (occupancy >= 0.25, as the pipeline) ==');
    seg = struct();
    % route A, as the pipeline wrote it (brainstem made the same way as stage 7)
    afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final NN -prefix A_bstem_ase.nii.gz', I.MaseInv, I.A_bstem, I.A_bstem);
    seg.Apipe = {P_(I.A_aparc_ase), P_(I.A_wmparc_ase), P_('A_bstem_ase.nii.gz')};
    % every route rebuilt from FreeSurfer with one NN step (A: the composite of its stored matrices)
    labs = {'fs_aparc','fs_wmparc','fs_bstem'};
    for k = 1:numel(routes)
        r = routes{k};  seg.(r) = {};
        for l = 1:numel(labs)
            out = sprintf('%s_%s_ase.nii.gz', r, labs{l});
            afni('3dAllineate -1Dmatrix_apply %s_inv.aff12.1D -source %s -master %s -final NN -prefix %s', r, I.(labs{l}), I.(labs{l}), out);
            seg.(r){end+1} = P_(out);
        end
    end
    % reference: main B with the labels moved 1 mm along anterior-posterior (the EPI phase-encode axis)
    [~, apAx] = max(abs(AT1(2,1:3)));
    seg.ref = {};
    for l = 1:numel(labs)
        f = seg.(mainB){l};  [v, inf_] = rd(f);  v = circshift(v, 1, apAx);
        out = P_(sprintf('ref_%s_ase.nii', labs{l}));  niftiwrite(v, out, inf_);  seg.ref{end+1} = out;
    end
    rlist = [{'Apipe'}, routes, {'ref'}];
    Mm = struct();
    for k = 1:numel(rlist)
        [sv_, lut_] = combineSegVols(seg.(rlist{k}), {lutFile}, 'keepOverlaps', 1, 'warnSkipped', false);
        if k == 1, LUTm = lut_; end
        assert(isequal(lut_.code, LUTm.code), 'LUTs differ between routes');
        Mm.(rlist{k}) = segToAseMembership(sv_, lut_, niftiinfo(seg.(rlist{k}){1}), iAse, 'maxSampling', [5 5 5], 'minFrac', 0.25);
    end
    dice = @(P, Q) 2*full(sum(min(P,Q),2)) ./ full(sum(P,2) + sum(Q,2));
    RO = table(LUTm.code, LUTm.name, full(sum(Mm.Apipe,2)), full(sum(Mm.(mainB),2)), ...
        dice(Mm.Apipe, Mm.A), dice(Mm.Apipe, Mm.(mainB)), dice(Mm.A, Mm.(mainB)), dice(Mm.(mainB), Mm.ref), ...
        'VariableNames', {'code','name','nVox_A','nVox_B','dice_Apipe_Acomp','dice_A_B','dice_Acomp_B','dice_B_ref1mm'});
    for d = setdiff(opts.dofs, opts.mainDof)
        RO.(sprintf('dice_A_B%d', d)) = dice(Mm.Apipe, Mm.(sprintf('B%d', d)));
        RO.(sprintf('dice_B_B%d', d)) = dice(Mm.(mainB), Mm.(sprintf('B%d', d)));
    end
    S.roi = RO;  ok = RO.nVox_A > 0 & RO.nVox_B > 0;
    note('check: route A rebuilt from its stored matrices vs the pipeline''s own ASE-space labels: median Dice %.4f (min %.4f)', ...
        median(RO.dice_Apipe_Acomp(ok)), min(RO.dice_Apipe_Acomp(ok)));
    note('A (as the pipeline wrote it) vs %s: median Dice %.3f (IQR %.3f-%.3f); %s vs %s moved 1 mm A-P: median %.3f; A-B agreement is worse than a 1 mm shift in %.0f%% of %d regions', ...
        mainB, median(RO.dice_A_B(ok)), prctile(RO.dice_A_B(ok),25), prctile(RO.dice_A_B(ok),75), mainB, mainB, median(RO.dice_B_ref1mm(ok)), ...
        100*mean(RO.dice_A_B(ok) < RO.dice_B_ref1mm(ok)), nnz(ok));
    note('transform only (route A rebuilt in one step from FreeSurfer) vs %s: median Dice %.3f (IQR %.3f-%.3f)', mainB, ...
        median(RO.dice_Acomp_B(ok)), prctile(RO.dice_Acomp_B(ok),25), prctile(RO.dice_Acomp_B(ok),75));
    for d = setdiff(opts.dofs, opts.mainDof)
        note('%s vs B%d: median Dice %.3f; A vs B%d: %.3f', mainB, d, median(RO.(sprintf('dice_B_B%d', d))(ok)), d, median(RO.(sprintf('dice_A_B%d', d))(ok)));
    end
    figDice(RO, mainB, C, fullfile(outDir, 'anat3_roi_dice.png'));

    %% ---------- 4. CSF fraction (quality bit 4) ----------
    note('== 4. Tissue fraction on the ASE grid (CSF = 1 - xvent) ==');
    afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final linear -prefix A_xvent_ase.nii.gz', I.MaseInv, I.A_xvent, I.ase_base);
    afni('3dcalc -a %s -expr "step(a)-step(equals(a,4)+equals(a,5)+equals(a,43)+equals(a,44)+equals(a,14)+equals(a,15))" -prefix fs_xvent.nii.gz', I.fs_aparc);
    afni('3dAllineate -1Dmatrix_apply %s_inv.aff12.1D -source fs_xvent.nii.gz -master %s -final linear -prefix B_xvent_ase.nii.gz', mainB, I.ase_base);
    afni('3dAllineate -1Dmatrix_apply A_inv.aff12.1D -source fs_xvent.nii.gz -master %s -final linear -prefix Acomp_xvent_ase.nii.gz', I.ase_base);
    xa = double(rd(P_('A_xvent_ase.nii.gz')));  xb = double(rd(P_('B_xvent_ase.nii.gz')));  xc = double(rd(P_('Acomp_xvent_ase.nii.gz')));
    m = xa > 0.05 | xb > 0.05;
    S.csf = struct('nVox', nnz(m), 'meanAbsDiff', mean(abs(xa(m)-xb(m))), 'p95AbsDiff', prctile(abs(xa(m)-xb(m)),95), ...
        'flagA', nnz(m & (1-xa) > 0.2), 'flagB', nnz(m & (1-xb) > 0.2), 'flagDisagree', nnz(m & xor((1-xa) > 0.2, (1-xb) > 0.2)));
    note('tissue fraction, A vs %s over %d ASE voxels: mean |diff| %.3f (95th %.3f); CSF > 0.2 flagged: A %d, B %d, flag differs in %d voxels (%.1f%%)', ...
        mainB, S.csf.nVox, S.csf.meanAbsDiff, S.csf.p95AbsDiff, S.csf.flagA, S.csf.flagB, S.csf.flagDisagree, 100*S.csf.flagDisagree/S.csf.nVox);
    S.csf.resamp = struct('meanAbsDiff', mean(abs(xa(m)-xc(m))), 'p95AbsDiff', prctile(abs(xa(m)-xc(m)),95), ...
        'flagDisagree', nnz(m & xor((1-xa) > 0.2, (1-xc) > 0.2)));
    S.csf.transf = struct('meanAbsDiff', mean(abs(xc(m)-xb(m))), 'p95AbsDiff', prctile(abs(xc(m)-xb(m)),95), ...
        'flagDisagree', nnz(m & xor((1-xc) > 0.2, (1-xb) > 0.2)));
    note('  of which resampling only (A as stored vs A rebuilt in one step): mean |diff| %.3f, flag differs in %.1f%%', ...
        S.csf.resamp.meanAbsDiff, 100*S.csf.resamp.flagDisagree/S.csf.nVox);
    note('  and transform only (A rebuilt vs %s): mean |diff| %.3f, flag differs in %.1f%%', mainB, ...
        S.csf.transf.meanAbsDiff, 100*S.csf.transf.flagDisagree/S.csf.nVox);

    %% ---------- 5. brainstem ----------
    note('== 5. Brainstem ==');
    bl = [173 174 175 178];  bn = ["Midbrain","Pons","Medulla","SCP"];
    BS = table(bl', bn', 'VariableNames', {'code','name'});
    for k = 1:numel(bl)
        mk = bst == bl(k);
        BS.nVox1mm(k) = nnz(mk);  BS.inSlabFrac(k) = nnz(mk & S.cover.(mainB))/nnz(mk);
        BS.inAsegBrainStem(k) = nnz(mk & aparc == 16)/nnz(mk);
        j = find(RO.code == bl(k));
        if isempty(j), BS.nVoxASE_A(k) = NaN; BS.nVoxASE_B(k) = NaN; BS.dice_A_B(k) = NaN;
        else, BS.nVoxASE_A(k) = RO.nVox_A(j); BS.nVoxASE_B(k) = RO.nVox_B(j); BS.dice_A_B(k) = RO.dice_A_B(j); end
    end
    S.bstem = BS;
    for k = 1:height(BS)
        note('%-8s: %5d mm^3, %.0f%% inside the ASE slab, %.0f%% inside aseg Brain-Stem; ASE voxels A %.1f / B %.1f; Dice A-B %.3f', ...
            BS.name(k), BS.nVox1mm(k), 100*BS.inSlabFrac(k), 100*BS.inAsegBrainStem(k), BS.nVoxASE_A(k), BS.nVoxASE_B(k), BS.dice_A_B(k));
    end
    figBstem(T1, bst, S.cover.(mainB), AT1, C, fullfile(outDir, 'anat5_brainstem.png'));

    %% ---------- 6. ROI values ----------
    if opts.doFits
        note('== 6. ROI values (stats-only fit, default settings) ==');
        h = jsondecode(fileread(P_(I.hdr)));
        if isfield(h, 'B0_T'), B0 = h.B0_T; else, B0 = h.ImagingFrequency / 42.577478518; end
        note('B0 = %.4f T', B0);
        fitR = [{'Apipe'}, arrayfun(@(d) sprintf('B%d', d), opts.dofs, 'UniformOutput', false), {'ref'}];
        V = struct();
        for k = 1:numel(fitR)
            r = fitR{k};  od = P_(['fit_' r]);  if ~isfolder(od), mkdir(od); end
            csv = fullfile(od, 'ase_ROI_Stats.csv');
            if ~(opts.reuse && isfile(csv))
                ASE_fitting(B0, P_(I.hdr), {P_(I.e1), P_(I.e2)}, 'segFiles', seg.(r), 'segLUTs', {lutFile}, ...
                    'writeOpts', {'stats'}, 'overwriteOpts', {'stats'}, 'outPref', fullfile(od, 'ase'), 'warnSegIssues', false);
            end
            V.(r) = readtable(csv, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
        end
        st = {'OEF','R2P','vCBV','R2'};
        rows = intersect(V.Apipe.Properties.RowNames, V.(mainB).Properties.RowNames, 'stable');
        VT = table(string(rows), 'VariableNames', {'name'});
        VT.nVox_A = V.Apipe{rows, 'nVox'};  VT.nVox_B = V.(mainB){rows, 'nVox'};
        for s = st
            VT.([s{1} '_A']) = V.Apipe{rows, s{1}};
            for k = 2:numel(fitR), VT.([s{1} '_' fitR{k}]) = V.(fitR{k}){rows, s{1}}; end
        end
        S.vals = VT;
        for s = st
            a = VT.([s{1} '_A']);  b = VT.([s{1} '_' mainB]);  rf = VT.([s{1} '_ref']);
            okv = isfinite(a) & isfinite(b) & isfinite(rf) & a ~= 0 & b ~= 0;
            rel = abs(b(okv)-a(okv))./abs(a(okv));  relRef = abs(rf(okv)-b(okv))./abs(b(okv));
            S.valSummary.(s{1}) = struct('medRelAB', median(rel), 'p90RelAB', prctile(rel,90), 'medRelRef', median(relRef), ...
                'p90RelRef', prctile(relRef,90), 'n', nnz(okv), 'medSigned', median((b(okv)-a(okv))./abs(a(okv))));
            note('%-4s: |%s - A|/A median %.1f%% (90th %.1f%%), signed median %+.1f%%; 1 mm A-P shift of %s: median %.1f%% (90th %.1f%%); n = %d ROIs', ...
                s{1}, mainB, 100*median(rel), 100*prctile(rel,90), 100*S.valSummary.(s{1}).medSigned, mainB, 100*median(relRef), 100*prctile(relRef,90), nnz(okv));
            for d = setdiff(opts.dofs, opts.mainDof)
                c = VT.([s{1} '_B' num2str(d)]);  ok2 = isfinite(b) & isfinite(c) & b ~= 0;
                note('      B%d vs %s: median %.1f%%', d, mainB, 100*median(abs(c(ok2)-b(ok2))./abs(b(ok2))));
            end
        end
        nf = @(r) nnz(isfinite(V.(r){rows,'OEF'}) & V.(r){rows,'OEF'} > 0);
        note('ROIs with finite, positive OEF: A %d, %s %d (of %d)', nf('Apipe'), mainB, nf(mainB), numel(rows));
        figVals(VT, mainB, C, fullfile(outDir, 'anat4_roi_values.png'));
    end

    %% ---------- write ----------
    S.log = logLines;
    fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
    fprintf(fid, '%s\n', logLines{:});
    fprintf(fid, '\n== Regions by A vs %s displacement (median, 95th, mm) ==\n', mainB);
    for k = 1:height(S.disp.byRegion)
        fprintf(fid, '%6d %-40s %6d %6.2f %6.2f\n', S.disp.byRegion.code(k), S.disp.byRegion.name(k), S.disp.byRegion.nVox1mm(k), ...
            S.disp.byRegion.medianDisp(k), S.disp.byRegion.p95Disp(k));
    end
    fclose(fid);
    writetable(S.roi, fullfile(outDir, 'roi_membership.csv'));
    if isfield(S, 'vals'), writetable(S.vals, fullfile(outDir, 'roi_values.csv')); end
    save(fullfile(outDir, 'summary.mat'), 'S');
    fprintf('Wrote figures and summary to %s\n', outDir);
end

%% ===================== helpers =====================
function [v, info] = rd(f)
    v = niftiread(f);  info = niftiinfo(f);
end

function afni(fmt, varargin)
    persistent reuse wd
    if strcmp(fmt, 'reuse'), reuse = varargin{1}; return; end
    if strcmp(fmt, 'dir'), wd = varargin{1}; return; end
    cmd = sprintf(fmt, varargin{:});
    tok = regexp(cmd, '-prefix (\S+)', 'tokens', 'once');
    if ~isempty(reuse) && reuse && ~isempty(tok) && isfile(fullfile(wd, tok{1})), return; end
    if ~isempty(tok), cmd = [cmd ' -overwrite']; end
    if contains(cmd, '> ') && ~isempty(reuse) && reuse
        t = regexp(cmd, '> (\S+)', 'tokens', 'once');  if isfile(fullfile(wd, t{1})), return; end
    end
    [st, out] = system(sprintf('export PATH="$HOME/abin:$PATH"; export AFNI_NIFTI_TYPE_WARN=NO; cd "%s" && %s', wd, cmd));
    if st ~= 0, error('AFNI command failed:\n%s\n%s', cmd, out); end
end

function out = afniOut(wd, cmd)
    [~, out] = system(['export PATH="$HOME/abin:$PATH"; export AFNI_NIFTI_TYPE_WARN=NO; cd "' wd '" && ' cmd ' 2>&1']);
end

function c = parseCosts(txt)
    c = struct();
    tk = regexp(txt, '\n\s*([A-Za-z+]+)\s*=\s*([-+0-9.eE]+)', 'tokens');
    for k = 1:numel(tk)
        nm = matlab.lang.makeValidName(strrep(tk{k}{1}, '+', 'P'));
        if ~isfield(c, nm), c.(nm) = str2double(tk{k}{2}); end
    end
end

function v = getc(c, name)
    nm = matlab.lang.makeValidName(strrep(name, '+', 'P'));
    if isfield(c, nm), v = c.(nm); else, v = NaN; end
end

function M = readAff(f)
    t = fileread(f);  t = regexprep(t, '#[^\n]*', '');
    v = sscanf(t, '%f');  M = reshape(v(1:12), 4, 3)';
end

function p = readParam(f)
% 3dAllineate -1Dparam_save: shifts, angles, then (if free) scales and shears; fixed ones are 1 and 0
    t = fileread(f);  t = regexprep(t, '#[^\n]*', '');  v = sscanf(t, '%f')';
    p = [0 0 0 0 0 0 1 1 1 0 0 0];  p(1:min(12,numel(v))) = v(1:min(12,numel(v)));
end

function T = readLUT(f)
    T = readtable(f, 'ReadVariableNames', false, 'TextType', 'string', 'FileType', 'text', 'Delimiter', '\t');
    T.Properties.VariableNames = {'code','name'};
    if ~isnumeric(T.code), T(1,:) = []; T.code = str2double(T.code); end
end

function [V, vs] = toRAS(V, A)
    R = A(1:3,1:3);  perm = zeros(1,3);  sg = zeros(1,3);
    for w = 1:3, [~, perm(w)] = max(abs(R(w,:))); sg(w) = sign(R(w,perm(w))); end
    V = permute(V, [perm 4:max(4,ndims(V))]);
    for w = 1:3, if sg(w) < 0, V = flip(V, w); end, end
    vs = sqrt(sum(R(:,perm).^2, 1));
end

function im = sl(V, view, k)
    switch view
        case 'ax',  im = flipud(V(:,:,k)');
        case 'cor', im = flipud(squeeze(V(:,k,:))');
        case 'sag', im = flipud(squeeze(V(k,:,:))');
    end
end

function styleAx(ax, C)
    set(ax, 'XColor', C.ink2, 'YColor', C.ink2, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.75);
end

function sv(f, file)
    exportgraphics(f, file, 'Resolution', 200, 'BackgroundColor', 'white');  close(f);
end

function [rr, cc] = cropBox(m, pad)
    [r, c] = find(m);
    if isempty(r), rr = 1:size(m,1); cc = 1:size(m,2); return; end
    rr = max(1, min(r)-pad):min(size(m,1), max(r)+pad);  cc = max(1, min(c)-pad):min(size(m,2), max(c)+pad);
end

function k = worstSlices(Dm, Vd)
% slices (ax = dim 3, cor = dim 2, sag = dim 1 of an RAS volume) with the largest mean displacement,
% among slices holding at least 30% of the largest slice's brain
    k = zeros(1,3);  dimOf = [3 2 1];
    for a = 1:3
        d = dimOf(a);  o = setdiff(1:3, d);
        s = sum(sum(double(Dm) .* Vd, o(1)), o(2));  n = sum(sum(double(Vd), o(1)), o(2));
        s = s(:);  n = n(:);  mval = s ./ max(n,1);  mval(n < 0.3*max(n)) = -Inf;
        [~, k(a)] = max(mval);
    end
end

%% ===================== figures =====================
function figDisp(dMap, dRAS, valid, T1, A, P, mainB, C, seqMap, file)
    D = toRAS(dMap, A);  B = toRAS(single(T1), A);  Vd = toRAS(valid, A);  Dap = toRAS(dRAS(:,:,:,2), A);
    k = worstSlices(D, Vd);  views = {'ax','cor','sag'};  vt = {'Axial','Coronal','Sagittal'};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1900 1150]);
    tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    cmax = max(2, ceil(prctile(dMap(valid), 99)));
    bl = prctile(B(B>0), 99);
    divMap = interp1([0 0.5 1], [235 104 52; 255 255 255; 42 120 214]/255, linspace(0,1,256));
    for row = 1:2
        for p = 1:3
            ax = nexttile(tl, (row-1)*4 + p);
            a = sl(B, views{p}, k(p)) / bl;  m = sl(Vd, views{p}, k(p));
            if row == 1, d = sl(D, views{p}, k(p)); idx = max(1, min(256, round(d/cmax*255)+1)); cm = seqMap;
            else, d = sl(Dap, views{p}, k(p)); idx = max(1, min(256, round((d/cmax+1)/2*255)+1)); cm = divMap; end
            rgbA = repmat(min(a,1), 1, 1, 3);  rgbD = reshape(cm(idx,:), [size(d) 3]);
            rgb = rgbA;  mm = repmat(m, 1, 1, 3);  rgb(mm) = rgbD(mm);
            [rr, cc] = cropBox(sl(B, views{p}, k(p)) > 0, 6);
            image(ax, rgb(rr,cc,:));  axis(ax, 'image', 'off');
            if row == 1, title(ax, vt{p}, 'FontWeight', 'normal', 'Color', C.ink); end
        end
        colormap(ax, cm);
        if row == 1, clim(ax, [0 cmax]); lbl = sprintf('distance A vs %s (mm)', mainB);
        else, clim(ax, [-cmax cmax]); lbl = sprintf('A − %s along anterior (+) / posterior (mm)', mainB); end
        cb = colorbar(ax, 'eastoutside');  cb.Label.String = lbl;  cb.Color = C.ink2;
    end
    ax = nexttile(tl, 4);  hold(ax, 'on');
    edges = 0:0.1:max(4, ceil(max(dMap(valid))));
    histogram(ax, dMap(valid), edges, 'FaceColor', C.B, 'EdgeColor', 'w', 'FaceAlpha', 1, 'Normalization', 'probability');
    yl = ylim(ax);
    for v = [1.72 4.8], plot(ax, [v v], yl, '--', 'Color', C.ink2, 'LineWidth', 1); end
    text(ax, 1.72, yl(2)*0.95, ' in-plane voxel 1.72 mm', 'Color', C.ink2, 'FontSize', 10);
    text(ax, 4.8, yl(2)*0.85, ' slice 4.8 mm', 'Color', C.ink2, 'FontSize', 10);
    xlabel(ax, 'distance (mm)');  ylabel(ax, 'fraction of brain voxels');  styleAx(ax, C);  xlim(ax, [0 edges(end)]);
    title(ax, sprintf('A vs %s: median %.2f · 95th %.2f · max %.2f mm', mainB, P.median(1), P.p95(1), P.max(1)), 'FontWeight', 'normal', 'Color', C.ink);
    ax = nexttile(tl, 8);  axis(ax, 'off');
    txt = compose('%-10s median %.2f   95th %.2f   max %.2f mm', string(P.pair), P.median, P.p95, P.max);
    text(ax, 0, 0.9, ['Pairwise distances (brain voxels in the slab)'; txt], 'VerticalAlignment', 'top', 'FontName', 'Menlo', 'FontSize', 11, 'Color', C.ink);
    title(tl, 'Where each route places the same brain voxel in the raw ASE data', 'FontSize', 14, 'Color', C.ink);
    subtitle(tl, sprintf(['Top: distance between route A (current, via the BOLD-aligned anatomy) and route %s (FreeSurfer directly), shown on the FreeSurfer T1 at the slices where they differ most. ' ...
        'Bottom: the anterior–posterior component (the EPI phase-encode axis).'], mainB), 'Color', C.ink2, 'FontSize', 11);
    sv(f, file);
end

function figOverlay(routes, mainB, dofs, T1, aseg, dMap, valid, A, W, C, file, lastB)
    show = unique([{'A'}, {mainB}, {lastB}], 'stable');
    lab = struct('A', 'Route A (current)');
    for k = 1:numel(routes), if startsWith(routes{k}, 'B'), lab.(routes{k}) = sprintf('Route B, %s parameters', routes{k}(2:end)); end, end
    B = toRAS(single(T1), A);  Wm = toRAS(ismember(aseg, [2 41]), A);  Ve = toRAS(ismember(aseg, [4 43]), A);
    Ctx = toRAS(ismember(aseg, [3 42]), A);  Vd = toRAS(valid, A);  Dm = toRAS(dMap, A);
    k = worstSlices(Dm, Vd);  views = {'ax','cor','sag'};  vt = {'Axial','Coronal','Sagittal'};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1500 1450]);
    tl = tiledlayout(f, 3, numel(show), 'TileSpacing', 'compact', 'Padding', 'compact');
    imgs = cell(1, numel(show));
    for c = 1:numel(show), imgs{c} = toRAS(single(niftiread(fullfile(W, ['base_' show{c} '_FS.nii.gz']))), A); end
    for p = 1:3
        [rr, cc] = cropBox(sl(B, views{p}, k(p)) > 0, 6);
        for c = 1:numel(show)
            ax = nexttile(tl);
            im = sl(imgs{c}, views{p}, k(p));  sc = prctile(im(sl(Vd, views{p}, k(p))), 99);
            imagesc(ax, min(im(rr,cc)/sc, 1), [0 1]);  colormap(ax, gray);  axis(ax, 'image', 'off');  hold(ax, 'on');
            w = double(sl(Wm, views{p}, k(p)));  v = double(sl(Ve, views{p}, k(p)));  g = double(sl(Ctx, views{p}, k(p)) | sl(Wm, views{p}, k(p)));
            cts = {w, C.B, 0.9; v, C.A, 0.9; g, [0.95 0.85 0.2], 0.6};
            for t = 1:3
                z = cts{t,1}(rr,cc);
                if any(z(:)) && ~all(z(:)), contour(ax, z, [0.5 0.5], 'EdgeColor', cts{t,2}, 'LineWidth', cts{t,3}); end
            end
            if p == 1, title(ax, lab.(show{c}), 'FontWeight', 'normal', 'Color', C.ink, 'FontSize', 13); end
            if c == 1, text(ax, -8, numel(rr)/2, vt{p}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Color', C.ink2, 'FontSize', 12); end
        end
    end
    title(tl, 'ASE base image (spin echo, frame 0) placed on the FreeSurfer anatomy by each route', 'FontSize', 14, 'Color', C.ink);
    subtitle(tl, sprintf('FreeSurfer contours: blue = white matter, yellow = pial (cortex + WM), orange = lateral ventricles. Slices where routes A and %s differ most.', lastB), ...
        'Color', C.ink2, 'FontSize', 12);
    sv(f, file);
end

function figDice(T, mainB, C, file)
    cls = repmat("Subcortical / other", height(T), 1);
    cls(T.code >= 1000 & T.code < 3000) = "Cortex";
    cls(T.code >= 3000 & T.code < 5003 | ismember(T.code, [2 41])) = "White matter";
    cls(ismember(T.code, [16 173 174 175 178])) = "Brainstem";
    groups = ["Cortex","White matter","Subcortical / other","Brainstem"];
    col = {[42 120 214]/255, [235 104 52]/255, [98 80 214]/255, [20 150 90]/255};
    ok = T.nVox_A > 0 & T.nVox_B > 0;
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1500 620]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(tl);  hold(ax, 'on');
    for g = 1:numel(groups)
        m = ok & cls == groups(g);
        scatter(ax, T.nVox_B(m), T.dice_A_B(m), 36, col{g}, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'DisplayName', groups(g));
    end
    set(ax, 'XScale', 'log');  ylim(ax, [min([T.dice_A_B(ok); T.dice_B_ref1mm(ok)]) - 0.02 1.005]);
    xlabel(ax, 'region size on the ASE grid (effective voxels)');  ylabel(ax, sprintf('weighted Dice, A vs %s', mainB));
    legend(ax, 'Location', 'southeast', 'Box', 'off');  styleAx(ax, C);
    title(ax, sprintf('Route A vs %s: median Dice %.3f', mainB, median(T.dice_A_B(ok))), 'FontWeight', 'normal', 'Color', C.ink);
    ax = nexttile(tl);  hold(ax, 'on');
    for g = 1:numel(groups)
        m = ok & cls == groups(g);
        scatter(ax, T.dice_B_ref1mm(m), T.dice_A_B(m), 30, col{g}, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
    end
    lim = [min([T.dice_A_B(ok); T.dice_B_ref1mm(ok)]) - 0.02 1.005];
    plot(ax, lim, lim, '--', 'Color', C.ink2, 'LineWidth', 1);
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);
    xlabel(ax, sprintf('Dice, %s vs %s with labels moved 1 mm anterior–posterior', mainB, mainB));
    ylabel(ax, sprintf('Dice, A vs %s', mainB));  styleAx(ax, C);
    title(ax, sprintf('Below the line: routes differ by more than a 1 mm shift (%.0f%% of regions)', 100*mean(T.dice_A_B(ok) < T.dice_B_ref1mm(ok))), ...
        'FontWeight', 'normal', 'Color', C.ink);
    title(tl, 'ROI membership on the ASE grid (pipeline method, occupancy ≥ 0.25)', 'FontSize', 13, 'Color', C.ink);
    sv(f, file);
end

function figVals(VT, mainB, C, file)
    st = {'OEF','R2P','vCBV'};  nm = {'OEF','R2′ (1/s)','vCBV'};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1800 560]);
    tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:3
        ax = nexttile(tl);  hold(ax, 'on');
        a = VT.([st{s} '_A']);  b = VT.([st{s} '_' mainB]);  r = VT.([st{s} '_ref']);  n = VT.nVox_B;
        ok = isfinite(a) & isfinite(b) & isfinite(r) & a ~= 0 & b ~= 0;
        scatter(ax, n(ok), 100*(r(ok)-b(ok))./abs(b(ok)), 18, C.ref, 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', sprintf('%s moved 1 mm A–P vs %s', mainB, mainB));
        scatter(ax, n(ok), 100*(b(ok)-a(ok))./abs(a(ok)), 26, C.B, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.4, 'DisplayName', sprintf('%s vs A', mainB));
        set(ax, 'XScale', 'log');  yline(ax, 0, '-', 'Color', C.ink2, 'HandleVisibility', 'off');
        lim = prctile(abs(100*(b(ok)-a(ok))./abs(a(ok))), 98);  ylim(ax, [-1 1]*max(5, 1.1*lim));
        xlabel(ax, 'region size on the ASE grid (effective voxels)');  ylabel(ax, ['change in ROI ' nm{s} ' (%)']);  styleAx(ax, C);
        title(ax, nm{s}, 'FontWeight', 'normal', 'Color', C.ink);
        if s == 1, legend(ax, 'Location', 'southeast', 'Box', 'off'); end
    end
    title(tl, sprintf('ROI values (fit of the region-averaged signal): route %s vs route A, against the effect of a 1 mm shift', mainB), 'FontSize', 13, 'Color', C.ink);
    sv(f, file);
end

function figBstem(T1, bst, cover, A, C, file)
    B = toRAS(single(T1), A);  L = toRAS(bst, A);  Cv = toRAS(cover, A);
    [x, y, ~] = ind2sub(size(L), find(L > 0));  kx = round(median(x));  ky = round(median(y));
    col = [235 104 52; 42 120 214; 20 150 90; 98 80 214]/255;  codes = [173 174 175 178];
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1300 640]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    vv = {'sag', kx; 'cor', ky};  vt = {'Sagittal (midline)','Coronal (through the brainstem)'};
    for p = 1:2
        ax = nexttile(tl);
        a = sl(B, vv{p,1}, vv{p,2});  a = min(a / prctile(B(B>0), 99), 1);  rgb = repmat(a, 1, 1, 3);
        l = sl(L, vv{p,1}, vv{p,2});
        for c = 1:4
            m = l == codes(c);
            for ch = 1:3, t = rgb(:,:,ch); t(m) = 0.35*t(m) + 0.65*col(c,ch); rgb(:,:,ch) = t; end
        end
        [rr, cc] = cropBox(a > 0.05, 4);
        image(ax, rgb(rr,cc,:));  axis(ax, 'image', 'off');  hold(ax, 'on');
        cv = double(sl(Cv, vv{p,1}, vv{p,2}));
        contour(ax, cv(rr,cc), [0.5 0.5], 'EdgeColor', [0.95 0.85 0.2], 'LineWidth', 1.5);
        title(ax, vt{p}, 'FontWeight', 'normal', 'Color', C.ink);
    end
    lg = arrayfun(@(c) patch(ax, NaN, NaN, col(c,:), 'EdgeColor', 'none'), 1:4);
    legend(ax, lg, {'Midbrain','Pons','Medulla','SCP'}, 'Location', 'eastoutside', 'Box', 'off', 'TextColor', C.ink, 'FontSize', 11);
    title(tl, 'FreeSurfer brainstem labels (brainstemSsLabels.v13) and the ASE slab', 'FontSize', 13, 'Color', C.ink);
    subtitle(tl, 'Yellow: edge of the ASE field of view (route B placement)', 'Color', C.ink2, 'FontSize', 11);
    sv(f, file);
end
