function S = comparePreprocRegistration(oldDir, newDir, outDir, opts)
% comparePreprocRegistration  Show that the new preprocessing (asePreProcess.sh) places the ASE
% data on the anatomy the same way the old cluster script (clusterFiles/batch_preprocess) did.
%
%   S = comparePreprocRegistration(oldDir, newDir, outDir)
%     oldDir : old session_ase (Dimon conversion, per-repetition 3dvolreg, old ase_base_al_mat)
%     newDir : new session_ase (dcm2niix conversion, echo-1 3dvolreg, new ase_base_al_mat)
%     outDir : where the PNGs, summary.txt and summary.mat go
%   opts.workDir : scratch folder for AFNI intermediates (default: tempname)
%   opts.series  : [echo1 echo2] series numbers (default [22 24])
%
% Figures (each also summarized in summary.txt):
%   1 conversion   : Dimon and dcm2niix hold the same pixels; only the header origin differs
%   2 displacement : for every brain voxel of the anatomy, the distance (mm) between the raw ASE
%                    sample the old pipeline maps there and the one the new pipeline maps there.
%                    Uses AFNI itself to apply each pipeline's own matrix to index images, so no
%                    coordinate convention is assumed.
%   3 overlay      : each pipeline's ASE base on the anatomy, white-matter / ventricle contours,
%                    checkerboard of old vs new
%   4 ROIs         : ROI membership on the ASE grid (pipeline method: occupancy >= 0.25) via the
%                    old vs the new inverse matrix; weighted Dice per region, against a 1 mm
%                    shift of the segmentation as a reference
%   5 motion       : old per-repetition 3dvolreg parameters vs the new echo-1 run, per acquired
%                    frame, both echoes; and the old vs new motion-corrected images
% Needs AFNI (~/abin or on PATH), Aux_Mat_Functions on the MATLAB path.
    arguments
        oldDir (1,:) char
        newDir (1,:) char
        outDir (1,:) char
        opts.workDir (1,:) char = tempname
        opts.series (1,2) double = [22 24]
        opts.reuse (1,1) logical = false   % skip AFNI steps whose output exists (debugging)
    end
    afni('reuse', opts.reuse);
    if ~isfolder(outDir), mkdir(outDir); end
    if ~isfolder(opts.workDir), mkdir(opts.workDir); end
    W = opts.workDir;
    ser = arrayfun(@(s) sprintf('%03d', s), opts.series, 'UniformOutput', false);
    logLines = {};
    say = @(varargin) fprintf('%s\n', sprintf(varargin{:}));
    function note(varargin)
        msg = sprintf(varargin{:}); logLines{end+1} = msg; fprintf('%s\n', msg);
    end

    C.old = [235 104 52]/255;  C.new = [42 120 214]/255;  C.ref = [150 149 145]/255;
    C.ink = [11 11 11]/255;    C.ink2 = [82 81 78]/255;
    seqMap = interp1([0 0.5 1], [255 255 255; 42 120 214; 16 42 92]/255, linspace(0,1,256));

    O = @(f) fullfile(oldDir, f);  N = @(f) fullfile(newDir, f);  P_ = @(f) fullfile(W, f);

    %% ---------- 1. conversion: same pixels, different header ----------
    note('== 1. Conversion (Dimon vs dcm2niix) ==');
    [rawO, iO] = rd(O(['ase_' ser{1} '.nii']));
    [rawN, iN] = rd(N(['ase_' ser{1} '.nii']));
    AO = niftiAffine(iO.raw);  AN = niftiAffine(iN.raw);
    [Pm, cm] = indexMap(AO, AN, size(rawO,1:3));
    rawOm = remap(rawO, Pm);
    dConv = max(abs(double(rawOm(:)) - double(rawN(:))));
    note('echo 1 (%s): index map old->new = %s, max |old - new| over all %d frames = %g', ser{1}, ...
        mat2str(Pm), size(rawN,4), dConv);
    ctrO = (size(rawO,1:3)' - 1)/2;  ctrN = Pm*ctrO + cm;
    off = AN*[ctrN;1] - AO*[ctrO;1];
    note('header position of the same pixel, new - old (RAS mm): [%.2f %.2f %.2f]', off(1:3));
    S.conv = struct('indexMap', Pm, 'maxAbsDiff', dConv, 'headerOffsetRAS', off(1:3)');
    if isfile(O(['ase_' ser{2} '.nii']))
        r2O = rd(O(['ase_' ser{2} '.nii']));  r2N = rd(N(['ase_' ser{2} '.nii']));
        S.conv.maxAbsDiffEcho2 = max(abs(double(reshape(remap(r2O,Pm),[],1)) - double(r2N(:))));
        note('echo 2 (%s): max |old - new| = %g', ser{2}, S.conv.maxAbsDiffEcho2);
        clear r2O r2N
    end
    % the old base must be raw frame 0 of the old conversion (else the cluster run overwrote it)
    [bO, ibO] = rd(O('ase_base.nii.gz'));  [bN, ibN] = rd(N('ase_base.nii.gz'));
    S.conv.oldBaseIsOldFrame0 = isequal(double(bO), double(rawO(:,:,:,1))) && ...
                                max(abs(niftiAffine(ibO.raw) - AO), [], 'all') < 1e-3;
    S.conv.newBaseIsNewFrame0 = isequal(double(bN), double(rawN(:,:,:,1)));
    note('old ase_base = old raw frame 0 (same header): %d;  new ase_base = new raw frame 0: %d', ...
        S.conv.oldBaseIsOldFrame0, S.conv.newBaseIsNewFrame0);

    fig1(rawO, rawN, rawOm, off, dConv, Pm, orientOf(O(['ase_' ser{1} '.nii'])), orientOf(N(['ase_' ser{1} '.nii'])), C, fullfile(outDir, 'reg1_conversion.png'));
    clear rawO rawOm

    %% ---------- 2. end-to-end displacement on the anatomy ----------
    note('== 2. ASE -> anatomy: old vs new transforms ==');
    anatF = N('anat_fs_brain_al.nii.gz');
    [anat, iA] = rd(anatF);
    anatO = rd(O('anat_fs_brain_al.nii.gz'));
    note('anatomy identical in both folders: %d', isequal(anat, anatO));  clear anatO
    brain = anat > 0;
    matO = O('ase_base_al_mat.aff12.1D');  matN = N('ase_base_al_mat.aff12.1D');
    MO = readAff(matO);  MN = readAff(matN);
    note('matrix difference, linear part max |d| = %.4f; translation d = [%.2f %.2f %.2f] mm', ...
        max(abs(MN(:,1:3) - MO(:,1:3)), [], 'all'), MN(:,4) - MO(:,4));
    for t = {'old', 'new'}
        tag = t{1};
        if strcmp(tag,'old'), base = O('ase_base.nii.gz'); mat = matO; else, base = N('ase_base.nii.gz'); mat = matN; end
        for e = {'i','j','k','1'}
            afni('3dcalc -a %s -expr %s -datum float -prefix %s -overwrite', q(base), e{1}, q(P_(['idx_' tag '_' e{1} '.nii.gz'])));
        end
        afni('3dTcat -prefix %s -overwrite %s %s %s %s', q(P_(['idx_' tag '.nii.gz'])), ...
            q(P_(['idx_' tag '_i.nii.gz'])), q(P_(['idx_' tag '_j.nii.gz'])), q(P_(['idx_' tag '_k.nii.gz'])), q(P_(['idx_' tag '_1.nii.gz'])));
        afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final linear -prefix %s -overwrite', ...
            q(mat), q(P_(['idx_' tag '.nii.gz'])), q(anatF), q(P_(['idx_' tag '_al.nii.gz'])));
        afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final wsinc5 -prefix %s -overwrite', ...
            q(mat), q(base), q(anatF), q(P_(['base_' tag '_al.nii.gz'])));
    end
    xO = single(rd(P_('idx_old_al.nii.gz')));  xN = single(rd(P_('idx_new_al.nii.gz')));
    valid = brain & xO(:,:,:,4) > 0.999 & xN(:,:,:,4) > 0.999;
    vox = sqrt(sum(AN(1:3,1:3).^2, 1));                   % new-grid voxel size per index axis
    d = zeros([nnz(valid) 3], 'single');
    xo = reshape(xO(:,:,:,1:3), [], 3);  xo = xo(valid(:),:);
    xn = reshape(xN(:,:,:,1:3), [], 3);  xn = xn(valid(:),:);
    xo2n = (double(xo)*Pm.') + cm.';                    % old indices in new index convention
    d(:,:) = single((xo2n - double(xn)) .* vox);        % mm along new i, j, k
    dmag = sqrt(sum(d.^2, 2));
    S.disp = struct('median', median(dmag), 'p95', prctile(dmag,95), 'max', max(dmag), ...
        'meanSigned_ijk_mm', mean(d,1), 'voxelSize_mm', vox, 'nVox', numel(dmag), ...
        'brainCoveredFrac', nnz(valid)/nnz(brain));
    note('displacement over %d brain voxels (%.1f%% of the brain inside both ASE fields of view): median %.2f mm, 95th %.2f mm, max %.2f mm', ...
        numel(dmag), 100*S.disp.brainCoveredFrac, S.disp.median, S.disp.p95, S.disp.max);
    note('mean signed shift along ASE axes i, j, k: [%.2f %.2f %.2f] mm; ASE voxel %.2f x %.2f x %.2f mm', mean(d,1), vox);
    dMap = zeros(size(brain), 'single');  dMap(valid) = dmag;
    fig2(dMap, valid, anat, iA, dmag, vox, C, seqMap, fullfile(outDir, 'reg2_displacement.png'));
    clear xO xN xo xn xo2n d

    %% ---------- 3. overlay ----------
    note('== 3. Overlay on the anatomy ==');
    segAF = N('anat_fs_aparc_aseg_al.nii.gz');
    segA = rd(segAF);
    segAO = rd(O('anat_fs_aparc_aseg_al.nii.gz'));
    note('aparc+aseg identical in both folders: %d', isequal(segA, segAO));  clear segAO
    alO = single(rd(P_('base_old_al.nii.gz')));  alN = single(rd(P_('base_new_al.nii.gz')));
    S.overlay.corrOldNew = corr(double(alO(brain)), double(alN(brain)));
    note('correlation of old- vs new-registered ASE base within the brain: %.4f', S.overlay.corrOldNew);
    % the transferred registration outputs must match a re-application of their own matrices
    for t = {'old','new'}
        if strcmp(t{1},'old'), f = O('ase_base_al.nii.gz'); r = alO; else, f = N('ase_base_al.nii.gz'); r = alN; end
        if isfile(f)
            try
                g = single(rd(f));
            catch
                note('%s stored ase_base_al could not be read (truncated transfer?); check skipped', t{1});
                continue
            end
            m = brain & g ~= 0;
            S.overlay.(['corrStored_' t{1}]) = corr(double(g(m)), double(r(m)));
            note('%s stored ase_base_al vs this matrix re-applied: r = %.4f (%.1f%% of brain nonzero in stored file)', ...
                t{1}, S.overlay.(['corrStored_' t{1}]), 100*nnz(m)/nnz(brain));
        end
    end
    wm = ismember(segA, [2 41]);  vent = ismember(segA, [4 5 43 44]);
    fig3(alO, alN, anat, wm, vent, brain, iA, S.overlay.corrOldNew, C, fullfile(outDir, 'reg3_overlay.png'));
    clear alO alN

    %% ---------- 4. ROI membership in ASE space ----------
    note('== 4. ROI membership on the ASE grid ==');
    invO = O('ase_base_al_inv_mat.aff12.1D');  invN = N('ase_base_al_inv_mat.aff12.1D');
    segFiles = {'anat_fs_aparc_aseg_al.nii.gz', 'anat_fs_wmparc_al.nii.gz'};
    R = struct('code',{},'dice',{},'diceRef',{},'nVox',{});
    nVoxN = prod(size(bN,1:3));
    lin = oldToNewLinear(size(bO,1:3), size(bN,1:3), Pm, cm);
    brainAse = false(size(bN,1:3));
    for s = 1:numel(segFiles)
        sf = N(segFiles{s});  [~, nm] = fileparts(erase(segFiles{s}, '.gz'));
        afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final NN -prefix %s -overwrite', q(invO), q(sf), q(sf), q(P_([nm '_old_ase.nii.gz'])));
        afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final NN -prefix %s -overwrite', q(invN), q(sf), q(sf), q(P_([nm '_new_ase.nii.gz'])));
        [sO, isO] = rd(P_([nm '_old_ase.nii.gz']));  [sN, isN] = rd(P_([nm '_new_ase.nii.gz']));
        codes = unique(sN(sN > 0));  LUT = table(double(codes(:)), 'VariableNames', {'code'});
        mArgs = {'maxSampling', [5 5 5], 'minFrac', 0.25};
        MOs = segToAseMembership(sO, LUT, isO, ibO, mArgs{:});
        MNs = segToAseMembership(sN, LUT, isN, ibN, mArgs{:});
        MRs = segToAseMembership(circshift(sN, 1, 1), LUT, isN, ibN, mArgs{:});   % 1 mm shift (1 seg voxel)
        [r_, c_, v_] = find(MOs);  MOs = sparse(r_, lin(c_), v_, height(LUT), nVoxN);
        dice = 2*full(sum(min(MOs,MNs),2)) ./ full(sum(MOs,2) + sum(MNs,2));
        diceRef = 2*full(sum(min(MRs,MNs),2)) ./ full(sum(MRs,2) + sum(MNs,2));
        for k = 1:height(LUT)
            R(end+1) = struct('code', LUT.code(k), 'dice', dice(k), 'diceRef', diceRef(k), 'nVox', full(sum(MNs(k,:)))); %#ok<AGROW>
        end
        brainAse(:) = brainAse(:) | full(any(MNs,1))';
        note('%s: %d regions, weighted Dice old vs new: median %.3f (min %.3f); 1 mm shift reference: median %.3f', ...
            segFiles{s}, height(LUT), median(dice,'omitnan'), min(dice), median(diceRef,'omitnan'));
    end
    S.roi = struct2table(R);
    fig4(S.roi, C, fullfile(outDir, 'reg4_roi_dice.png'));

    %% ---------- 5. motion correction ----------
    note('== 5. Motion correction ==');
    dtes = fileread(firstFile(O(['dcm_hdr_' ser{1} '_dtes.txt']), N(['dcm_hdr_' ser{1} '_dtes.txt'])));
    T = cellfun(@(c) str2double(c{1}), regexp(strsplit(strtrim(dtes), newline), 'e\d_DTE=(\d+)us', 'tokens', 'once'));
    [~, ord] = sortrows([T(:) (1:numel(T))']);          % old: sort -n by dTE, then frame
    idx0 = ord - 1;  r1 = idx0(1:2:end);  r2 = idx0(2:2:end);
    S.motion.r1 = r1;  S.motion.r2 = r2;
    zs = (numel(T)-7):numel(T);                          % last 8 frames = z-shim block (0-based: 90..97)
    note('old split: r1 = %d frames, r2 = %d; z-shim frames (0-based) in r1: %s, in r2: %s', numel(r1), numel(r2), ...
        mat2str(intersect(r1, zs-1)'), mat2str(intersect(r2, zs-1)'));
    if isfile(O(['ase_' ser{1} '_r1.nii.gz']))
        x = rd(O(['ase_' ser{1} '_r1.nii.gz']));  [rawO] = rd(O(['ase_' ser{1} '.nii']));
        S.motion.splitVerified = isequal(x, rawO(:,:,:,r1+1));
        note('old ase_%s_r1 = raw frames r1 (sorted split reproduced): %d', ser{1}, S.motion.splitVerified);
        clear x rawO
    end
    % current pipeline (motionMode = echo1), regenerated here so echo 2 reflects today's code
    afni('3dvolreg -Fourier -twopass -twodup -base %s -dfile %s -1Dmatrix_save %s -prefix %s -overwrite %s', ...
        q(N('ase_base.nii.gz')), q(P_('e1_mc.1D')), q(P_('e1_mat.aff12.1D')), q(P_('e1_proc.nii.gz')), q(N(['ase_' ser{1} '.nii'])));
    afni('3dAllineate -1Dmatrix_apply %s -source %s -master %s -final wsinc5 -prefix %s -overwrite', ...
        q(P_('e1_mat.aff12.1D')), q(N(['ase_' ser{2} '.nii'])), q(N('ase_base.nii.gz')), q(P_('e2_proc.nii.gz')));
    mcN = load(P_('e1_mc.1D'));
    if isfile(N('ase_e1_mc.1D'))
        mcS = load(N('ase_e1_mc.1D'));
        note('re-run echo-1 motion parameters vs stored ase_e1_mc.1D: max |d| = %.4f', max(abs(mcN(:,2:7) - mcS(:,2:7)), [], 'all'));
    end
    mc = struct();
    for e = 1:2
        for rr = 1:2
            f = O(sprintf('ase_%s_mc_r%d.1D', ser{e}, rr));
            if isfile(f), mc.(sprintf('e%d_r%d', e, rr)) = load(f); end
        end
    end
    pn = {'roll (°)','pitch (°)','yaw (°)','dS (mm)','dL (mm)','dP (mm)'};
    for e = 1:2
        for rr = 1:2
            k = sprintf('e%d_r%d', e, rr);
            if ~isfield(mc,k), continue; end
            fr = S.motion.(sprintf('r%d', rr));
            dd = mc.(k)(:,2:7) - mcN(fr+1, 2:7);
            S.motion.(['rmsDiff_' k]) = sqrt(mean(dd.^2, 1));
            note('old %s vs new echo-1 params at the same frames, RMS diff [%s]', k, ...
                strjoin(compose('%s %.3f', string(pn'), S.motion.(['rmsDiff_' k])'), ', '));
        end
    end
    % image agreement: old vs new motion-corrected frames, inside the ASE brain mask
    nrm = NaN(numel(T), 2);
    for e = 1:2
        if e == 1, newProc = single(rd(P_('e1_proc.nii.gz'))); else, newProc = single(rd(P_('e2_proc.nii.gz'))); end
        for rr = 1:2
            f = O(sprintf('ase_%s_reg_r%d.nii.gz', ser{e}, rr));
            if ~isfile(f), continue; end
            try
                oldReg = remap(single(rd(f)), Pm);
            catch
                note('could not read %s (transfer incomplete?); skipped', f);  continue
            end
            fr = S.motion.(sprintf('r%d', rr));
            for n = 1:numel(fr)
                a = oldReg(:,:,:,n);  b = newProc(:,:,:,fr(n)+1);
                nrm(fr(n)+1, e) = norm(a(brainAse) - b(brainAse)) / norm(b(brainAse));
            end
        end
        clear newProc oldReg
    end
    S.motion.imageNRMS = nrm;
    note('old vs new motion-corrected images, normalized RMS difference in brain: echo 1 median %.4f (max %.4f), echo 2 median %.4f (max %.4f)', ...
        median(nrm(:,1),'omitnan'), max(nrm(:,1)), median(nrm(:,2),'omitnan'), max(nrm(:,2)));
    fig5(mc, mcN, S.motion, nrm, pn, zs-1, C, fullfile(outDir, 'reg5_motion.png'));

    %% ---------- write ----------
    S.log = logLines;
    S.takeaways = takeaways(S);
    fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
    fprintf(fid, '%s\n', S.takeaways{:});
    fprintf(fid, '\n=============== DETAILS (run log) ===============\n');
    fprintf(fid, '%s\n', logLines{:});  fclose(fid);
    save(fullfile(outDir, 'summary.mat'), 'S');
    say('Wrote figures and summary to %s', outDir);
end

function L = takeaways(S)
% Main conclusions for slides, filled in from the computed numbers
    L = {'=============== MAIN TAKEAWAYS ===============', ''};
    L{end+1} = ['OVERALL: The new preprocessing places the ASE data on the anatomy within ' ...
        sprintf('%.1f mm of the old (max over the brain); where they differ (echo-2 motion), the new pipeline removes an error in the old.', S.disp.max)];
    L{end+1} = '';
    c = S.conv;
    e2 = '';  if isfield(c, 'maxAbsDiffEcho2'), e2 = sprintf(' (echo 2: %g)', c.maxAbsDiffEcho2); end
    L{end+1} = sprintf(['Fig 1 (reg1_conversion): Dimon (old) and dcm2niix (new) write identical pixels: max |difference| = %g ' ...
        'over all frames%s. Only the header position differs, by [%.1f %.1f %.1f] mm (R, A, S).'], c.maxAbsDiff, e2, c.headerOffsetRAS);
    d = S.disp;
    L{end+1} = sprintf(['Fig 2 (reg2_displacement): Each pipeline''s registration cancels its own header offset. For every brain voxel, ' ...
        'old and new point to raw ASE positions %.2f mm apart (median; 95th %.2f, max %.2f mm) against a %.2f x %.2f x %.1f mm ASE voxel.'], ...
        d.median, d.p95, d.max, d.voxelSize_mm);
    L{end+1} = sprintf(['Fig 3 (reg3_overlay): The old- and new-registered ASE images are indistinguishable on the anatomy ' ...
        '(r = %.4f within the brain); white-matter and ventricle contours fall in the same place.'], S.overlay.corrOldNew);
    if isfield(S, 'roi')
        T = S.roi;  ok = T.nVox > 0;
        L{end+1} = sprintf(['Fig 4 (reg4_roi_dice): ROI masks differ between pipelines less than they would if the segmentation were ' ...
            'nudged by 1 mm, in %.0f%% of regions (median Dice old vs new %.2f; 1 mm shift %.2f). Small regions score lower because they are ' ...
            'mostly boundary. The remaining difference is likely the 1 mm nearest-neighbour rounding in the ROI method, not registration ' ...
            '(not yet tested separately).'], 100*mean(T.dice(ok) > T.diceRef(ok)), median(T.dice(ok)), median(T.diceRef(ok)));
    end
    if isfield(S, 'motion')
        m = S.motion;  s5 = 'Fig 5 (reg5_motion): Echo-1 motion parameters are the same in both pipelines';
        if isfield(m, 'rmsDiff_e1_r2'), s5 = [s5 sprintf(' (RMS difference <= %.3f deg / mm)', max([m.rmsDiff_e1_r1 m.rmsDiff_e1_r2]))]; end
        s5 = [s5 '.'];
        if isfield(m, 'rmsDiff_e2_r1')
            r2 = max([m.rmsDiff_e2_r1; m.rmsDiff_e2_r2], [], 1);
            s5 = [s5 sprintf([' The old pipeline registered echo 2 on its own against the echo-1 base, which put it off the echo-1 motion ' ...
                'by up to %.2f deg pitch and %.2f mm dS (RMS); the new pipeline applies the echo-1 transforms to echo 2, since both echoes ' ...
                'come from the same excitation.'], r2(2), r2(4))];
        end
        if isfield(m, 'imageNRMS') && any(isfinite(m.imageNRMS(:)))
            s5 = [s5 sprintf(' Motion-corrected images, old vs new: median RMS difference %.2f%% (echo 1), %.2f%% (echo 2) of the signal.', ...
                100*median(m.imageNRMS(:,1), 'omitnan'), 100*median(m.imageNRMS(:,2), 'omitnan'))];
        end
        L{end+1} = s5;
    end
    L{end+1} = '';
    L{end+1} = ['CAVEAT (both pipelines, see docs/motion-correction-plan.md): the echo-1 motion estimates follow the dTE order ' ...
        '(dTE explains 83% of dP) and the z-shim block shows up to 2.5 deg of apparent pitch, i.e. 3dvolreg partly reads ASE contrast as motion. ' ...
        'This does not separate old from new.'];
    L{end+1} = ['NOTE: ROI membership uses the fixed segToAseMembership (2026-11-20); before the fix, ROI masks were shifted ' ...
        'about one slice and matched AFNI resampling in only 38% of labelled voxels.'];
end

%% ===================== helpers =====================
function [v, info] = rd(f)
    v = niftiread(f);  info = niftiinfo(f);
end

function s = q(p), s = ['"' p '"']; end

function f = firstFile(varargin)
    for k = 1:numel(varargin), if isfile(varargin{k}), f = varargin{k}; return; end, end
    error('None of the files exist: %s', strjoin(varargin, ', '));
end

function afni(fmt, varargin)
    persistent reuse
    if strcmp(fmt, 'reuse'), reuse = varargin{1}; return; end
    cmd = sprintf(fmt, varargin{:});
    tok = regexp(cmd, '-prefix "([^"]+)"', 'tokens', 'once');
    if ~isempty(reuse) && reuse && ~isempty(tok) && isfile(tok{1}), return; end
    [st, out] = system(['export PATH="$HOME/abin:$PATH"; export AFNI_NIFTI_TYPE_WARN=NO; ' cmd]);
    if st ~= 0, error('AFNI command failed:\n%s\n%s', cmd, out); end
end

function M = readAff(f)
    t = fileread(f);  t = regexprep(t, '#[^\n]*', '');
    v = sscanf(t, '%f');  M = reshape(v(1:12), 4, 3)';
end

function [P, c] = indexMap(AO, AN, dimsO)
% signed permutation P and offset c with idx_new = P*idx_old + c (0-based) for the same pixel,
% from the direction parts of the two affines (the origins differ, so c comes from the flips)
    RO = AO(1:3,1:3) ./ sqrt(sum(AO(1:3,1:3).^2,1));
    RN = AN(1:3,1:3) ./ sqrt(sum(AN(1:3,1:3).^2,1));
    P = round(RN \ RO);
    assert(all(abs(P(:)) <= 1) && all(sum(abs(P),1) == 1) && all(sum(abs(P),2) == 1), ...
        'Old and new grids are not related by axis flips/permutations: %s', mat2str(RN\RO, 3));
    c = zeros(3,1);
    for d = 1:3
        s = find(P(d,:));
        if P(d,s) < 0, c(d) = dimsO(s) - 1; end
    end
end

function B = remap(A, P)
% reorder an old-grid array (3D or 4D) into the new grid's index order
    src = zeros(1,3);  for d = 1:3, src(d) = find(P(d,:)); end
    B = permute(A, [src 4:max(4,ndims(A))]);
    for d = 1:3, if P(d,src(d)) < 0, B = flip(B, d); end, end
end

function lin = oldToNewLinear(dO, dN, P, c)
    [I,J,K] = ndgrid(0:dO(1)-1, 0:dO(2)-1, 0:dO(3)-1);
    x = P*[I(:) J(:) K(:)]' + c;
    lin = sub2ind(dN, x(1,:)+1, x(2,:)+1, x(3,:)+1)';
end

function [V, vs] = toRAS(V, A)
% reorder a volume so dims 1..3 run along +R, +A, +S (for display)
    R = A(1:3,1:3);  perm = zeros(1,3);  sg = zeros(1,3);
    for w = 1:3, [~, perm(w)] = max(abs(R(w,:))); sg(w) = sign(R(w,perm(w))); end
    V = permute(V, [perm 4:max(4,ndims(V))]);
    for w = 1:3, if sg(w) < 0, V = flip(V, w); end, end
    vs = sqrt(sum(R(:,perm).^2, 1));
end

function im = sl(V, view, k)
% 2-D display slice of an RAS-ordered volume: axial (A up), coronal / sagittal (S up)
    switch view
        case 'ax',  im = flipud(V(:,:,k)');
        case 'cor', im = flipud(squeeze(V(:,k,:))');
        case 'sag', im = flipud(squeeze(V(k,:,:))');
    end
end

function k = ctrSlices(mask)
    [x,y,z] = ind2sub(size(mask), find(mask));
    k = round([median(z) median(y) median(x)]);          % ax, cor, sag
end

function styleAx(ax, C)
    set(ax, 'XColor', C.ink2, 'YColor', C.ink2, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.75);
end

function sv(f, file)
    exportgraphics(f, file, 'Resolution', 200, 'BackgroundColor', 'white');  close(f);
end

%% ===================== figures =====================
function fig1(rawO, rawN, rawOm, off, dConv, P, orO, orN, C, file)
    k = round(size(rawN,3)/2);
    g = @(V) rot90(double(V(:,:,k,1)));
    imgs = {g(rawO), g(rawN), g(rawOm), abs(g(rawOm) - g(rawN))};
    ttl = {sprintf('Old: Dimon (%s), as stored', orO), sprintf('New: dcm2niix (%s), as stored', orN), ...
           'Old, reordered to the new grid', '|reordered old − new|'};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1700 520]);
    tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    cl = [0 prctile(imgs{2}(:), 99.5)];
    for p = 1:4
        ax = nexttile(tl);
        if p < 4, imagesc(ax, imgs{p}, cl); colormap(ax, gray); else, imagesc(ax, imgs{p}, [0 1]); colormap(ax, gray); end
        axis(ax, 'image', 'off');  title(ax, ttl{p}, 'FontWeight', 'normal', 'Color', C.ink, 'FontSize', 12, 'Interpreter', 'none');
    end
    title(tl, sprintf(['Conversion: identical pixels (max |difference| over all frames = %g); ' ...
        'only the header position differs, by [%.1f, %.1f, %.1f] mm (R, A, S)'], dConv, off(1:3)), ...
        'FontSize', 13, 'Color', C.ink, 'Interpreter', 'none');
    subtitle(tl, sprintf('Echo 1, slice %d, frame 0. Index map old→new: %s', k, mat2str(P)), 'Color', C.ink2, 'FontSize', 11);
    sv(f, file);
end

function s = orientOf(f)
    [~, s] = system(sprintf('export PATH="$HOME/abin:$PATH"; 3dinfo -orient "%s"', f));  s = strtrim(s);
end

function fig2(dMap, valid, anat, iA, dmag, vox, C, seqMap, file)
    A = niftiAffine(iA.raw);
    D = toRAS(dMap, A);  B = toRAS(single(anat), A);  Vd = toRAS(valid, A);
    k = ctrSlices(Vd);  views = {'ax','cor','sag'};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1800 560]);
    tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    cmax = 0.5;                                          % mm; histogram shows the voxel scale
    bl = prctile(B(B>0), 99);
    for p = 1:3
        ax = nexttile(tl);
        a = sl(B, views{p}, k(p)) / bl;  d = sl(D, views{p}, k(p));  m = sl(Vd, views{p}, k(p));
        rgbA = repmat(min(a,1), 1, 1, 3);                      % anatomy in grey outside the mask
        idx = max(1, min(256, round(d / cmax * 255) + 1));
        rgbD = reshape(seqMap(idx,:), [size(d) 3]);
        rgb = rgbA;  mm = repmat(m, 1, 1, 3);  rgb(mm) = rgbD(mm);
        image(ax, rgb);  axis(ax, 'image', 'off');
        vt = {'Axial','Coronal','Sagittal'};  title(ax, vt{p}, 'FontWeight', 'normal', 'Color', C.ink);
    end
    colormap(ax, seqMap);  clim(ax, [0 cmax]);
    cb = colorbar(ax, 'eastoutside');  cb.Label.String = 'old − new distance (mm)';  cb.Color = C.ink2;
    ax = nexttile(tl);  hold(ax, 'on');
    edges = 0:0.05:max(5, ceil(max(dmag)*10)/10);
    histogram(ax, dmag, edges, 'FaceColor', C.new, 'EdgeColor', 'w', 'FaceAlpha', 1, 'Normalization', 'probability');
    yl = ylim(ax);
    for v = [vox(1) vox(3)]
        plot(ax, [v v], yl, '--', 'Color', C.ink2, 'LineWidth', 1);
    end
    text(ax, vox(1), yl(2)*0.95, sprintf(' in-plane voxel %.2f mm', vox(1)), 'Color', C.ink2, 'FontSize', 10);
    text(ax, vox(3), yl(2)*0.85, sprintf(' slice %.1f mm', vox(3)), 'Color', C.ink2, 'FontSize', 10);
    xlabel(ax, 'distance (mm)');  ylabel(ax, 'fraction of brain voxels');  styleAx(ax, C);  xlim(ax, [0 edges(end)]);
    title(ax, sprintf('median %.2f mm · 95th %.2f · max %.2f', median(dmag), prctile(dmag,95), max(dmag)), ...
        'FontWeight', 'normal', 'Color', C.ink);
    title(tl, 'Where each pipeline puts the same ASE sample on the anatomy', 'FontSize', 13, 'Color', C.ink);
    subtitle(tl, 'Distance between the old and new pipelines'' mapping of every brain voxel back to the raw ASE data (each pipeline''s own header + matrix)', ...
        'Color', C.ink2, 'FontSize', 11);
    sv(f, file);
end

function fig3(alO, alN, anat, wm, vent, brain, iA, r, C, file)
    A = niftiAffine(iA.raw);
    O = toRAS(alO, A);  N = toRAS(alN, A);  W = toRAS(wm, A);  V = toRAS(vent, A);  Bm = toRAS(brain, A);
    k = ctrSlices(Bm);  views = {'ax','cor','sag'};  vt = {'Axial','Coronal','Sagittal'};
    sO = prctile(alO(brain), 99);  sN = prctile(alN(brain), 99);
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1500 1450]);
    tl = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    cols = {'Old pipeline','New pipeline','Checkerboard of old and new'};
    for p = 1:3
        o = sl(O, views{p}, k(p)) / sO;  n = sl(N, views{p}, k(p)) / sN;
        w = sl(W, views{p}, k(p));  v = sl(V, views{p}, k(p));
        [rr, cc] = find(sl(Bm, views{p}, k(p)));  pad = 8;      % crop to the brain
        rr = max(1, min(rr)-pad):min(size(o,1), max(rr)+pad);  cc = max(1, min(cc)-pad):min(size(o,2), max(cc)+pad);
        o = o(rr,cc);  n = n(rr,cc);  w = w(rr,cc);  v = v(rr,cc);
        [cy, cx] = ndgrid(1:size(o,1), 1:size(o,2));  cb = mod(floor(cy/16) + floor(cx/16), 2) == 0;
        ims = {o, n, o.*cb + n.*~cb};
        for c = 1:3
            ax = nexttile(tl);  imagesc(ax, min(ims{c},1), [0 1]);  colormap(ax, gray);  axis(ax, 'image', 'off');  hold(ax, 'on');
            contour(ax, double(w), [0.5 0.5], 'EdgeColor', C.new, 'LineWidth', 0.9);
            contour(ax, double(v), [0.5 0.5], 'EdgeColor', C.old, 'LineWidth', 0.9);
            if p == 1, title(ax, cols{c}, 'FontWeight', 'normal', 'Color', C.ink, 'FontSize', 13, 'Interpreter', 'none'); end
            if c == 1, text(ax, -8, size(o,1)/2, vt{p}, 'Rotation', 90, 'HorizontalAlignment', 'center', 'Color', C.ink2, 'FontSize', 12); end
        end
    end
    title(tl, sprintf('ASE base image on the anatomy, old vs new registration (r = %.4f within the brain)', r), 'FontSize', 14, 'Color', C.ink);
    subtitle(tl, 'Contours from FreeSurfer: blue = cerebral white matter, orange = lateral ventricles', 'Color', C.ink2, 'FontSize', 12);
    sv(f, file);
end

function fig4(T, C, file)
    cls = repmat("Subcortical / other", height(T), 1);
    cls(T.code >= 1000 & T.code < 3000) = "Cortex";
    cls(T.code >= 3000 & T.code < 5003 | ismember(T.code, [2 41])) = "White matter";
    groups = ["Cortex","White matter","Subcortical / other"];
    col = {[42 120 214]/255, [235 104 52]/255, [98 80 214]/255};
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1500 620]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax = nexttile(tl);  hold(ax, 'on');
    ok = T.nVox > 0;
    for g = 1:3
        m = ok & cls == groups(g);
        scatter(ax, T.nVox(m), T.dice(m), 36, col{g}, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'DisplayName', groups(g));
    end
    set(ax, 'XScale', 'log');  ylim(ax, [min([0.8; T.dice(ok); T.diceRef(ok)]) - 0.02 1.005]);
    xlabel(ax, 'region size on the ASE grid (effective voxels)');  ylabel(ax, 'weighted Dice, old vs new');
    legend(ax, 'Location', 'southeast', 'Box', 'off');  styleAx(ax, C);
    title(ax, sprintf('Old vs new inverse transform: median Dice %.3f', median(T.dice(ok))), 'FontWeight', 'normal', 'Color', C.ink);
    ax = nexttile(tl);  hold(ax, 'on');
    scatter(ax, T.diceRef(ok), T.dice(ok), 30, C.new, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
    lim = [min([T.dice(ok); T.diceRef(ok)]) - 0.02 1.005];
    plot(ax, lim, lim, '--', 'Color', C.ink2, 'LineWidth', 1);
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);
    xlabel(ax, 'Dice, new vs new with segmentation shifted 1 mm');  ylabel(ax, 'Dice, old vs new');  styleAx(ax, C);
    above = mean(T.dice(ok) > T.diceRef(ok));
    title(ax, sprintf('Old vs new agree better than a 1 mm shift in %.0f%% of regions', 100*above), 'FontWeight', 'normal', 'Color', C.ink);
    title(tl, 'ROI membership on the ASE grid (pipeline method, occupancy ≥ 0.25): old vs new transform', 'FontSize', 13, 'Color', C.ink);
    sv(f, file);
end

function fig5(mc, mcN, M, nrm, pn, zs, C, file)
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1800 1150]);
    tl = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    fr = 0:size(mcN,1)-1;
    for p = 1:6
        ax = nexttile(tl);  hold(ax, 'on');
        patch(ax, [zs(1)-0.5 zs(end)+0.5 zs(end)+0.5 zs(1)-0.5], [-1e3 -1e3 1e3 1e3], [0.95 0.95 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, fr, mcN(:,p+1), '-', 'Color', C.new, 'LineWidth', 2, 'DisplayName', 'New: echo 1, one run (applied to echo 2)');
        mk = {'o','s'};
        for rr = 1:2
            k = sprintf('e1_r%d', rr);
            if isfield(mc,k), plot(ax, M.(sprintf('r%d',rr)), mc.(k)(:,p+1), mk{rr}, 'Color', C.old, 'MarkerSize', 5, 'LineWidth', 1, ...
                    'DisplayName', sprintf('Old: echo 1, repetition split %d', rr)); end
            k = sprintf('e2_r%d', rr);
            if isfield(mc,k), plot(ax, M.(sprintf('r%d',rr)), mc.(k)(:,p+1), mk{rr}, 'Color', C.ref, 'MarkerSize', 4, 'LineWidth', 0.75, ...
                    'DisplayName', sprintf('Old: echo 2 (own 3dvolreg), split %d', rr)); end
        end
        v = [mcN(:,p+1); cell2mat(cellfun(@(c) mc.(c)(:,p+1), fieldnames(mc), 'UniformOutput', false))];
        pad = 0.1*max(range(v), 0.05);  ylim(ax, [min(v)-pad max(v)+pad]);  xlim(ax, [-1 fr(end)+1]);
        ylabel(ax, pn{p});  styleAx(ax, C);
        if p > 3, xlabel(ax, 'acquired frame'); end
        if p == 1, legend(ax, 'Location', 'northoutside', 'Box', 'off', 'NumColumns', 1, 'FontSize', 10); end
    end
    ax = nexttile(tl, [1 3]);  hold(ax, 'on');
    patch(ax, [zs(1)-0.5 zs(end)+0.5 zs(end)+0.5 zs(1)-0.5], [0 0 1e3 1e3], [0.95 0.95 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, fr, 100*nrm(:,1), 'o-', 'Color', C.old, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'Echo 1');
    plot(ax, fr, 100*nrm(:,2), 's-', 'Color', C.ref, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'Echo 2');
    yl = [0 max(1, 1.1*max(100*nrm(:)))];  ylim(ax, yl);  xlim(ax, [-1 fr(end)+1]);
    xlabel(ax, 'acquired frame');  ylabel(ax, 'old − new, % of signal (RMS)');  styleAx(ax, C);
    legend(ax, 'Location', 'northwest', 'Box', 'off');
    title(ax, 'Motion-corrected images, old vs new, inside the brain (shaded: z-shim block)', 'FontWeight', 'normal', 'Color', C.ink);
    title(tl, 'Motion correction: old per-repetition 3dvolreg vs new single echo-1 run (same base: raw frame 0)', 'FontSize', 13, 'Color', C.ink);
    sv(f, file);
end
