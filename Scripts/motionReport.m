function R = motionReport(sessDir, opts)
% motionReport  Per-subject motion-correction decision figure (see docs/motion-correction-plan.md).
%
%   R = motionReport(sessDir)            sessDir = <subjID>/session_ase, after asePreProcess.sh
%   R = motionReport(sessDir, 'anchorMaxDTE', 8, 'threshLow', 0.3, 'threshHigh', 1, ...)
%
% Motion is estimated only on "anchor" frames, the sweep frames with dTE <= anchorMaxDTE (ms),
% which have nearly spin-echo contrast. Transforms for all other frames are interpolated
% linearly in acquisition order; frames after the last anchor (including the whole z-shim block)
% hold the last anchor's transform. Echo 2 would reuse echo 1's transforms.
%
% Decision metric: the largest displacement of any brain voxel, over all frames, in voxel units
% (each displacement component divided by the voxel size along that axis, so 1 = one voxel in
% any direction). With linear interpolation the maximum always falls on an anchor.
%   metric <  threshLow               -> no correction
%   threshLow <= metric < threshHigh  -> correct with the interpolated transforms
%   metric >= threshHigh              -> highAction: 'correct' or 'exclude' (frames above
%                                        threshHigh are marked)
%
% The figure also shows, in grey, the current pipeline's per-frame 3dvolreg of every frame, whose
% short/long-tau contrast changes are partly read as motion.
%
% Brain mask: sessDir/<maskName> if present, else 3dAutomask of ase_base.nii.gz (in workDir).
% Inputs read from sessDir: ase_hdr.json, ase_base.nii.gz, the raw echo-1 NIfTI (ase_<series>.nii),
% and the saved echo-1 matrices ase_e1_mc_mat.aff12.1D if present (otherwise 3dvolreg is rerun on
% echo 1 in opts.workDir, ~2 min; per-frame results are independent, so the anchor transforms are
% the same as registering the anchors alone).
% Needs AFNI (~/abin or on PATH) and Aux_Mat_Functions on the MATLAB path.
    arguments
        sessDir (1,:) char
        opts.anchorMaxDTE (1,1) double {mustBeNonnegative} = 8
        opts.threshLow (1,1) double {mustBePositive} = 0.3     % voxels; placeholder
        opts.threshHigh (1,1) double {mustBePositive} = 1.0    % voxels; placeholder
        opts.highAction (1,:) char {mustBeMember(opts.highAction, {'correct','exclude'})} = 'correct'
        opts.outFile (1,:) char = ''
        opts.workDir (1,:) char = ''
        opts.volregOpts (1,:) char = '-Fourier -twopass -twodup'   % asePreProcess.sh volregOpts
        opts.matName (1,:) char = 'ase_e1_mc_mat.aff12.1D'          % asePreProcess.sh motionMatName
        opts.maskName (1,:) char = 'ase_brain_mask.nii.gz'          % brain mask on the ASE grid, if any
    end
    assert(opts.threshHigh > opts.threshLow, 'threshHigh must be larger than threshLow.');
    F = @(f) fullfile(sessDir, f);
    h = jsondecode(fileread(F('ase_hdr.json')));
    subj = '';  if isfield(h, 'subject'), subj = h.subject; end
    if isempty(opts.outFile), opts.outFile = F('motion_report.png'); end
    if isempty(opts.workDir), opts.workDir = fullfile(tempdir, ['motionReport_' subj]); end
    W = opts.workDir;  if ~isfolder(W), mkdir(W); end

    % ---------- frames ----------
    [dTE, ~, ~, info] = readASEHeader(F('ase_hdr.json'));
    nF = numel(dTE);  zs = info.zShimInds(:);  sw = setdiff((1:nF)', zs);
    an = sw(dTE(sw) <= opts.anchorMaxDTE);
    assert(~isempty(an), 'No sweep frame has dTE <= %g ms.', opts.anchorMaxDTE);
    e = h.echos;  if iscell(e), e = e{1}; else, e = e(1); end
    raw = F(sprintf('ase_%s.nii', e.series));
    TR = NaN;
    js = F(sprintf('dcm_hdr_%s.json', e.series));
    if isfile(js), d = jsondecode(fileread(js)); if isfield(d, 'RepetitionTime'), TR = d.RepetitionTime; end, end

    % ---------- echo-1 transforms (every frame; anchors are a subset) ----------
    matFile = F(opts.matName);
    if ~isfile(matFile)
        matFile = fullfile(W, 'e1_mat.aff12.1D');
        if ~isfile(matFile)
            fprintf('[%s] running 3dvolreg on echo 1 (%s) -> %s\n', subj, raw, W);
            afni('3dvolreg %s -base %s -dfile %s -1Dmatrix_save %s -prefix NULL %s', opts.volregOpts, ...
                q(F('ase_base.nii.gz')), q(fullfile(W, 'e1_par.1D')), q(matFile), q(raw));
        end
    end
    M = readMats(matFile);
    assert(size(M,1) == nF, '%s has %d matrices, header lists %d frames.', matFile, size(M,1), nF);
    Mi = interpMats(M, an, nF);

    % ---------- brain mask and coordinates ----------
    % the pipeline's brain mask on the ASE grid if it exists, else 3dAutomask of the base
    mk = F(opts.maskName);
    if ~isfile(mk)
        mk = fullfile(W, 'automask.nii.gz');
        if ~isfile(mk), afni('3dAutomask -prefix %s -overwrite %s', q(mk), q(F('ase_base.nii.gz'))); end
    end
    mask = niftiread(mk) > 0;  base = double(niftiread(F('ase_base.nii.gz')));
    ii = niftiinfo(F('ase_base.nii.gz'));  A = niftiAffine(ii.raw);
    [i1, i2, i3] = ind2sub(size(mask), find(mask));
    ras = A * [[i1 i2 i3]-1 ones(numel(i1),1)]';                % 0-based indices -> RAS mm
    xyz = (diag([-1 -1 1]) * ras(1:3,:))';                      % AFNI DICOM (RAI) mm
    toVox = A(1:3,1:3) \ diag([-1 -1 1]);                       % DICOM displacement -> voxel units
    vox = sqrt(sum(A(1:3,1:3).^2));

    dFull = maxDisp(M, xyz, toVox);
    dInt  = maxDisp(Mi, xyz, toVox);
    metric = max(dInt);  [~, worst] = max(dInt);
    dFullMM = maxDisp(M, xyz, eye(3));  dIntMM = maxDisp(Mi, xyz, eye(3));   % same, in mm
    metricMM = max(dIntMM);
    if metric < opts.threshLow,       decision = 'none';
    elseif metric < opts.threshHigh,  decision = 'correct';
    else,                             decision = opts.highAction;
    end
    excl = false(nF,1);
    if strcmp(decision, 'exclude'), excl = dInt >= opts.threshHigh; end
    [~, dMapV] = maxDisp(Mi(worst,:), xyz, toVox);
    dMap = nan(size(mask));  dMap(mask) = dMapV;
    % 3dvolreg parameters for display (pipeline's -dfile, else ours, else from the matrices)
    parFile = '';
    for c = {F('ase_e1_mc.1D'), fullfile(W, 'e1_par.1D')}, if isfile(c{1}), parFile = c{1}; break, end, end
    if ~isempty(parFile)
        pFull = load(parFile);  pFull = pFull(:, 2:7);   % -dfile: n roll pitch yaw dS dL dP rms rms
    else
        pFull = matToParams(M);
    end
    pInt = interpMats(pFull, an, nF);   % same linear interpolation as the matrices

    R = struct('subject', subj, 'decision', decision, 'metric_vox', metric, 'metric_mm', metricMM, 'worstFrame', worst-1, ...
        'anchors', an-1, 'anchorMaxDTE', opts.anchorMaxDTE, 'maxGap', max(diff([an; sw(end)+1])), ...
        'threshLow', opts.threshLow, 'threshHigh', opts.threshHigh, 'highAction', opts.highAction, ...
        'excludedFrames', find(excl)-1, 'dispInterp_vox', dInt, 'dispFull_vox', dFull, 'dispInterp_mm', dIntMM, 'dispFull_mm', dFullMM, ...
        'mask', mk, ...
        'fullSweepMax_vox', max(dFull(sw)), 'fullZshimMax_vox', max(dFull(zs)), 'voxel_mm', vox, 'TR', TR);
    fprintf('[%s] anchors dTE<=%g ms: %d frames, max gap %d; max displacement %.3f vox (frame %d) -> %s\n', ...
        subj, opts.anchorMaxDTE, numel(an), R.maxGap, metric, worst-1, decision);

    drawFigure(R, dTE, an, sw, zs, pFull, pInt, dMap, base, mask, A, opts.outFile);
end

%% ===================== figure =====================
function drawFigure(R, dTE, an, sw, zs, pFull, pInt, dMap, base, mask, A, file)
    C.int = [42 120 214]/255;  C.full = [150 149 145]/255;  C.ink = [11 11 11]/255;  C.ink2 = [82 81 78]/255;
    C.shade = [0.95 0.95 0.94];  C.excl = [214 69 65]/255;
    nF = numel(dTE);  fr = (0:nF-1)';
    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1700 1300]);
    % three independent rows (normalised [left bottom width height]), so the timeline's right-axis
    % labels do not narrow the rows below
    rowPos = {[0.01 0.645 0.98 0.315], [0.01 0.335 0.98 0.30], [0.01 0.01 0.98 0.315]};

    % --- displacement timeline (voxels; the decision metric) ---
    % own nested layout, so the right-axis labels don't narrow the rows below
    t0 = tiledlayout(f, 1, 1, 'Padding', 'tight', 'OuterPosition', rowPos{1});
    ax = nexttile(t0);  yyaxis(ax, 'left');  hold(ax, 'on');
    yl = [0 1.08*max([R.dispFull_vox; R.threshHigh*1.15])];
    zband(ax, zs, yl, C);
    for t = [R.threshLow R.threshHigh]
        plot(ax, [-1 nF], [t t], '--', 'Color', C.ink2, 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    plot(ax, fr, R.dispFull_vox, 'o', 'Color', C.full, 'MarkerSize', 4, 'LineWidth', 1, ...
        'DisplayName', 'every frame registered (current pipeline; contrast partly read as motion)');
    plot(ax, fr, R.dispInterp_vox, '-', 'Color', C.int, 'LineWidth', 2, ...
        'DisplayName', sprintf('anchors (dTE \\leq %g ms) + interpolation', R.anchorMaxDTE));
    plot(ax, fr(an), R.dispInterp_vox(an), 'o', 'MarkerFaceColor', C.int, 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'anchor frames');
    if ~isempty(R.excludedFrames)
        x = R.excludedFrames;
        plot(ax, x, R.dispInterp_vox(x+1), 'x', 'Color', C.excl, 'MarkerSize', 9, 'LineWidth', 2, 'DisplayName', 'excluded');
    end
    xlim(ax, [-1 nF]);  ylim(ax, yl);  ylabel(ax, 'max displacement (voxels)');  styleAx(ax, C);
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 10);
    title(ax, 'Largest displacement of any brain voxel, per frame, relative to frame 0 (shaded: z-shim block)', ...
        'FontWeight', 'normal', 'Color', C.ink);
    xt = 0:10:nF-1;  xticks(ax, xt);
    if isnan(R.TR)
        xlabel(ax, 'acquired frame');
    else
        xticklabels(ax, arrayfun(@(x) sprintf('%d\\newline\\color[rgb]{0.5,0.5,0.5}%.0f s', x, x*R.TR), xt, 'UniformOutput', false));
        set(ax, 'TickLabelInterpreter', 'tex');
        xlabel(ax, sprintf('acquired frame  /  \\color[rgb]{0.5,0.5,0.5}acquisition time (s, TR %.1f s)', R.TR));
    end
    % threshold labels outside the plot, on the right axis
    yyaxis(ax, 'right');  ylim(ax, yl);
    set(ax, 'YTick', [R.threshLow R.threshHigh], 'YTickLabel', ...
        {sprintf('no correction below %.2g', R.threshLow), sprintf('%s above %.2g', R.highAction, R.threshHigh)});
    ax.YAxis(1).Color = C.ink2;  ax.YAxis(2).Color = C.ink2;
    yyaxis(ax, 'left');

    % --- six parameters: rotations share one y range, translations another ---
    t2 = tiledlayout(f, 1, 6, 'TileSpacing', 'compact', 'Padding', 'tight', 'OuterPosition', rowPos{2});
    pn = {'roll (°)','pitch (°)','yaw (°)','dS (mm)','dL (mm)','dP (mm)'};
    grp = {1:3, 4:6};  ylG = zeros(2,2);
    for g = 1:2
        % y range from the sweep only: the z-shim block's per-frame estimates are clipped (arrows)
        v = [reshape(pFull(sw,grp{g}),[],1); reshape(pInt(:,grp{g}),[],1)];
        pad = 0.08*max(range(v), 0.05);  ylG(g,:) = [min(v)-pad max(v)+pad];
    end
    for p = 1:6
        yl = ylG(1 + (p > 3), :);
        ax = nexttile(t2);  hold(ax, 'on');
        zband(ax, zs, yl, C);
        plot(ax, fr, pFull(:,p), '.', 'Color', C.full, 'MarkerSize', 8);
        hi = zs(pFull(zs,p) > yl(2));  lo = zs(pFull(zs,p) < yl(1));
        plot(ax, fr(hi), repmat(yl(2), size(hi)), '^', 'Color', C.full, 'MarkerSize', 4, 'MarkerFaceColor', C.full);
        plot(ax, fr(lo), repmat(yl(1), size(lo)), 'v', 'Color', C.full, 'MarkerSize', 4, 'MarkerFaceColor', C.full);
        plot(ax, fr, pInt(:,p), '-', 'Color', C.int, 'LineWidth', 2);
        plot(ax, fr(an), pInt(an,p), 'o', 'MarkerFaceColor', C.int, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
        xlim(ax, [-1 nF]);  ylim(ax, yl);  title(ax, pn{p}, 'FontWeight', 'normal', 'Color', C.ink);
        xlabel(ax, 'frame');  styleAx(ax, C);
    end

    % --- displacement map at the worst frame: three planes through the brain centroid ---
    tB = tiledlayout(f, 1, 6, 'TileSpacing', 'compact', 'Padding', 'tight', 'OuterPosition', rowPos{3});
    t3 = tiledlayout(tB, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
    t3.Layout.Tile = 1;  t3.Layout.TileSpan = [1 4];
    cmapD = interp1(linspace(0,1,6), [40 11 84; 101 21 110; 159 42 99; 212 72 66; 245 125 21; 250 193 39]/255, ...
        linspace(0,1,256));
    cmax = max(dMap(:), [], 'omitnan');
    [P, side] = threePlanes(dMap, base, mask, A);
    names = {'axial', 'coronal', 'sagittal'};
    gMax = prctile(base(mask), 99);
    for k = 1:3
        ax = nexttile(t3);  pl = P(k);
        g = min(max(pl.base / gMax, 0), 1);
        rgb = repmat(g, 1, 1, 3) * 0.9;
        idx = max(1, min(256, round(pl.d / cmax * 255) + 1));
        cm = reshape(cmapD(idx(:), :), [size(idx) 3]);
        m = repmat(pl.mask, 1, 1, 3);  rgb(m) = cm(m);
        image(ax, 'XData', pl.x, 'YData', pl.y, 'CData', rgb);  hold(ax, 'on');
        % equal square field of view (mm) for all three planes, centred on the slice
        xc = pl.cen(1);  yc = pl.cen(2);
        set(ax, 'XLim', xc + [-1 1]*side/2, 'YLim', yc + [-1 1]*side/2, 'YDir', 'reverse', ...
            'Color', 'k', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'DataAspectRatio', [1 1 1]);
        % crosshair: where the other two planes cut this one
        xline(ax, pl.cx, '-', 'Color', [1 1 1 0.6], 'LineWidth', 0.75);
        yline(ax, pl.cy, '-', 'Color', [1 1 1 0.6], 'LineWidth', 0.75);
        text(ax, 0.03, 0.04, pl.lab{1}, 'Units', 'normalized', 'Color', 'w', 'FontSize', 11, 'VerticalAlignment', 'bottom');
        text(ax, 0.97, 0.04, pl.lab{2}, 'Units', 'normalized', 'Color', 'w', 'FontSize', 11, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        text(ax, 0.03, 0.96, names{k}, 'Units', 'normalized', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
    end
    colormap(ax, cmapD);  clim(ax, [0 cmax]);
    cb = colorbar(ax, 'eastoutside');  cb.Label.String = 'displacement (voxels)';  cb.Color = C.ink2;
    title(t3, sprintf('Displacement at the largest-motion frame (%d, dTE %g ms)', R.worstFrame, dTE(R.worstFrame+1)), ...
        'FontSize', 11, 'Color', C.ink);

    % --- decision panel ---
    ax = nexttile(tB, 5, [1 2]);
    set(ax, 'Box', 'on', 'XTick', [], 'YTick', [], 'XColor', C.ink2, 'YColor', C.ink2, 'LineWidth', 1, ...
        'XLim', [0 1], 'YLim', [0 1]);
    switch R.decision
        case 'none',    head = 'NO CORRECTION';   sub = 'Motion is below the lower threshold.';
        case 'correct', head = 'CORRECT';         sub = 'Apply anchor-interpolated transforms to both echoes.';
        case 'exclude', head = 'EXCLUDE FRAMES';  sub = sprintf('%d frames are above the upper threshold.', numel(R.excludedFrames));
    end
    grey = '\color[rgb]{0.32,0.32,0.31}';
    L = {
        sprintf('\\bf%s\\rm', head), sub, '', ...
        sprintf('Max displacement:  \\bf%.2f voxels\\rm  (frame %d)', R.metric_vox, R.worstFrame), ...
        sprintf('Largest shift in mm:  %.2f mm  (voxel %.2f \\times %.2f \\times %.2f mm)', R.metric_mm, R.voxel_mm), ...
        sprintf('Thresholds:  %.2g / %.2g voxels  (above: %s)', R.threshLow, R.threshHigh, R.highAction), '', ...
        sprintf('Anchors:  dTE \\leq %g ms,  %d frames', R.anchorMaxDTE, numel(R.anchors)), ...
        sprintf('Largest gap between anchors:  %d frames%s', R.maxGap, tsec(R.maxGap, R.TR)), ...
        sprintf('Brain mask:  %s', strrep(maskLabel(R.mask), '_', '\_')), '', ...
        [grey 'Every frame registered (current pipeline):'], ...
        sprintf('%s    sweep max %.2f vox, z-shim max %.2f vox', grey, R.fullSweepMax_vox, R.fullZshimMax_vox), ...
        [grey 'Z-shim values off the parameter scales are shown as arrows.']};
    text(ax, 0.04, 0.95, L, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10.5, 'Color', C.ink, ...
        'Interpreter', 'tex');

    annotation(f, 'textbox', [0 0.965 1 0.033], 'String', sprintf('Motion decision: %s', R.subject), ...
        'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 15, 'Color', C.ink);
    exportgraphics(f, file, 'Resolution', 200);  close(f);
    fprintf('wrote %s\n', file);
end

function s = maskLabel(f)
    [~, n, e] = fileparts(f);  s = [n e];
    if strcmp(s, 'automask.nii.gz'), s = '3dAutomask of ase_base (no pipeline mask)'; end
end

function s = tsec(n, TR)
    if isnan(TR), s = ''; else, s = sprintf(' (%.0f s)', n*TR); end
end

function zband(ax, zs, yl, C)
    patch(ax, [zs(1)-1.5 zs(end)-0.5 zs(end)-0.5 zs(1)-1.5], [yl(1) yl(1) yl(2) yl(2)], C.shade, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function styleAx(ax, C)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'XColor', C.ink2, 'YColor', C.ink2, 'FontSize', 10, ...
        'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1);
end

function [P, side] = threePlanes(dMap, base, mask, A)
    % axial, coronal and sagittal slices through the mask centroid, with mm coordinates (so the
    % 4.8 mm slices are drawn at the right aspect) and the crosshair position of the other planes
    [~, ax] = max(abs(A(1:3,1:3)), [], 2);
    vs = sqrt(sum(A(1:3,1:3).^2));  vs = vs([ax(2) ax(1) ax(3)]);   % display dims (A->P, R->L, I->S)
    [dMap, base, mask] = deal(toDisplay(dMap, A), toDisplay(base, A), toDisplay(mask, A));
    dMap(~mask) = 0;
    [r, c, s] = ind2sub(size(mask), find(mask));  ctr = round(mean([r c s]));
    n = size(mask);  sI = n(3):-1:1;                                  % superior at the top
    yS = @(k) (n(3) - k + 1) * vs(3);                                 % display y (mm) of slice k
    sl = {@(V) V(:, :, ctr(3)), ...                                   % axial: A at top, R on left
          @(V) squeeze(V(ctr(1), :, sI))', ...                        % coronal: S at top, R on left
          @(V) squeeze(V(:, ctr(2), sI))'};                           % sagittal: S at top, A on left
    xs = {[1 n(2)]*vs(2), [1 n(2)]*vs(2), [1 n(1)]*vs(1)};
    ys = {[1 n(1)]*vs(1), [1 n(3)]*vs(3), [1 n(3)]*vs(3)};
    cx = [ctr(2)*vs(2), ctr(2)*vs(2), ctr(1)*vs(1)];                  % sagittal / sagittal / coronal cut
    cy = [ctr(1)*vs(1), yS(ctr(3)), yS(ctr(3))];                      % coronal / axial / axial cut
    lab = {{'R','L'}, {'R','L'}, {'A','P'}};
    % field of view: the head's largest extent over the three display axes, plus a margin
    ext = [range(r)*vs(1), range(c)*vs(2), range(s)*vs(3)];  side = 1.2 * max(ext);
    cen = {[mean(c)*vs(2), mean(r)*vs(1)], [mean(c)*vs(2), yS(mean(s))], [mean(r)*vs(1), yS(mean(s))]};
    for k = 3:-1:1
        P(k) = struct('d', sl{k}(dMap), 'base', sl{k}(base), 'mask', logical(sl{k}(mask)), ...
            'x', xs{k}, 'y', ys{k}, 'cen', cen{k}, 'cx', cx(k), 'cy', cy(k), 'lab', {lab{k}});
    end
end

function V = toDisplay(V, A)
    % permute/flip so dims are (anterior->posterior, right->left, inferior->superior)
    [~, ax] = max(abs(A(1:3,1:3)), [], 2);    % world axis r -> voxel dim ax(r)
    V = permute(V, [ax(2) ax(1) ax(3)]);
    if A(2, ax(2)) > 0, V = flip(V, 1); end    % anterior at the top
    if A(1, ax(1)) > 0, V = flip(V, 2); end    % patient right on the left
    if A(3, ax(3)) < 0, V = flip(V, 3); end
end

%% ===================== transforms =====================
function M = readMats(f)
    t = regexprep(fileread(f), '#[^\n]*', '');  v = sscanf(t, '%f');  M = reshape(v, 12, [])';
end

function Mi = interpMats(M, an, nF)
    % linear in acquisition order between anchors; hold before the first and after the last
    Mi = zeros(nF, size(M,2));  f = (1:nF)';
    for c = 1:size(M,2)
        if isscalar(an), Mi(:,c) = M(an,c); continue, end
        Mi(:,c) = interp1(an, M(an,c), f, 'linear');
        Mi(f < an(1), c) = M(an(1), c);  Mi(f > an(end), c) = M(an(end), c);
    end
end

function [d, dv] = maxDisp(M, xyz, toVox)
    % max over brain voxels of the voxel-normalised displacement, per matrix row
    d = zeros(size(M,1), 1);  X = [xyz ones(size(xyz,1),1)]';
    for k = 1:size(M,1)
        T = reshape(M(k,:), 4, 3)';  dv = sqrt(sum((toVox * (T*X - X(1:3,:))).^2))';  d(k) = max(dv);
    end
end

function P = matToParams(M)
    % small-angle rotations (deg) and shifts (mm) from DICOM base->source matrices, in 3dvolreg
    % order roll (about I-S), pitch (about R-L), yaw (about A-P), dS, dL, dP. For display only:
    % signs follow the matrix, which may be opposite to 3dvolreg's -dfile convention.
    P = zeros(size(M,1), 6);
    for k = 1:size(M,1)
        T = reshape(M(k,:), 4, 3)';  Rm = T(:,1:3);  S = (Rm - Rm')/2;
        P(k,:) = [rad2deg([S(2,1) S(3,2) S(1,3)]) T(3,4) T(1,4) T(2,4)];
    end
end

%% ===================== AFNI =====================
function s = q(p), s = ['"' p '"']; end

function afni(fmt, varargin)
    cmd = sprintf(fmt, varargin{:});
    [st, out] = system(['export PATH="$HOME/abin:$PATH"; export AFNI_NIFTI_TYPE_WARN=NO; ' cmd]);
    if st ~= 0, error('AFNI command failed:\n%s\n%s', cmd, out); end
end
