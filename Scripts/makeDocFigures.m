function makeDocFigures(hdrFile, outDir, opts)
% makeDocFigures  Figures for docs/pipeline.md.
%
%   makeDocFigures(hdrFile, outDir)
%     hdrFile : a subject's ase_hdr.json (frame timing is read from it)
%     outDir  : output folder (e.g. docs/figures)
%
% Writes frame_layout.png: every acquired frame on the tau axis, coloured by how the default fit
% uses it (short regime / excluded / long regime / z-shim block), with t_crit and the buffer.
% Needs Aux_Mat_Functions on the path (readASEHeader).
    arguments
        hdrFile (1,:) char
        outDir (1,:) char
        opts.tCrit_buffer (1,2) double = [3 12]
        opts.StO2_max (1,1) double = 0.6
        opts.Hct (1,1) double = 0.42*0.85
        opts.dChi0 (1,1) double = 0.18e-6
    end
    if ~isfolder(outDir), mkdir(outDir); end
    [dTEs, zShim, TEs, info] = readASEHeader(hdrFile);
    gamma = 2.67502e8;
    tc = 1000/((4/3)*pi*gamma*opts.dChi0*opts.Hct*(1-opts.StO2_max)*info.B0);   % ms
    tcrit = round(1.5*tc, 3);
    lo = tcrit - opts.tCrit_buffer(1);  hi = tcrit + opts.tCrit_buffer(2);

    % colours (validated categorical slots) and ink
    cShort = [42 120 214]/255;  cLong = [235 104 52]/255;  cZ = [98 80 214]/255;
    cExcl = [150 149 145]/255;  ink = [11 11 11]/255;  ink2 = [82 81 78]/255;
    bandShort = [0.93 0.95 0.99];  bandExcl = [0.95 0.95 0.94];  bandLong = [0.99 0.94 0.92];

    f = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1600 660]);
    ax = axes(f,'Position',[0.135 0.16 0.62 0.72]); hold(ax,'on');
    xl = [-2 80];  yl = [0.3 3.9];
    patch(ax,[xl(1) lo lo xl(1)], yl([1 1 2 2]), bandShort, 'EdgeColor','none');
    patch(ax,[lo hi hi lo], yl([1 1 2 2]), bandExcl, 'EdgeColor','none');
    patch(ax,[hi xl(2) xl(2) hi], yl([1 1 2 2]), bandLong, 'EdgeColor','none');

    sweep = info.repInds(1,:);                      % one rep (reps share the same taus)
    rowY = [3 2];  hS = []; hL = []; hE = []; hZ = [];
    for e = 1:numel(TEs)
        tau = dTEs(sweep) + (TEs(e) - TEs(1));
        isS = tau <= lo;  isL = tau > hi;  isE = ~isS & ~isL;
        y = rowY(e)*ones(size(tau));
        if any(isS), hS = plot(ax, tau(isS), y(isS), 'o', 'MarkerSize',9, 'MarkerFaceColor',cShort, 'MarkerEdgeColor','w', 'LineWidth',1); end
        if any(isE), hE = plot(ax, tau(isE), y(isE), 'o', 'MarkerSize',9, 'MarkerFaceColor','w', 'MarkerEdgeColor',cExcl, 'LineWidth',1.5); end
        if any(isL), hL = plot(ax, tau(isL), y(isL), 'o', 'MarkerSize',9, 'MarkerFaceColor',cLong, 'MarkerEdgeColor','w', 'LineWidth',1); end
        % z-shim block: trailing frames, all at the echo's longest tau, spread vertically
        zi = info.zShimInds;  tz = dTEs(zi(1)) + (TEs(e) - TEs(1));
        yz = 1 + linspace(-0.32, 0.32, numel(zi)) + (e-1)*0;      % row 1
        hZ = plot(ax, tz*ones(size(zi)), yz, 's', 'MarkerSize',8, 'MarkerFaceColor',cZ, 'MarkerEdgeColor','w');
        text(ax, tz, 1.48, sprintf('echo %d', e), 'HorizontalAlignment','center', 'FontSize',12, 'Color',ink2);
    end

    xline(ax, tcrit, '--', 'Color', ink, 'LineWidth', 1.5);
    text(ax, tcrit + 0.6, 3.74, sprintf('t_{crit} = 1.5 t_c = %.1f ms', tcrit), 'HorizontalAlignment','left', 'FontSize',13, 'Color',ink);
    ya = 3.42;
    annotation_arrow(ax, tcrit, lo, ya, ink2);  text(ax, lo - 0.6, ya, sprintf('−%g ms', opts.tCrit_buffer(1)), 'HorizontalAlignment','right','FontSize',12,'Color',ink2);
    annotation_arrow(ax, tcrit, hi, ya, ink2);  text(ax, (tcrit+hi)/2 + 0.6, ya+0.14, sprintf('+%g ms', opts.tCrit_buffer(2)), 'HorizontalAlignment','center','FontSize',12,'Color',ink2);
    text(ax, lo/2, 0.5, 'short regime', 'HorizontalAlignment','center','FontSize',13,'Color',cShort*0.8,'FontWeight','bold');
    text(ax, (tcrit+hi)/2, 0.5, 'excluded', 'HorizontalAlignment','center','FontSize',13,'Color',ink2,'FontWeight','bold');
    text(ax, (hi+xl(2))/2, 0.5, 'long regime', 'HorizontalAlignment','center','FontSize',13,'Color',cLong*0.85,'FontWeight','bold');
    text(ax, TEs(2)-TEs(1), 1.75, sprintf('echo 2 starts at τ = TE_2 − TE_1 = %.0f ms', TEs(2)-TEs(1)), 'FontSize',12, 'Color',ink2, 'HorizontalAlignment','left');

    ax.YTick = [1 2 3];  ax.YTickLabel = {'z-shim block', sprintf('echo 2 (TE %.1f ms)', TEs(2)), sprintf('echo 1 (TE %.1f ms)', TEs(1))};
    ax.FontSize = 13;  ax.XColor = ink2;  ax.YColor = ink2;  ax.TickLength = [0 0];
    ax.XGrid = 'on';  ax.GridColor = [0.85 0.85 0.85];  ax.GridAlpha = 1;  ax.Layer = 'top';
    xlim(ax, xl);  ylim(ax, yl);  box(ax,'off');
    xlabel(ax, 'τ, dephasing time relative to the spin echo (ms)', 'FontSize',14, 'Color',ink);
    title(ax, 'How the default fit uses each acquired frame', 'FontSize',16, 'FontWeight','normal', 'Color',ink);
    lg = legend(ax, [hS hE hL hZ], {'short-regime fit','excluded (transition)','long-regime fit','z-shim (Gz estimate only)'}, ...
        'Location','eastoutside', 'Box','off', 'FontSize',12);
    lg.Position(1) = 0.77;
    exportgraphics(f, fullfile(outDir,'frame_layout.png'), 'Resolution', 150, 'BackgroundColor','white');
    close(f);
    fprintf('t_c = %.2f ms, t_crit = %.2f ms, short <= %.2f, long > %.2f; %d z-shim frames per echo\n', ...
        tc, tcrit, lo, hi, numel(info.zShimInds));
    nz = nnz(zShim(info.zShimInds,1) ~= 0);
    fprintf('z-shim moments per echo: %d nonzero + %d zero-moment (all excluded from the regime fits by default)\n', ...
        nz, numel(info.zShimInds) - nz);
end

function annotation_arrow(ax, x0, x1, y, col)
    plot(ax, [x0 x1], [y y], '-', 'Color', col, 'LineWidth', 1.2);
    plot(ax, x1, y, '|', 'Color', col, 'MarkerSize', 10, 'LineWidth', 1.5);
end
