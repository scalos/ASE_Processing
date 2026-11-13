classdef ASE_Viewer<handle

    properties
        aseData = []
        aseInfo = []
        maps = struct;
        fitCoeffsLong = []
        fitCoeffsShort = []
        params = struct;
        frames = struct;   % per-frame info for aseData: taus, TEs, echo, use, shortInds, longInds
        roiM = [];         % sparse nReg x nAseVox fractional membership (same as the fit uses)
        roiNames = {}
        roiStats = table;
        window;
        windowName = "ASE_Viewer"
        viewParams = struct;
        uiObjs = struct;
    end

    methods
        function obj = ASE_Viewer(manifestFile)
            arguments
                manifestFile = []
            end
            if ~isempty(manifestFile)
                obj.loadManifest(manifestFile);
            end
            
        end

        function loadManifest(obj,manifestFile,opts)
            arguments
                obj
                manifestFile
                opts.verbose = true;
            end
            v = opts.verbose;
            assert(isfile(manifestFile),'ERROR: Manifest file does not exist %s',manifestFile)
            try
                manifest = jsondecode(fileread(manifestFile));
            catch ME
                fprintf('ERROR: Unable to read manifest file!\n');
                rethrow(ME);
            end
            fitOpts = struct;
            if isfield(manifest,'fitOpts'); fitOpts = manifest.fitOpts; end

            %% ASE data: rebuild exactly as ASE_fitting combined it
            aseFiles = cellstr(manifest.sorted_ASE_Files);
            if v
                fprintf('Loading ase data from:\n');
                fprintf('\t%s\n',aseFiles{:});
            end
            [~,~,~,hdrInfo] = readASEHeader(manifest.hdr_File);
            [~,obj.aseInfo] = loadNifti(aseFiles{1});
            avg = getOr(manifest,'averageReps',getOr(fitOpts,'averageReps',false));
            if avg; method = 'avg'; else; method = 'stack'; end
            obj.aseData = single(combineAseReps(aseFiles,hdrInfo,"method",method));

            %% per-frame info
            nt = size(obj.aseData,4);
            f = struct;
            f.taus = manifest.taus_full(:);
            f.TEs  = manifest.TEs_full(:);
            f.echo = manifest.frameEcho(:);
            f.use  = logical(manifest.usePts_full(:));
            f.shortInds = manifest.shortInds_good(:);
            f.longInds  = manifest.longInds_good(:);
            f.nEchos = numel(aseFiles);
            assert(numel(f.taus)==nt, ...
                'Manifest lists %d frames but the combined ASE data has %d.',numel(f.taus),nt);
            obj.frames = f;

            %% parameters
            obj.params.B0 = manifest.B0;
            obj.params.critTime = manifest.t_crit;
            obj.params.tCrit_buffer = getOr(fitOpts,'tCrit_buffer',0);
            obj.params.tauRange = getOr(fitOpts,'tauRange',[]);
            obj.params.normErr = logical(getOr(fitOpts,'normErr',true));   % imErr = RMS / mean signal

            %% maps (paths keyed by map name in manifest.paths)
            paths = struct;
            if isfield(manifest,'paths') && isstruct(manifest.paths); paths = manifest.paths; end
            obj.maps = struct;
            defs = ASE_Viewer.mapDefs();
            mapFields = defs(:,1).';            % every displayable map present in the manifest
            for idx = 1:numel(mapFields)
                mField = mapFields{idx};
                if isfield(paths,mField)
                    if v && isfield(paths,mField);fprintf('Loading map %s from %s...\n',mField,paths.(mField));end
                    m = loadNifti(paths.(mField),"niiName",mField,'failType','warn');
                    if ~isempty(m); obj.maps.(mField) = m; end
                end
            end
            if isfield(paths,"CoeffsLong")
               if v;fprintf('Loading long regime fitting coeffs from %s...\n',paths.CoeffsLong);end
               obj.fitCoeffsLong = loadNifti(paths.CoeffsLong,"niiName",'Coeff Mat Long','failType','warn');
            end
            if isfield(paths,"CoeffsShort")
               if v;fprintf('Loading short regime fitting coeffs from %s...\n',paths.CoeffsShort);end
               obj.fitCoeffsShort = loadNifti(paths.CoeffsShort,"niiName",'Coeff Mat Short','failType','warn');
            end
            obj.roiStats = table;
            if isfield(paths,"ROI_Stats")
                if v;fprintf('Loading ROI stats from %s...\n',paths.ROI_Stats);end
                try
                    obj.roiStats = readtable(paths.ROI_Stats,'ReadRowNames',true,'VariableNamingRule','preserve');
                catch ME
                    warning('Unable to load roi stats.\n@: %s\nME: (%s; %s)', ...
                        paths.ROI_Stats,ME.identifier,ME.message);
                end
            end

            %% ROIs: fractional membership on the ASE grid, same as the fit
            obj.roiM = [];  obj.roiNames = {};
            segFiles = getOr(fitOpts,'segFiles',{});
            segLUTs  = getOr(fitOpts,'segLUTs',{});
            if ~isempty(segFiles) && ~isempty(segLUTs)
                segFiles = cellstr(segFiles);  segLUTs = cellstr(segLUTs);
                if v; fprintf('Building ROI membership from %d segmentation(s)...\n',numel(segFiles)); end
                orig_state = warning;
                if ~getOr(fitOpts,'warnSegIssues',true), warning('off','all'); end
                [segVol,LUT] = combineSegVols(segFiles,segLUTs,"verbose",false,"keepOverlaps",1);
                warning(orig_state);
                segInfo = niftiinfo(segFiles{1});
                obj.roiM = segToAseMembership(segVol,LUT,segInfo,obj.aseInfo, ...
                    "maxSampling",[5,5,5],"minFrac",getOr(fitOpts,'minOverlapFrac',0));
                obj.roiNames = cellstr(LUT.name);
            end

            function val = getOr(st,fld,default)
                if isstruct(st) && isfield(st,fld) && ~isempty(st.(fld)); val = st.(fld); else; val = default; end
            end
        end

        function newFocPoint(obj,src,~)
            if ~strcmp(obj.uiObjs.regionLB.Value,'Crosshair'); return; end
            cfg = obj.viewParams.viewCfg;                          % config for the current view
            sz  = size(obj.aseData);
        
            p = src.CurrentPoint;
            h = min(max(round(p(1,1)),1), sz(cfg.h));   % index along raw dim cfg.h
            v = min(max(round(p(1,2)),1), sz(cfg.v));   % index along raw dim cfg.v
        
            obj.viewParams.currFoc(cfg.h) = h;
            obj.viewParams.currFoc(cfg.v) = v;                             % cfg.f is left unchanged
        
            obj.syncSpinners;                           % write idx back to the 3 spinners
            obj.updateWindow;
        end
        

        function initWindow(obj)
            if ~isempty(obj.window)
                if isvalid(obj.window)
                    delete(obj.window);
                end
            end
            obj.window = uifigure;
            obj.window.Position = [0 61 755 860];
            mainGrid = uigridlayout(obj.window,[3,2],"ColumnWidth",{'2.1x','7x','2x'},"RowHeight",{'3x','2x'});
            gridNW = uigridlayout(mainGrid,[4,1],'ColumnWidth',{'1x'},'RowHeight',{'2x','2x','2x','0.5x','2x'});
            gridNW.Layout.Row = 1;
            gridNW.Layout.Column = 1;
            gridNE = uigridlayout(mainGrid,[6,2],'ColumnWidth',{'1x','1x'},'RowHeight',{'1x','3x','1x','3x','1x','3x'});
            gridNE.Layout.Row = 1;
            gridNE.Layout.Column = 3;
            upperAx = axes(mainGrid);
            upperAx.Layout.Row = 1;
            upperAx.Layout.Column = 2;
            lowerAx = axes(mainGrid);
            lowerAx.Layout.Row = 2;
            lowerAx.Layout.Column =[1,3];

            
            
        
            % Build the toolbar with the native tools you want available in "explore" mode
            

            
            aseImg = imagesc(upperAx,[]);
            axis(upperAx,[0,size(obj.aseData,2)+1,0,size(obj.aseData,1)+1]);
            axis(upperAx,'off');
            set(upperAx,'YDir','normal');
            
            upperAx.PickableParts = 'all';
            set(aseImg, 'HitTest', 'off');

            colorbar(upperAx);

            obj.uiObjs.volImg = aseImg;

            startSlice = floor(size(obj.aseData,3)/2);
            
            initCoords = [floor(size(obj.aseData,1)/2),floor(size(obj.aseData,2)/2)];
            
            
            set(upperAx, 'ButtonDownFcn', @(src,evt) obj.newFocPoint(src,evt));
            axtoolbar(upperAx, {'zoomin','zoomout','pan','restoreview'});
            obj.uiObjs.upperAx = upperAx;
            obj.uiObjs.lowerAx = lowerAx;

            obj.viewParams.aseCurrSlice = startSlice;
            obj.viewParams.aseCurrdTE = 1;
            obj.viewParams.aseFocCoords = initCoords;
            obj.viewParams.showEx = 1;
            obj.viewParams.resetLims = true;
            obj.viewParams.selectedROI = [];
            views = viewConfigFromAffine(niftiAffine(obj.aseInfo.raw),"radio");
            views.axial.t = 4;
            views.sagittal.t = 4;
            views.coronal.t = 4;
            obj.viewParams.views = views;
            obj.viewParams.spinTimer = timer('ExecutionMode','singleShot','StartDelay',0.05,'TimerFcn',@(src,evt) obj.updateWindow);
            obj.viewParams.viewCfg = views.axial;
            
            

            obj.uiObjs.roiHandle = gobjects(0);
            obj.uiObjs.crosshair = gobjects(0);


            volumeDDLbl = uilabel(gridNW,'Text','Volumes','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            volumeDDLbl.Layout.Row = 1;
            volumeDDLbl.Layout.Column = 1;
            obj.uiObjs.volumeDDLbl = volumeDDLbl;
            
            volumeDD = uidropdown(gridNW,'Value','ASE','Items',{'ASE'},'FontSize',18);
            volumeDD.Layout.Row = 2;
            volumeDD.Layout.Column = 1;
            volumeDD.ValueChangedFcn = @(src,evt) obj.updateWindow();
            addStyle(volumeDD,uistyle("Interpreter", "latex"));
            obj.uiObjs.volumeDD = volumeDD;

            defs = obj.mapDefs();
            volumes = {'ASE'};
            for idx = 1:size(defs,1)
                if isfield(obj.maps,defs{idx,1})
                    volumes{end+1} = defs{idx,2}; %#ok
                end
            end
            obj.uiObjs.volumeDD.Items = volumes;
        
            regionSelectionLbl = uilabel(gridNW,'Text',sprintf('Region\nSelection'),'FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            regionSelectionLbl.Layout.Row = 3;
            regionSelectionLbl.Layout.Column = 1;
            obj.uiObjs.regionSelectionLbl = regionSelectionLbl;

            regionEF = uieditfield(gridNW);
            regionEF.ValueChangingFcn = @(src,evt) obj.filterRegions(src,evt);
            regionEF.Layout.Row = 4;
            regionEF.Layout.Column = 1;
            obj.uiObjs.regionEF = regionEF;

            regionLB = uilistbox(gridNW,'Value','Crosshair','Items',{'Crosshair'},'FontSize',10);
            regionLB.ValueChangedFcn = @(src,evt) obj.updateWindow();
            regionLB.Layout.Row = 5;
            regionLB.Layout.Column = 1;
            addStyle(regionLB,uistyle("Interpreter", "latex"));
            

            regionOptions = {'Crosshair'};
            if ~isempty(obj.roiNames)
                nROIs = numel(obj.roiNames);
                regionOptions = vertcat(regionOptions,cell(nROIs,1));
                for rIdx = 1:nROIs
                    regionOptions{rIdx+1} = obj.roiNames{rIdx};
                end
            end
            obj.viewParams.regionOptions = regionOptions;
            regionLB.Items = regionOptions;
            obj.uiObjs.regionLB = regionLB;

            sliceSpinnerLbl = uilabel(gridNE,'Text','Slice','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            sliceSpinnerLbl.Layout.Row = 1;
            sliceSpinnerLbl.Layout.Column = 1;
            obj.uiObjs.SliceSpinnerLbl = sliceSpinnerLbl;

            sliceSpinner = uispinner(gridNE,'Tag','f');
            sliceSpinner.Layout.Row = 2;
            sliceSpinner.Layout.Column = 1;
            sliceSpinner.ValueChangedFcn = @(src,evt) obj.spinnerChg(src,evt);
            obj.uiObjs.sliceSpinner = sliceSpinner;
        
            dTESpinnerLbl = uilabel(gridNE,'Text','Frame','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            dTESpinnerLbl.Layout.Row = 1;
            dTESpinnerLbl.Layout.Column = 2;
            obj.uiObjs.dTESpinnerLbl = dTESpinnerLbl;

            dTESpinner = uispinner(gridNE,'Tag','t');
            dTESpinner.Layout.Row = 2;
            dTESpinner.Layout.Column = 2;
            dTESpinner.ValueChangedFcn = @(src,evt) obj.spinnerChg(src,evt);
            obj.uiObjs.dTESpinner = dTESpinner;

            dim1SpinnerLbl = uilabel(gridNE,'Text',"←→",'FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            dim1SpinnerLbl.Layout.Row = 3;
            dim1SpinnerLbl.Layout.Column = 1;
            obj.uiObjs.dim1SpinnerLbl = dim1SpinnerLbl;

            dim1Spinner = uispinner(gridNE,'Tag','h');
            dim1Spinner.Layout.Row = 4;
            dim1Spinner.Layout.Column = 1;
            dim1Spinner.ValueChangedFcn = @(src,evt) obj.spinnerChg(src,evt);
            obj.uiObjs.dim1Spinner = dim1Spinner;
            
            dim2SpinnerLbl = uilabel(gridNE,'Text',"↑↓",'FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            dim2SpinnerLbl.Layout.Row = 3;
            dim2SpinnerLbl.Layout.Column = 2;
            obj.uiObjs.dim2SpinnerLbl = dim2SpinnerLbl;
        
            dim2Spinner = uispinner(gridNE,'Tag','v');
            dim2Spinner.Layout.Row = 4;
            dim2Spinner.Layout.Column = 2;
            dim2Spinner.ValueChangedFcn = @(src,evt) obj.spinnerChg(src,evt);
            obj.uiObjs.dim2Spinner = dim2Spinner;

            viewLBlbl = uilabel(gridNE,'Text','View','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
            viewLBlbl.Layout.Row = 5;
            viewLBlbl.Layout.Column = [1,2];
            obj.uiObjs.viewLBlbl = viewLBlbl;

            viewLB = uilistbox(gridNE,"Value",'A', ...
                                      'Items',{'Axial','Sagittal','Coronal'}, ...
                                      'ItemsData',{'A','S','C'});
            viewLB.Layout.Row = 6;
            viewLB.Layout.Column = [1,2];
            viewLB.ValueChangedFcn = @(src,evt) obj.viewSelectionChg(src,evt);
            obj.uiObjs.viewLB = viewLB;

            %init nav spinners
            obj.viewParams.currFoc = [floor(size(obj.aseData,1)/2),...
                                      floor(size(obj.aseData,2)/2),...
                                      floor(size(obj.aseData,3)/2),...
                                      1];
            obj.syncSpinners;
            obj.updateWindow();
        end

        function viewSelectionChg(obj,src,~)
            switch src.Value
                case 'A'
                    obj.viewParams.viewCfg = obj.viewParams.views.axial;
                case 'S'
                    obj.viewParams.viewCfg = obj.viewParams.views.sagittal;
                case 'C'
                    obj.viewParams.viewCfg = obj.viewParams.views.coronal;

            end
            obj.viewParams.resetLims = true;   % new view has a different in-plane size
            obj.syncSpinners;                  % spinners now map to different raw dims
            obj.updateWindow;
        end

        function filterRegions(obj,~,evt)
            currTxt = evt.Value;
            rOpts = obj.viewParams.regionOptions;
            msk = contains(rOpts,currTxt,'IgnoreCase',true);
            obj.uiObjs.regionLB.Items = rOpts(msk);
        end

        function spinnerChg(obj,src,~)
            cfg = obj.viewParams.viewCfg;
            currFoc = obj.viewParams.currFoc;
            spinnerDir = src.Tag;
            if ~isfield(cfg,spinnerDir)
                error('Tag %s doesn''t match cfg orientation (v,h,f,t)',spinnerDir);
            end
            
            currFoc(cfg.(spinnerDir)) = src.Value;
            obj.viewParams.currFoc = currFoc;
            if strcmp(obj.viewParams.spinTimer.Running,'on')
                stop(obj.viewParams.spinTimer);
            end
            start(obj.viewParams.spinTimer);
        end

        function syncSpinners(obj)
            cfg = obj.viewParams.viewCfg;
            sz = size(obj.aseData);  sz(end+1:4) = 1;
            currFoc = min(max(round(obj.viewParams.currFoc),1),sz(1:4));
            obj.viewParams.currFoc = currFoc;

            % Limits before Value: the new value may be outside the previous view's limits
            setSpin(obj.uiObjs.dim1Spinner,  currFoc(cfg.h), sz(cfg.h));
            setSpin(obj.uiObjs.dim2Spinner,  currFoc(cfg.v), sz(cfg.v));
            setSpin(obj.uiObjs.sliceSpinner, currFoc(cfg.f), sz(cfg.f));
            setSpin(obj.uiObjs.dTESpinner,   currFoc(4),     sz(4));

            function setSpin(sp,val,n)
                sp.Limits = [1,max(n,1+eps)];   % Limits must be increasing (n = 1 edge case)
                sp.Value = val;
            end
        end

        function toggleSpinners(obj,state)
            arguments
                obj
                state {mustBeMember(state,{'on','off'})};
            end
            obj.uiObjs.dim1Spinner.Enable = state;
            obj.uiObjs.dim2Spinner.Enable = state;
            obj.uiObjs.sliceSpinner.Enable = state;
            obj.uiObjs.dTESpinner.Enable = state;
        end


        function updateWindow(obj)
            %First make sure ui windows reflect currently loaded data:
            %%% Volume Drop Down:
            if ~isempty(obj.window)
                if ~isvalid(obj.window)
                    return;
                end
            else
                return;
            end
            obj.toggleSpinners("off");

            cfg = obj.viewParams.viewCfg;
            currFoc = obj.viewParams.currFoc;

            volAx = obj.uiObjs.upperAx;
            
            volSel = obj.uiObjs.volumeDD.Value;
            if strcmp(volSel,'ASE')
                f = obj.frames;
                title(volAx,sprintf('ASE frame %d (echo %d, $\\tau$ = %0.1f ms)', ...
                    currFoc(4),f.echo(currFoc(4)),f.taus(currFoc(4))),'Interpreter','latex','FontSize',16);
                imgData = obj.getSlice(obj.aseData(:,:,:,currFoc(4)),cfg,currFoc);
                set(volAx,'Colormap',gray(256));
                crossColor = [0.0667,0.4431,0.7451];
                clim(volAx,'auto');
            else
                defs = obj.mapDefs();
                d = defs(strcmp(defs(:,2),volSel),:);
                title(volAx,d{3},'Interpreter','latex','FontSize',20);
                imgData = obj.getSlice(obj.maps.(d{1}),cfg,currFoc).*d{4};
                set(volAx,'Colormap',d{5});
                if isempty(d{6}); clim(volAx,'auto'); else; clim(volAx,d{6}); end
                crossColor = d{7};
            end
            hold(volAx,'on');
            set(obj.uiObjs.volImg,'CData',imgData,'XData',[1,size(imgData,2)],'YData',[1,size(imgData,1)]);
            if ~isfield(obj.viewParams,'resetLims') || obj.viewParams.resetLims
                % fit axes to this view's slice (only on view change, so zoom is kept otherwise)
                xlim(volAx,[0.5,size(imgData,2)+0.5]);
                ylim(volAx,[0.5,size(imgData,1)+0.5]);
                obj.viewParams.resetLims = false;
            end
            set(volAx, 'XDir', cfg.xdir, 'YDir', cfg.ydir);   % set AFTER imagesc, which resets YDir
            pixdim = obj.aseInfo.PixelDimensions;
            daspect(volAx, [pixdim(cfg.v) pixdim(cfg.h) 1]);
            lbl = obj.edgeLabels(aff2axcodes(niftiAffine(obj.aseInfo.raw)),cfg);
            
            if isfield(obj.uiObjs,"codeLbls")
                codeLbls = obj.uiObjs.codeLbls;
                for idx = 1:numel(codeLbls)
                    delete(codeLbls{idx});
                end
            end
            lblColor = crossColor;
            lblSize = 12;
            lbl_l = text(volAx, 0.02, 0.5,  lbl.left,   'Units','normalized', ...
                'Color',lblColor, 'FontWeight','bold','FontSize',lblSize);

            lbl_r = text(volAx, 0.98, 0.5,  lbl.right,  'Units','normalized', ...
                'Color',lblColor, 'FontWeight','bold', 'HorizontalAlignment','right','FontSize',lblSize);

            lbl_t = text(volAx, 0.5,  0.97, lbl.top,    'Units','normalized', ...
                'Color',lblColor, 'FontWeight','bold', 'HorizontalAlignment','center','FontSize',lblSize);

            lbl_b = text(volAx, 0.5,  0.03, lbl.bottom, 'Units','normalized', ...
                'Color',lblColor, 'FontWeight','bold', 'HorizontalAlignment','center','FontSize',lblSize);

            obj.uiObjs.codeLbls = {lbl_l,lbl_r,lbl_t,lbl_b};


            f = obj.frames;
            roiData = [];
            roiW = [];
            roiEx = [];
            selectedRegion = obj.uiObjs.regionLB.Value;
            if ~strcmp(selectedRegion,'Crosshair')
                rIdx = find(strcmp(obj.roiNames,selectedRegion),1);
                w = obj.roiM(rIdx,:);
                % voxel quality mask (if the fit used one): excluded voxels are drawn in a
                % different colour and left out of the ROI curve, as in the fit
                roiEx = [];
                if isfield(obj.maps,'qualityMask')
                    keepV = obj.maps.qualityMask(:).' > 0;
                    exIdx = find(w & ~keepV);
                    roiEx = zeros(size(obj.aseData,1:3),'single');
                    roiEx(exIdx) = full(w(exIdx));
                    w = w .* keepV;
                end
                vIdx = find(w);
                roiW = zeros(size(obj.aseData,1:3),'single');
                roiW(vIdx) = full(w(vIdx));
                % weighted mean per frame ignoring invalid values (same as the fit's regionMeans)
                A = reshape(obj.aseData,[],size(obj.aseData,4));
                X = double(A(vIdx,:));
                ok = isfinite(X) & X > 0;
                X(~ok) = 0;
                wv = full(w(vIdx));
                roiData = (wv*X) ./ (wv*double(ok));
                roiData = roiData(:);
            end
            if isfield(obj.uiObjs,'roiExHandle') && ~isempty(obj.uiObjs.roiExHandle) && isvalid(obj.uiObjs.roiExHandle)
                delete(obj.uiObjs.roiExHandle);
            end
            obj.uiObjs.roiExHandle = gobjects(0);

            fullFoc = [obj.uiObjs.dim1Spinner.Value,...
                       obj.uiObjs.dim2Spinner.Value,...
                       obj.uiObjs.sliceSpinner.Value,...
                       obj.uiObjs.dTESpinner.Value];
            if ~isempty(obj.uiObjs.roiHandle) && isvalid(obj.uiObjs.roiHandle)
                delete(obj.uiObjs.roiHandle);
            end
            if isempty(roiW)
                redrawCross = false;
                crosshair = obj.uiObjs.crosshair;
                if isempty(obj.uiObjs.crosshair)
                    redrawCross = true;
                elseif ~isvalid(obj.uiObjs.crosshair)
                    redrawCross = true;
                end
                if redrawCross
                    crosshair = images.roi.Crosshair(obj.uiObjs.upperAx, ...
                                                    'LineWidth',1, ...
                                                    'InteractionsAllowed','none');
                end
                crosshair.Position = [fullFoc(1),fullFoc(2)];
                crosshair.Color = crossColor;
                obj.uiObjs.crosshair = crosshair;
            else
                if ~isempty(obj.uiObjs.crosshair)&&isvalid(obj.uiObjs.crosshair)
                    delete(obj.uiObjs.crosshair);
                end
                roi_slice = double(obj.getSlice(roiW,cfg,currFoc));
                roi_TC = cat(3,ones(size(roi_slice)),zeros(size(roi_slice)),zeros(size(roi_slice)));
                obj.uiObjs.roiHandle = image(obj.uiObjs.upperAx,roi_TC,'HitTest','off');
                set(obj.uiObjs.roiHandle,'AlphaData',roi_slice.*0.5);
                if ~isempty(roiEx)
                    ex_slice = double(obj.getSlice(roiEx,cfg,currFoc));
                    ex_TC = cat(3,zeros(size(ex_slice)),ones(size(ex_slice)),ones(size(ex_slice)));   % cyan
                    obj.uiObjs.roiExHandle = image(obj.uiObjs.upperAx,ex_TC,'HitTest','off');
                    set(obj.uiObjs.roiExHandle,'AlphaData',ex_slice.*0.6);
                end
            end
            hold(obj.uiObjs.upperAx,'off');

            dTE_ax = obj.uiObjs.lowerAx;
            cla(dTE_ax);
            hold(dTE_ax,'on');
            shortColor = [0.0667,0.4431,0.7451];
            longColor = [0.5216,0.0863,0.8196];
            exColor = [0.8667,0.3294,0];
            echoMarkers = {'o','s','^','d','v','p'};
            echoLines = {'--',':','-.','--',':','-.'};

            % fit coefficients: long = [a, (R2,) R2'], short = [b, (k)] in s / s^2 units
            if ~isempty(roiData)
                focData = log(roiData);
                [BL,BS,R2,R2P,vCBV,OEF] = obj.roiCoeffs(selectedRegion);
                Signal = roiData(currFoc(4));
                % fit RMS from the ROI fit; imErr = ROI-weighted mean of the voxel imErr map
                [rmsL,rmsS] = obj.roiRms(selectedRegion);
                imErrVal = NaN;
                if isfield(obj.maps,'imErr')
                    e = double(obj.maps.imErr(vIdx));  e = e(:).';  ok = isfinite(e);
                    if any(ok); imErrVal = sum(wv(ok).*e(ok)) / sum(wv(ok)); end
                end
                imErrLbl = 'mean imErr';
            else
                spatIdxCell = num2cell(currFoc(1:3));
                focData = log(double(squeeze(obj.aseData(spatIdxCell{:},:))));
                BL = []; BS = [];
                if ~isempty(obj.fitCoeffsLong);  BL = double(squeeze(obj.fitCoeffsLong(spatIdxCell{:},:))); end
                if ~isempty(obj.fitCoeffsShort); BS = double(squeeze(obj.fitCoeffsShort(spatIdxCell{:},:))); end
                R2   = obj.mapVal('R2',spatIdxCell);
                R2P  = obj.mapVal('R2P',spatIdxCell);
                vCBV = obj.mapVal('vCBV',spatIdxCell);
                OEF  = obj.mapVal('OEF',spatIdxCell);
                Signal = double(obj.aseData(spatIdxCell{:},currFoc(4)));
                rmsL = obj.mapVal('RMS_L',spatIdxCell);
                rmsS = obj.mapVal('RMS_S',spatIdxCell);
                imErrVal = obj.mapVal('imErr',spatIdxCell);
                imErrLbl = 'imErr';
            end
            focData = focData(:);
            focData(~isfinite(focData)) = NaN;

            nt = numel(f.taus);
            isShort = false(nt,1); isShort(f.shortInds) = true;
            isLong  = false(nt,1); isLong(f.longInds) = true;
            isShort = isShort & f.use;  isLong = isLong & f.use;
            isEx = ~(isShort | isLong);
            usedData = focData(isShort | isLong);

            if all(isnan(usedData))
                xlim(dTE_ax,[0,1]);
                ylim(dTE_ax,[0,1]);
                Signal = NaN; R2 = NaN; R2P = NaN; vCBV = NaN; OEF = NaN;
                rmsL = NaN; rmsS = NaN; imErrVal = NaN;
                text(dTE_ax,0.5,0.5,'No Valid Points!','FontSize',18,'HorizontalAlignment','center');
            else
                lgdHandles = gobjects(0);
                for e = 1:f.nEchos
                    mk = echoMarkers{mod(e-1,numel(echoMarkers))+1};
                    ls = echoLines{mod(e-1,numel(echoLines))+1};
                    inE = f.echo == e;
                    TEe = f.TEs(find(inE,1))/1000;   % s
                    sfx = '';
                    if f.nEchos > 1; sfx = sprintf(' (echo %d)',e); end
                    m = inE & isShort;
                    if any(m)
                        lgdHandles(end+1) = scatter(dTE_ax,f.taus(m),focData(m),36,shortColor,'filled', ...
                            'Marker',mk,'DisplayName',['Short Regime' sfx]); %#ok
                    end
                    m = inE & isLong;
                    if any(m)
                        lgdHandles(end+1) = scatter(dTE_ax,f.taus(m),focData(m),36,longColor,'filled', ...
                            'Marker',mk,'DisplayName',['Long Regime' sfx]); %#ok
                    end
                    m = inE & isEx;
                    if obj.viewParams.showEx && any(m)
                        lgdHandles(end+1) = scatter(dTE_ax,f.taus(m),focData(m),36,exColor, ...
                            'Marker','x','LineWidth',2,'DisplayName',['Excluded' sfx]); %#ok
                    end

                    % fit curves for this echo
                    [yL,yS] = obj.fitCurves(BL,BS,f.nEchos,TEe);
                    m = inE & isLong;
                    if ~isempty(yL) && any(m)
                        t = linspace(min(f.taus(m)),max(f.taus(m)),100);
                        lgdHandles(end+1) = plot(dTE_ax,t,yL(t),ls,'Color',longColor, ...
                            'DisplayName',['Long Regime Fit' sfx]); %#ok
                    end
                    m = inE & isShort;
                    if ~isempty(yS) && any(m)
                        t = linspace(min(f.taus(m)),max(f.taus(m)),100);
                        lgdHandles(end+1) = plot(dTE_ax,t,yS(t),ls,'Color',shortColor, ...
                            'DisplayName',['Short Regime Fit' sfx]); %#ok
                    end
                end

                yBds = [min(usedData),max(usedData)];
                yRange = max(abs(diff(yBds)),eps);
                padding = 0.1;
                critTime = obj.params.critTime;
                critBuffer = obj.params.tCrit_buffer;
                xlim(dTE_ax,[min(f.taus),max(f.taus)]);
                ylim(dTE_ax,[yBds(1)-yRange*padding,yBds(2)+yRange*padding]);
                lgdHandles(end+1) = xline(dTE_ax,critTime,'-','1.5x\ $t_{crit}$', ...
                                                 'Interpreter','latex', ...
                                                 'Color','black', ...
                                                 'FontSize',12,...
                                                 'DisplayName',sprintf('1.5x Critical Time (%0.1fms)',critTime));
                critRange = [max(0,critTime-critBuffer),critTime+critBuffer];
                lgdHandles(end+1) = xregion(dTE_ax,critRange, ...
                                                    'LineStyle','--', ...
                                                    'FaceColor',[1,0,0], ...
                                                    'EdgeColor',[0,0,0], ...
                                                    'EdgeAlpha',1, ...
                                                    'FaceAlpha',0.2,...
                                                    'DisplayName',sprintf('Buffer Window (%0.1fms)',2*critBuffer));
                lgdHandles(end+1) = xline(dTE_ax,f.taus(currFoc(4)),'-','$\tau_{current}$', ...
                                                   'Interpreter','latex', ...
                                                   'Color',[0.2314,0.6667,0.1961], ...
                                                   'FontSize',12,...
                                                   'DisplayName',sprintf('Current frame %d (echo %d, %0.1fms)', ...
                                                        currFoc(4),f.echo(currFoc(4)),f.taus(currFoc(4))));
                legend(dTE_ax,lgdHandles,'Interpreter','none','Location','eastoutside');
            end
            hold(dTE_ax,'off');
            xlabel(dTE_ax,'$\tau$ (dTE + TE$_e$ - TE$_1$) [ms]','Interpreter','latex');
            ylabel(dTE_ax,'log(S) [a.u.]','Interpreter','latex');

            dTE_ax.Title.FontName = 'FixedWidth';
            % two columns: map values | fit error values (lines padded to equal width so
            % the centred title keeps the columns aligned)
            if obj.params.normErr
                imErrStr = sprintf('%-6.2f [%%]',imErrVal*100);    % normalized: fraction of mean signal
            else
                imErrStr = sprintf('%-6.1f [a.u.]',imErrVal);
            end
            ttl = splitlines(sprintf([' Signal: %-6.1f [a.u.]    %10s: %-7.4f\n' ...
                                 '    R2: %-6.1f [Hz]      %10s: %-7.4f\n' ...
                                 '   R2'': %-6.1f [Hz]      %10s: %s\n' ...
                                 '  vCBV: %-6.1f [%%]\n' ...
                                 '   OEF: %-6.1f [%%]'], ...
                                 Signal,'RMS long',rmsL, ...
                                 R2,'RMS short',rmsS, ...
                                 R2P,imErrLbl,imErrStr, ...
                                 vCBV*100,OEF*100));
            ttl = char(ttl);
            ttl(ttl == ' ') = char(160);   % non-breaking spaces: trailing padding isn't trimmed
            title(dTE_ax,cellstr(ttl));
            dTE_ax.Title.Interpreter = 'none';
            obj.uiObjs.dim1Spinner.Value = fullFoc(1);
            obj.uiObjs.dim2Spinner.Value = fullFoc(2);
            obj.toggleSpinners("on");
        end
    end

    methods (Access = private)
        function val = mapVal(obj,name,spatIdxCell)
            if isfield(obj.maps,name); val = double(obj.maps.(name)(spatIdxCell{:})); else; val = NaN; end
        end

        function [rmsL,rmsS] = roiRms(obj,regionName)
            % ROI-fit RMS (log space) from the ROI stats table (NaN if unavailable)
            rmsL = NaN; rmsS = NaN;
            T = obj.roiStats;
            if isempty(T) || ~any(strcmp(T.Properties.RowNames,regionName)); return; end
            vn = T.Properties.VariableNames;
            if any(strcmp(vn,'RMS_L')); rmsL = T{regionName,'RMS_L'}; end
            if any(strcmp(vn,'RMS_S')); rmsS = T{regionName,'RMS_S'}; end
        end

        function [BL,BS,R2,R2P,vCBV,OEF] = roiCoeffs(obj,regionName)
            % coefficients/values for one region from the ROI stats table (NaN if unavailable)
            BL = []; BS = []; R2 = NaN; R2P = NaN; vCBV = NaN; OEF = NaN;
            T = obj.roiStats;
            if isempty(T) || ~any(strcmp(T.Properties.RowNames,regionName)); return; end
            row = T(regionName,:);
            vn = T.Properties.VariableNames;
            getv = @(n) row.(n);
            colsL = vn(startsWith(vn,'B_L_'));  colsS = vn(startsWith(vn,'B_S_'));
            BL = cellfun(getv,colsL);  BS = cellfun(getv,colsS);
            if any(strcmp(vn,'R2'));   R2   = row.R2;   end
            if any(strcmp(vn,'R2P'));  R2P  = row.R2P;  end
            if any(strcmp(vn,'vCBV')); vCBV = row.vCBV; end
            if any(strcmp(vn,'OEF'));  OEF  = row.OEF;  end
        end
    end

    methods (Static)
        function defs = mapDefs()
            % {field, dropdown label, title, display scale, colormap, clim ([] = auto), crosshair color}
            defs = {
                'R2P',  '$R_2''$',  'Computed $R_2''$ Volume [Hz]',      1,   parula(256), [0,15],  [0,0,0];
                'vCBV', 'vCBV',     'Computed vCBV Volume [\%]',          100, jet(256),    [0,15],  [1,1,1];
                'OEF',  'OEF',      'Computed OEF Volume [\%]',           100, rgb2cmap([0,0,1;1,1,1;1,0,0],256), [0,100], [0,0,0];
                'R2',   '$R_2$',    'Computed $R_2$ Volume [Hz]',         1,   parula(256), [0,30],  [0,0,0];
                'imErr','imErr',    'Fit Error',                          1,   hot(256),    [],      [0,0.6,1];
                'RMS_L','$RMS_l$',  'Long Regime Fit RMS (log space)',    1,   rgb2cmap([0,1,0;1,0,0],256), [0,0.1], [0,0,0];
                'RMS_S','$RMS_s$',  'Short Regime Fit RMS (log space)',   1,   rgb2cmap([0,1,0;1,0,0],256), [0,0.1], [0,0,0];
                % voxel-quality outputs (present only when quality features were enabled)
                'qualityMask', 'Quality mask', 'Voxel kept (1) / excluded (0)',          1, gray(256),   [0,1],    [1,0,0];
                'qualityFlags','Quality flags','Quality bit flags (see computeVoxelQuality)', 1, lines(256), [0,1023], [1,1,1];
                'vCBV_quad',   'vCBV (quad)',  'vCBV from short-regime curvature [\%]', 100, jet(256),    [0,15],   [1,1,1];
                'vCBV_z',      'vCBV z',       'z of vCBV$_{int}$ $-$ vCBV$_{quad}$',     1, rgb2cmap([0,0,1;1,1,1;1,0,0],256), [-5,5], [0,0,0];
                'R2P_e1',      '$R_2''$ e1',   'Echo-1 $R_2''$ [Hz]',                     1, parula(256), [0,15],   [0,0,0];
                'R2P_e2',      '$R_2''$ e2',   'Echo-2 $R_2''$ [Hz]',                     1, parula(256), [0,15],   [0,0,0];
                'R2P_Fstat',   '$R_2''$ F',    'Per-echo $R_2''$ F statistic',            1, hot(256),    [0,50],   [0,0.6,1];
                'Gz_e1',       '$G_z$ e1',     'Through-slice gradient, echo 1 [$\mu$T/m]', 1, rgb2cmap([0,0,1;1,1,1;1,0,0],256), [-100,100], [0,0,0];
                'Gz_e2',       '$G_z$ e2',     'Through-slice gradient, echo 2 [$\mu$T/m]', 1, rgb2cmap([0,0,1;1,1,1;1,0,0],256), [-100,100], [0,0,0];
                'zshimFracLow','Z-shim low F', 'Fraction of long frames with $F<F_{min}$', 1, hot(256),    [0,1],    [0,0.6,1];
                'sigma_e1',    '$\sigma$ e1', 'Repetition noise, echo 1',                1, hot(256),    [],       [0,0.6,1];
                'sigma_e2',    '$\sigma$ e2', 'Repetition noise, echo 2',                1, hot(256),    [],       [0,0.6,1];
                'SNR_e1',      'SNR e1',       'SNR (spin-echo signal / $\sigma$), echo 1', 1, parula(256), [0,200], [0,0,0];
                'SNR_e2',      'SNR e2',       'SNR (spin-echo signal / $\sigma$), echo 2', 1, parula(256), [0,200], [0,0,0];
                'csfFrac',     'CSF fraction', 'CSF fraction (1 $-$ tissue mask)',         1, gray(256),   [0,1],    [1,0,0];
            };
        end

        function [yL,yS] = fitCurves(BL,BS,nEchos,TEe)
            % log-signal model curves (tau in ms) for one echo, matching fitASERegimes:
            %   long : ln S = a - R2*TE - R2'*|tau|        (single echo: a - R2'*|tau|)
            %   short: ln S = b + k*tau^2 - R2*TE          (single echo: b + k*tau^2)
            yL = []; yS = [];
            if isempty(BL) || any(isnan(BL)); return; end
            if nEchos > 1
                R2 = BL(2); R2P = BL(3); te = R2*TEe;
            else
                R2P = BL(2); te = 0;
            end
            yL = @(t) BL(1) - te - R2P*abs(t)/1000;
            if isempty(BS) || isnan(BS(1)); return; end
            k = 0; if numel(BS) > 1 && ~isnan(BS(2)); k = BS(2); end
            yS = @(t) BS(1) + k*(t/1000).^2 - te;
        end

        function img = getSlice(vol3d, cfg, idx)
            % idx = [i j k] current raw voxel indices; returns matrix with rows = v, cols = h
            sub = {':',':',':'};
            sub{cfg.f} = idx(cfg.f);
            S = squeeze(vol3d(sub{:}));      % remaining dims are in ascending order
            if cfg.h < cfg.v
                img = S.';                      % rows=lower dim, so transpose to put h horizontal
            else
                img = S;                        % already rows=v (higher dim on columns = h)
            end
        end

        function lbl = edgeLabels(codes, cfg)
        % Returns a struct with left/right/top/bottom orientation letters for one view.
        %   codes : 1x3 char from aff2axcodes (direction each raw dim's index increases toward)
        %   cfg   : one entry from viewConfigFromAffine (needs h, v, xdir, ydir)
            hi  = codes(cfg.h);   hlo = opposite(hi);   % letters at high / low index, horizontal dim
            vi  = codes(cfg.v);   vlo = opposite(vi);   % letters at high / low index, vertical dim
        
            if strcmp(cfg.xdir, 'normal')
                lbl.right = hi;   lbl.left = hlo;       % index increases toward screen right
            else
                lbl.right = hlo;  lbl.left = hi;        % reversed axis: index increases toward screen left
            end
        
            if strcmp(cfg.ydir, 'normal')
                lbl.top = vi;     lbl.bottom = vlo;     % index increases toward screen top
            else
                lbl.top = vlo;    lbl.bottom = vi;      % reversed axis: index increases toward screen bottom
            end

            function o = opposite(c)
                from = 'LRAPSI';
                to   = 'RLPAIS';
                o = to(from == c);
            end
        end
    end
end