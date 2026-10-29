classdef ASE_Viewer<handle

    properties
        aseData = []
        aseInfo = []
        maps = struct;
        fitCoeffsLong = []
        fitCoeffsShort = []
        params = struct;
        rois = struct;
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
                manifestData_raw = fileread(manifestFile);
                manifest = jsondecode(manifestData_raw);
            catch ME
                fprintf('ERROR: Unable to read manifest file!');
                rethrow(ME);
            end
            aseFiles = manifest.sorted_ASE_File;
            if v
                fprintf('Loading ase data from:\n');
                for idx = 1:numel(aseFiles)
                     fprintf("%s\n",aseFiles{idx});
                end
            end
            [~,obj.aseInfo] = loadNifti(aseFiles{1}); 
            obj.aseData = combineAseReps(aseFiles);
            
            paramFields = {'B0','taus','tauRange','tCrit_buffer','tauPartitions','critTime'};
            if v;fprintf('Loading parameters...\n');end
            for idx = 1:numel(paramFields)
                pField = paramFields{idx};
                if v;fprintf('\t-%s\n',pField);end
                obj.params.(pField) = manifest.(pField);
            end
            
            %if computed, found in manifest.paths
            mapFields = {'R2P','vCBV','OEF'};
            for idx = 1:numel(mapFields)
                mField = mapFields{idx};
                if isfield(manifest.paths,mField)
                    if v;fprintf('Loading Map %s from %s...\n',mField,manifest.paths.(mField));end
                    obj.maps.(mField) = loadNifti(manifest.paths.(mField));
                end
                
            end
            if isfield(manifest.paths,"coeffsLong")
               if v;fprintf('Loading long regime fitting coeffs from %s...\n',manifest.paths.coeffsLong);end
               obj.fitCoeffsLong = loadNifti(manifest.paths.coeffsLong, ...
                   "niiName",'Coeff Mat Long','failType','warn'); 
               
            end
            if isfield(manifest.paths,"coeffsShort")
               if v;fprintf('Loading short regime fitting coeffs from %s...\n',manifest.paths.coeffsLong);end
               obj.fitCoeffsShort = loadNifti(manifest.paths.coeffsShort, ...
                   "niiName",'Coeff Mat Short','failType','warn');   
            end
            if isfield(manifest.paths,"statsCSV")
                if v;fprintf('Loading ROI stats from %s...\n',manifest.paths.statsCSV);end
                try
                    obj.roiStats = readtable(manifest.paths.statsCSV);
                catch ME
                    obj.roiStats = table;
                    warning('Unable to load roi stats.\n@: %s\nME: (%s; %s)', ...
                        manifest.paths.statsCSV,ME.identifier,ME.message);
                end
            end

            %load ROIs:
            if isfield(manifest,"ROI_Files")
                if ~isempty(manifest.ROI_Files)
                    %cell list was used to load ROI's
                    [obj.rois,obj.roiNames] = getROIsFromCell(manifest.ROI_Files,"verbose",v);
                end
            end
            if isfield(manifest,"segFiles") && isfield(manifest,"segLUTs")
                if ~isempty(manifest.segFiles) && ~isempty(manifest.segLUTs)
                    %segemented volume(s) was/were used to load ROI's
                    [obj.rois,obj.roiNames] = getROIsFromLUT(manifest.segFiles,manifest.segLUTs,"verbose",v);
                end
            end
            
        end

        function newFocPoint(obj,src,~)
            if ~strcmp(obj.uiObjs.regionLB.Value,'Crosshair'); return; end
            cfg = obj.viewParams.viewCfg;                          % config for the current view
            sz  = size(obj.aseData);
        
            p = get(src,'CurrentPoint');
            h = min(max(round(p(1,1)),1), sz(cfg.h));   % index along raw dim cfg.h
            v = min(max(round(p(1,2)),1), sz(cfg.v));   % index along raw dim cfg.v
        
            obj.viewParams.currFoc(cfg.h) = h;
            obj.viewParams.currFoc(cfg.v) = v;                             % cfg.f is left unchanged
        
            obj.syncSpinners;                           % write idx back to the 3 spinners
            obj.updateWindow;
            
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

            mapNames = struct();
            mapNames.R2P = '$R_2$''';
            mapNames.vCBV = 'vCBV';
            mapNames.OEF = 'OEF';
            fNames = fieldnames(mapNames);
            volumes = {'ASE'};
            for idx = 1:numel(fNames)
                mName = fNames{idx};
                if isfield(obj.maps,mName)
                    volumes{end+1} = mapNames.(mName); %#ok
                end
            end
            if ~isempty(obj.fitCoeffsLong)
                volumes{end+1} = '$RMS_l$';
            end
            if ~isempty(obj.fitCoeffsShort)
                volumes{end+1} = '$RMS_s$';
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
        
            dTESpinnerLbl = uilabel(gridNE,'Text','dTE','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
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

            viewLBLbl = uilabel(gridNE,'Text','View','FontSize',20,'HorizontalAlignment','center','Interpreter','latex','VerticalAlignment','bottom');
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
            sz = size(obj.aseData);
            currFoc = obj.viewParams.currFoc;

            obj.uiObjs.dim1Spinner.Value = currFoc(cfg.h);
            obj.uiObjs.dim1Spinner.Limits = [1,sz(cfg.h)];

            obj.uiObjs.dim2Spinner.Value = currFoc(cfg.v);
            obj.uiObjs.dim2Spinner.Limits = [1,sz(cfg.v)];

            obj.uiObjs.sliceSpinner.Value = currFoc(cfg.f);
            obj.uiObjs.sliceSpinner.Limits = [1,sz(cfg.f)];

            obj.uiObjs.dTESpinner.Value = currFoc(4);
            obj.uiObjs.dTESpinner.Limits = [1,sz(4)];
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

            crossColor = [0,0,0];
            volAx = obj.uiObjs.upperAx;
            
            switch obj.uiObjs.volumeDD.Value
                case 'ASE' 
                    title(volAx,'Post-Processed ASE Volume','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(obj.aseData(:,:,:,currFoc(4)),cfg,currFoc);
                    %imgData = squeeze(obj.aseData(:,:,currSlice,currdTE));
                    set(volAx,'Colormap',gray(256));
                    crossColor = [0.0667,0.4431,0.7451];
                    clim(volAx,'auto');
                case '$R_2$'''
                    title(obj.uiObjs.upperAx,'Computed $R_2^*$ Volume','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(obj.maps.R2P,cfg,currFoc);
                    set(obj.uiObjs.upperAx,'Colormap',parula(256));
                    clim(obj.uiObjs.upperAx,[0,15]);

                case 'vCBV'
                    title(volAx,'Computed vCBV Volume','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(obj.maps.vCBV,cfg,currFoc).*100;
                    set(volAx,'Colormap',jet(256));
                    clim(volAx,[0,15]);
                    crossColor = [1,1,1];
                case 'OEF'
                    title(volAx,'Computed OEF Volume','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(obj.maps.OEF,cfg,currFoc).*100;
                    set(volAx,'Colormap',rgb2cmap([0,0,1;1,1,1;1,0,0],256));
                    clim(volAx,[0,100]);
                case '$RMS_s$'
                    title(volAx,'Computed Short Regime RMS','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(squeeze(obj.fitCoeffsShort(:,:,:,3)),cfg,currFoc);
                    set(volAx,'Colormap',rgb2cmap([0,1,0;1,0,0],256));
                    clim(volAx,[0,0.1]);
                case '$RMS_l$'
                    title(volAx,'Computed Long Regime RMS','Interpreter','latex','FontSize',20);
                    imgData = obj.getSlice(squeeze(obj.fitCoeffsLong(:,:,:,3)),cfg,currFoc);
                    set(volAx,'Colormap',rgb2cmap([0,1,0;1,0,0],256));
                    clim(volAx,[0,0.1]);
                otherwise
                    fprintf("oops\n");
            end
            hold(volAx,'on');
            set(obj.uiObjs.volImg,'CData',imgData);
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


            roiData = [];
            roi = [];
            selectedRegion = obj.uiObjs.regionLB.Value;
            if ~strcmp(selectedRegion,'Crosshair')
                roi = obj.rois{strcmp(obj.roiNames,selectedRegion)};
                roiInds = find(roi);
                A = reshape(obj.aseData, [], size(obj.aseData,4));
                p = sum(A(roiInds,:), 1); %#ok
                s = sum(roi,'all');
                roiData = squeeze(p/s);
            end

            fullFoc = [obj.uiObjs.dim1Spinner.Value,...
                       obj.uiObjs.dim2Spinner.Value,...
                       obj.uiObjs.sliceSpinner.Value,...
                       obj.uiObjs.dTESpinner.Value];
            if ~isempty(obj.uiObjs.roiHandle) && isvalid(obj.uiObjs.roiHandle)
                delete(obj.uiObjs.roiHandle);
            end
            if isempty(roiData)
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
                roi_slice = roi(:,:,fullFoc(3));
                roi_slice = roi_slice';
                roi_TC = mat2TC(roi_slice,repmat([1,0,0],[256,1]));
                obj.uiObjs.roiHandle = imagesc(obj.uiObjs.upperAx,roi_TC);
                set(obj.uiObjs.roiHandle,'AlphaData',roi_slice.*0.3);
            end
            hold(obj.uiObjs.upperAx,'off');

            dTE_ax = obj.uiObjs.lowerAx;
            cla(dTE_ax);

            taus = obj.params.taus;
            tauParts = obj.params.tauPartitions;
            shortInds = tauParts.minTauIdx:tauParts.critInd_lower;
            longInds = tauParts.critInd_upper+1:tauParts.maxTauIdx;
            unusedInds = [1:tauParts.minTauIdx-1,...
                          tauParts.critInd_lower+1:tauParts.critInd_upper,...
                          tauParts.maxTauIdx+1:numel(taus)];
            hold(obj.uiObjs.lowerAx,'on');
            shortColor = [0.0667,0.4431,0.7451];
            longColor = [0.5216,0.0863,0.8196];
            exColor = [0.8667,0.3294,0];
            
            if ~isempty(roiData)
                focData = log(roiData);
                %fullFoc(1:3) = 1;
                roiIdx = find(strcmp(obj.roiNames,selectedRegion));
                roiRow = obj.roiStats(roiIdx,:);
                M_S = roiRow.M_S;
                B_S = roiRow.B_S;
                RMS_S = roiRow.RMS_S;
                M_L = roiRow.M_L;
                B_L = roiRow.B_L;
                RMS_L = roiRow.RMS_L;
                R2P = roiRow.R2P;
                vCBV = roiRow.vCBV;
                OEF = roiRow.OEF;
                Signal = roiData(currFoc(4));
            else
                spatIdxCell = num2cell(currFoc(1:3));
                focData = log(squeeze(obj.aseData(spatIdxCell{:},:)));
                cMatLong = obj.fitCoeffsLong;
                B_L = cMatLong(spatIdxCell{:},1);
                M_L = cMatLong(spatIdxCell{:},2);
                RMS_L = cMatLong(spatIdxCell{:},3);
                cMatShort = obj.fitCoeffsShort;
                B_S = cMatShort(spatIdxCell{:},1);
                M_S = cMatShort(spatIdxCell{:},2);
                RMS_S = cMatShort(spatIdxCell{:},3);
                
                R2P = obj.maps.R2P(spatIdxCell{:});
                vCBV = obj.maps.vCBV(spatIdxCell{:});
                OEF = obj.maps.OEF(spatIdxCell{:});
                Signal = obj.aseData(spatIdxCell{:},currFoc(4));
            end
            if all(focData<0)
                xlim(dTE_ax,[0,1]);
                ylim(dTE_ax,[0,1]);
                Signal = NaN;
                R2P = NaN;
                vCBV = NaN;
                OEF = NaN;
                text(dTE_ax,0.5,0.5,'No Valid Points!','FontSize',18,'HorizontalAlignment','center');
            else
                lgdHandles = [];
                lgdHandles(end+1) = scatter(dTE_ax,taus(shortInds),focData(shortInds),'filled','DisplayName','Short Regime','MarkerFaceColor',shortColor);
                lgdHandles(end+1) = scatter(dTE_ax,taus(longInds),focData(longInds),'filled','DisplayName','Long Regime','MarkerFaceColor',longColor);
    
                if obj.viewParams.showEx
                    lgdHandles(end+1) = scatter(dTE_ax,taus(unusedInds),focData(unusedInds),'Marker','x','DisplayName','Excluded','MarkerEdgeColor',exColor,'LineWidth',2);
                end
                
                M_L = M_L/1e3; %s
                M_S = M_S/1e6;  %s
                lgdHandles(end+1) = plot(dTE_ax,taus,M_L.*taus+B_L,'--','DisplayName','Long Regime Fit','Color',longColor);
                lgdHandles(end+1) = plot(dTE_ax,taus(shortInds),M_S.*taus(shortInds).^2+B_S,'--','DisplayName','Short Regime Fit','Color',shortColor);
    
                usedData = focData([shortInds,longInds]);
    
                yBds = [min(usedData),B_L];
                yRange = abs(yBds(2)-yBds(1));
                padding = 0.1;
                critTime = obj.params.critTime;
                critBuffer = obj.params.tCrit_buffer;
                xlim(dTE_ax,[taus(1),taus(end)]);
                ylim(dTE_ax,[yBds(1)-yRange*padding,max(yBds(2)+yRange*padding,max(usedData))]);
                lgdHandles(end+1) = xline(dTE_ax,critTime,'-','1.5x\ $t_{crit}$', ...
                                                 'Interpreter','latex', ...
                                                 'Color','black', ...
                                                 'FontSize',12,...
                                                 'DisplayName',sprintf('1.5x Critical Time (%0.1fms)',critTime));
                critRange = [critTime-critBuffer,critTime+critBuffer];
                lgdHandles(end+1) = xregion(dTE_ax,critRange, ...
                                                    'LineStyle','--', ...
                                                    'FaceColor',[1,0,0], ...
                                                    'EdgeColor',[0,0,0], ...
                                                    'EdgeAlpha',1, ...
                                                    'FaceAlpha',0.2,...
                                                    'DisplayName',sprintf('Buffer Window (%0.1fms)',2*critBuffer));
                lgdHandles(end+1) = xline(dTE_ax,taus(currFoc(4)),'-','$\tau_{current}$', ...
                                                   'Interpreter','latex', ...
                                                   'Color',[0.2314,0.6667,0.1961], ...
                                                   'FontSize',12,...
                                                   'DisplayName',sprintf('Current dTE (%0.1fms)',taus(currFoc(4))));
                legend(dTE_ax,lgdHandles,'Interpreter','latex','Location','northeast');
                
            end
            hold(dTE_ax,'off');
            xlabel(dTE_ax,'$2\tau$ [ms]','Interpreter','latex');
            ylabel(dTE_ax,'log(S) [a.u.]','Interpreter','latex');



            dTE_ax.Title.FontName = 'FixedWidth';
            title(dTE_ax,sprintf([' Signal: %-6.1f [a.u.]\n' ...
                                 '  R2'': %-6.1f [Hz]\n' ...
                                 'vCBV: %-6.1f [%%]\n' ...
                                 ' OEF: %-6.1f [%%]'], ...
                                 Signal,R2P,vCBV*100,OEF*100));
            obj.uiObjs.dim1Spinner.Value = fullFoc(1);
            obj.uiObjs.dim2Spinner.Value = fullFoc(2);
            obj.toggleSpinners("on");
        end
    end

    methods (Static)
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