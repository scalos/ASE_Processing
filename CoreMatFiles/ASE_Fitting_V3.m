function ASE_Fitting_V3(B0,hdr_File,sorted_ASE_File,opts)
    arguments
        B0 %Field strength
        hdr_File %Path to ase header file
        sorted_ASE_File %path to ase .nii(.gz) file. Can be cell array of files to be averaged
        opts.ROI_Files = {} %path to roi .nii(.gz) file. Can be cell array of multiple ROI's. Expected to already be aligned
        opts.segFiles = {}  %path to segentation file (or files). should match segLUT
        opts.segLUTs = {}    %path to segmentation LUT. Should match all segFiles
        opts.outPref = './ase_fitting_results'
        opts.Hct = 0.42
        opts.Hct_correction = 0.85;
        opts.dChi_0 = 0.18*1e-6;   % sus. diff. between oxy and deoxy blood
        opts.StO2_max = 0.60        % max expected StO2
        opts.verbose = false;
        opts.tauRange = []; %ms
        opts.tCrit_buffer = 3; %ms
        opts.writeOpts {mustBeMember(opts.writeOpts,{'maps','coeffs','stats'})} = {'maps','coeffs','stats'}; %outputs
        opts.overwriteOpts {mustBeMember(opts.overwriteOpts,{'maps','coeffs','stats'})} = {'maps','coeffs','stats'}; %overwriting
    end
    v = opts.verbose;
    if v; fprintf('Verifying Files...\n');end

    %ensure critical parameters
    assert(~isempty(hdr_File),'No header file was passed!');
    assert(~isempty(sorted_ASE_File),'No ASE file was passed!');

    assert(~isempty(opts.writeOpts),'No write options selected');
    if ~isempty(opts.ROI_Files)&&~isempty(opts.segFiles)
        error('Cannot set both ROI_Files and segFiles');
    end

    if xor(isempty(opts.segFiles),isempty(opts.segLUTs))
        error('A LUT and segmented volume pair is required if using segmented ROI option.');
    end

    %warn but allow empty ROI's
    if (isempty(opts.ROI_Files) && isempty(opts.segFiles)) && any(contains(opts.writeOpts,'stats'))
        %stats were asked for but either no ROI files were passed, or no
        %seg files were passed. either way there is no ROI information.
        warning('ROI Statistics output was requested but no ROI''s were passed. No statistics will be written.');
    end

    %ensure all lists are cells:
    if ~iscell(sorted_ASE_File); sorted_ASE_File = {sorted_ASE_File}; end
    if ~iscell(opts.ROI_Files); opts.ROI_Files = {opts.ROI_Files}; end
    if ~iscell(opts.segFiles); opts.segFiles = {opts.segFiles}; end
    if ~iscell(opts.segLUTs); opts.segLUTs = {opts.segLUTs}; end

    %if using segmented ROIs:
    segROI = ~isempty(opts.segLUTs);

    %parse header file:
    dTEs = readASEHeader(hdr_File);
    dTEs_sorted = sort(dTEs,1,"ascend");
    nReps = sum(dTEs_sorted==0);
    %write taus as dTE in ms:
    taus = dTEs_sorted(1:nReps:end)/1000;
    
    
    sorted_ASE_Data = combineAseReps(sorted_ASE_File);
    if v; fprintf('Averaged %d matching ASE files to use for fitting\n',numel(sorted_ASE_File)); end

    ROIs = {};
    ROI_Names = {};
    if any(contains(opts.writeOpts,'stats'))
        if segROI
            [ROIs,ROI_Names] = getROIsFromLUT(opts.segFiles,opts.segLUTs,"verbose",v);
        else
            [ROIs,ROI_Names] = getROIsFromCell(ROI_Files,"verbose",v);
        end
    end

    %%%%%%%%%%%% FITTING %%%%%%%%%%%%%%%%%%%%%%%%%

    
    gamma = 2.67502*1e8;            % rad/s/T
    dChi_0 = opts.dChi_0;             % sus. diff. between oxy and deoxy blood
    Hct = opts.Hct*opts.Hct_correction;               % corrected hematocrit
    StO2_max = opts.StO2_max;               % max expected StO2

    ase_sorted = double(sorted_ASE_Data);     
    tauRange = opts.tauRange;
    if isempty(tauRange)
        tauRange = [taus(1),taus(end)];
    end
    [~,minTauIdx] = min(abs(taus-tauRange(1)));
    [~,maxTauIdx] = min(abs(taus-tauRange(2)));
    nt = size(ase_sorted,4);
    assert(nt==numel(taus),"ASE has %d time points but tau series has %d",nt,numel(taus));
    if v; fprintf('Using times [%0.2fms:%0.2fms] out of %0.2fms (pts [%d:%d] of %d)\n',taus(minTauIdx),taus(maxTauIdx),taus(end),minTauIdx,maxTauIdx,numel(taus));end
    

    t_crit = 1/(gamma*(4/3)*pi*dChi_0*Hct*(1-StO2_max)*B0)*1000;
    t_crit = 1.5*t_crit;
    t_crit_lower = max(0,t_crit-opts.tCrit_buffer);
    t_crit_upper = t_crit+opts.tCrit_buffer;
    [~,critInd] = min(abs(taus-t_crit));
    [~,critInd_lower] = min(abs(taus-t_crit_lower));
    [~,critInd_upper] = min(abs(taus-t_crit_upper));
    shortTimes = [taus(minTauIdx),taus(critInd_lower)];
    longTimes = [taus(critInd_upper+1),taus(maxTauIdx)];

    if v; fprintf(['Computed critical time (1.5x): %0.2fms (at idx %d)\n' ...
                   '\t- %0.1fms buffer => [%0.1fms-%0.1fms] = crit-window\n' ...
                   '\t=> - Short Regime: [%d:%d]ms (frames [%d:%d])\n' ...
                   '\t   - Long Regime:  [%d:%d]ms (frames [%d:%d])\n'], ...
                   t_crit,critInd, ...
                   opts.tCrit_buffer,taus(critInd_lower),taus(critInd_upper),...
                   shortTimes(1),shortTimes(2),minTauIdx,critInd_lower, ...
                   longTimes(1),longTimes(2),critInd_upper,maxTauIdx);
    end

    shortInds = minTauIdx:critInd_lower;
    longInds = critInd_upper+1:maxTauIdx;
    ase_short = ase_sorted(:,:,:,shortInds);
    
    ase_long = ase_sorted(:,:,:,longInds);
    taus_long = taus(longInds);
    taus_short = taus(shortInds);

    if any(contains(opts.writeOpts,{'maps','coeffs'}))
        %If either maps or fit coeffs are asked for, do global fits
        if v;fprintf('Computing long-regime fits...\n');end
        cMat_long = fitLongRegime(ase_long,taus_long/1000);
    
        if v;fprintf('Computing short-regime fits...\n');end
        cMat_short = fitShortRegime(ase_short,taus_short/1000);
    
    
        tmpHdr = niftiinfo(sorted_ASE_File{1});
        if any(contains(opts.writeOpts,'coeffs'))
            %write fitting coefficients and RMS values:
            overwriteCoeffs = any(contains(opts.overwriteOpts,'coeffs'));
            cLong_path = sprintf("%s_CoeffsLong",opts.outPref);
            cLong_path_resolved = writeVec(tmpHdr,cMat_long,cLong_path,'verbose',v,...
                                                 'overwrite',overwriteCoeffs,...
                                                 'vecName','Long_Coeffs',...
                                                 'vecDesc','Lin. Fit Coeffs.',...
                                                 'intent','LinFitParams');
            cShort_path = sprintf("%s_CoeffsShort",opts.outPref);
            cShort_path_resolved = writeVec(tmpHdr,cMat_short,cShort_path,'verbose',v,...
                                                  'overwrite',overwriteCoeffs,...
                                                  'vecName','Short_Coeffs',...
                                                  'vecDesc','Lin. Fit Coeffs.',...
                                                  'intent','LinFitParams');
        end
    
        
        R2P = -cMat_long(:,:,:,2); %Hz
        b = cMat_long(:,:,:,1);
        vCBV = b-cMat_short(:,:,:,1);
        OEF = R2P./(vCBV*(4/3)*pi*gamma*dChi_0*Hct*B0);
    
        if any(contains(opts.writeOpts,'maps'))
            %write maps:
            overwriteMaps = any(contains(opts.overwriteOpts,'maps'));
            R2P_path = sprintf("%s_R2P",opts.outPref);
            R2P_path_resolved = writeMap(tmpHdr,R2P,R2P_path,'verbose',v,...
                                         'overwrite',overwriteMaps,...
                                         'mapName','R2P');
            
            vCBV_path = sprintf("%s_vCBV",opts.outPref);
            vCBV_path_resolved = writeMap(tmpHdr,vCBV,vCBV_path,'verbose',v,...
                                           'overwrite',overwriteMaps,...
                                           'mapName','vCBV');
        
            
            OEF_path = sprintf("%s_OEF",opts.outPref);
            OEF_path_resolved = writeMap(tmpHdr,OEF,OEF_path,'verbose',v,...
                                         'overwrite',overwriteMaps,...
                                         'mapName','OEF');
        end
    end
    
    if any(contains(opts.writeOpts,'stats'))
        %compute ROI stats
        roi_OEF = zeros(numel(ROIs),1);
        roi_R2P = zeros(numel(ROIs),1);
        roi_vCBV = zeros(numel(ROIs),1);
        roi_cMatLong = zeros(numel(ROIs),3);
        roi_cMatShort = zeros(numel(ROIs),3);
        dt = 0;
        totalROITime = 0;
        maxNameLen = max(cellfun(@(x) numel(x),ROI_Names));
        if ~isempty(ROIs)
            for idx = 1:numel(ROIs)
                roiTimer = tic;
                ROI = ROIs{idx};
                if isempty(ROI); continue; end
                if v
                    percentProg = round(idx/numel(ROIs)*100);
                    totalROITime = totalROITime + dt;
                    avgTime = totalROITime / idx;
                    estSec = avgTime * (numel(ROIs)-idx);
                    timeLeft = string(seconds(estSec),'hh:mm:ss');
                    fprintf('Computing stats for ROI: %*s | %3d%% (%3d/%3d) %s remaining\n', ...
                            maxNameLen,ROI_Names{idx},percentProg,idx,numel(ROIs),timeLeft);
                end
                
                ase_short_ROI = sum(ase_short.*ROI,[1,2,3])/sum(ROI,'all');
                ase_long_ROI = sum(ase_long.*ROI,[1,2,3])/sum(ROI,'all');
                cMat_long_ROI = fitLongRegime(ase_long_ROI,taus_long/1000);
                cMat_short_ROI = fitShortRegime(ase_short_ROI,taus_short/1000);
                R2P_ROI = -cMat_long_ROI(2); %Hz
                b_ROI = cMat_long_ROI(1);
                vCBV_ROI = b_ROI-cMat_short_ROI(1);
                OEF_ROI = R2P_ROI/(vCBV_ROI*(4/3)*pi*gamma*dChi_0*Hct*B0);
                roi_OEF(idx) = OEF_ROI;
                roi_R2P(idx) = R2P_ROI;
                roi_vCBV(idx) = vCBV_ROI;
                roi_cMatLong(idx,:) = cMat_long_ROI;
                roi_cMatShort(idx,:) = cMat_short_ROI;
                dt = toc(roiTimer);
            end
        end

        %Purge bad rows:
        goodRows_Bin = zeros(numel(ROI_Names),1);
        for idx = 1:numel(ROI_Names)
            if ~isempty(ROI_Names{idx})
                goodRows_Bin(idx) = 1;
            end
        end
        goodRows = find(goodRows_Bin==1);
        badRows = setdiff(1:numel(ROI_Names),goodRows);
        ROI_Names(badRows) = [];
        roi_OEF = roi_OEF(goodRows);
        roi_R2P = roi_R2P(goodRows);
        roi_vCBV = roi_vCBV(goodRows);
        roi_cMatShort = roi_cMatShort(goodRows,:);
        roi_cMatLong = roi_cMatLong(goodRows,:);

        roiData = table(roi_OEF,roi_R2P,roi_vCBV,roi_cMatShort(:,2),roi_cMatShort(:,1),roi_cMatShort(:,3),...
                                                 roi_cMatLong(:,2),roi_cMatLong(:,1),roi_cMatLong(:,3),...
                        'VariableNames',{'OEF','R2P','vCBV','M_S','B_S','RMS_S','M_L','B_L','RMS_L'},...
                        'RowNames',ROI_Names);
        overwriteStats = any(contains(opts.overwriteOpts,'stats'));
        ROI_Tbl_Name = findGoodPath(sprintf("%s_ROI_Stats.csv",opts.outPref),overwriteStats);
        writetable(roiData, ROI_Tbl_Name, 'WriteRowNames', true, 'Delimiter', ',');
        if v;fprintf('Wrote ROI stats to %s\n',ROI_Tbl_Name);end
    end

    %%write manifest%%
    manifest = struct();
    runTimestamp = string(datetime('now',"Format","yyyyMMdd_HHmmSS"));
    manifest.runTimestamp = runTimestamp;
    manifest.B0 = B0;
    manifest.hdr_File = hdr_File;
    manifest.sorted_ASE_File = sorted_ASE_File;
    manifest.tauRange = opts.tauRange;
    manifest.tCrit_buffer = opts.tCrit_buffer;
    manifest.taus = taus;
    manifest.tauPartitions = struct('minTauIdx',minTauIdx,'critInd_lower',critInd_lower, ...
                                     'critInd',critInd,'critInd_upper',critInd_upper,'maxTauIdx',maxTauIdx);
    manifest.critTime = t_crit;
    
    manifest.paths = struct();
    if any(contains(opts.writeOpts,'maps'))
        manifest.paths.R2P  = R2P_path_resolved;   % actual path returned by findGoodPath inside writeMap
        manifest.paths.vCBV = vCBV_path_resolved;
        manifest.paths.OEF  = OEF_path_resolved;
    end
    if any(contains(opts.writeOpts,'coeffs'))
        manifest.paths.coeffsLong  = cLong_path_resolved;
        manifest.paths.coeffsShort = cShort_path_resolved;
    end
    if any(contains(opts.writeOpts,'stats'))
        manifest.paths.statsCSV = ROI_Tbl_Name;    % findGoodPath's return value, already captured
    end
    
    manifest.segFiles = opts.segFiles;
    manifest.segLUTs  = opts.segLUTs;
    manifest.ROI_Files = opts.ROI_Files;
    
    manifestPath = sprintf("%s_manifest_%s.json", opts.outPref, runTimestamp);
    fid = fopen(manifestPath, 'w');
    fprintf(fid, '%s', jsonencode(manifest,'PrettyPrint',true));
    fclose(fid);
    
    if v; fprintf('Wrote manifest to %s\n', manifestPath); end
    

    if v;fprintf('DONE!\n');end
end





