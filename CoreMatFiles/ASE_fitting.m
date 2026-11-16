function ASE_fitting(B0,hdr_file,ase_files,opts)
    %expected data formats:
        % ase_files:
        %   - cell array of ase echos. same sort as applied to hdr_files
        %
        % hdr_Files: single header file containing info for all echos.
        %   - should be readable by readASEHeader.m and contain:
        %       - dTEs (per echo)
        %       - zshimming (per echo) 
        %       - TE (global)
    arguments
        % Required inputs:
        B0 %Field strength
        hdr_file %Path to ase header file 
        ase_files (:,:) %path to ase .nii(.gz) file.

        % Fitting options:
        opts.tauRange = []; %ms
        opts.tCrit_buffer = [3 12]; %ms; [below above] t_crit (scalar = symmetric). Above > below because
                                    %   the long regime is reached late when true t_c > assumed (StO2_max); see docs Q3
        opts.averageReps = false %whether to average repetitions before fitting 
        opts.dBCorr = false %Through-slice (z-shim based) dB correction; settings in opts.quality.zshim_*
        opts.tauUsePts = []; %use to restrict tau points used in all fitting operations.  ndTE x nEchos
        opts.physBd_R2 = [0,inf];
        opts.physBd_R2P = [0,inf];
        opts.physBd_vCBV = [0,1];
        opts.physBd_OEF = [0,1];
        opts.physBdClipVal = [] %if empty: clip to bound, if not: set to value;
        opts.normErr = true;

        % Voxel-quality options (see qualityDefaults.m; every feature off by default):
        opts.quality = struct()

        % ROI computation options:
        opts.segFiles = {}  %path to segentation file (or files). should match segLUT
        opts.segLUTs = {}    %path to segmentation LUT. Should match all segFiles
        opts.minOverlapFrac = 0.25 %fraction ase voxel must overlap with region to be added to that ROI
        opts.warnSegIssues = true;
        opts.qualityMaskFiles = {} %tissue-fraction masks on the ASE grid (e.g. brain minus
                                   %ventricles); CSF fraction = 1 - fraction (see quality.useCSF)
        
        % Fitting constants:
        opts.Hct = 0.42
        opts.Hct_correction = 0.85;
        opts.dChi_0 = 0.18*1e-6;   % sus. diff. between oxy and deoxy blood
        opts.StO2_max = 0.60        % max expected StO2 (sets t_crit only). REVISIT: 0.6 is a typical, not a
                                    %   maximal, venous Y; see docs/quality-mask-report.md (Q3, StO2_max)

        % Write options:
        opts.outPref = './ase_fitting_results'
        opts.writeOpts {mustBeMember(opts.writeOpts,{'maps','coeffs','stats'})} = {'maps','coeffs','stats'}; %outputs
        opts.overwriteOpts {mustBeMember(opts.overwriteOpts,{'maps','coeffs','stats'})} = {'maps','coeffs','stats'}; %overwriting
        opts.writeRawMaps = false;

        % Other options:
        opts.verbose = false;
    end
    v = opts.verbose;
    q = qualityDefaults(opts.quality);
    opts.quality = q;          % resolved block is recorded in the manifest
    % which quality diagnostics are needed
    needExt   = q.fitExtensions || (q.enabled && (q.useQuadSign || q.useVcbvZ || q.useR2PEcho));
    needNoise = q.repNoise || needExt || q.noiseWeighting || (q.enabled && q.useSNR);

    %% Load files and verify inputs
    if v; fprintf('Verifying Files...\n');end

    %ensure critical parameters
    assert(~isempty(hdr_file),'No header file was passed!');
    assert(~isempty(ase_files),'No ASE file(s) were passed!');
    assert(~isempty(opts.writeOpts),'No write options selected');
    
    if xor(isempty(opts.segFiles),isempty(opts.segLUTs))
        error("Both segFiles and segLUTs need to be provided for ROI analysis");
    end

    %warn but allow empty ROI's
    if (isempty(opts.segFiles)) && any(contains(opts.writeOpts,'stats'))
        %stats were asked for but either no ROI files were passed, or no
        %seg files were passed. either way there is no ROI information.
        warning('ROI Statistics output was requested but no ROI''s were passed. No statistics will be written.');
    end

    %ensure all lists are cells:
    if ~iscell(ase_files); ase_files = {ase_files}; end
    if ~iscell(opts.segFiles); opts.segFiles = {opts.segFiles}; end
    if ~iscell(opts.segLUTs); opts.segLUTs = {opts.segLUTs}; end

    %parse header file:
    [taus_raw,zShim,TEs,info] = readASEHeader(hdr_file);
    nReps  = info.nReps;
    nEchos = info.nEchos;
    nSweep = info.nSweep;
    ndTEs  = numel(taus_raw);

    %% Frame-use mask (ndTEs x nEchos): user mask AND z-shim == 0
    tauUsePts = opts.tauUsePts;
    if isempty(tauUsePts)
        tauUsePts = ones(ndTEs,nEchos);
    end
    tauUsePts = double(tauUsePts);               % must be numeric so 0 -> NaN works later
    assert(isequal(size(tauUsePts),[ndTEs,nEchos]), ...
        'tauUsePts must be ndTEs x nEchos.');

    % frames with a nonzero z-moment never enter the regime fits (with or without dB
    % correction); they are only used to estimate the through-slice gradient Gz.
    tauUsePts(zShim~=0) = 0;
    needZshim = opts.dBCorr || q.zshimEstimate || ...
                (q.enabled && (q.useZshimEdge || q.useZshimEcho || q.useZshimFrac));
    if needZshim && any(isnan(info.SliceThickness_mm))
        if opts.dBCorr
            error(['dBCorr requires SliceThickness_mm in the header (%s). ' ...
                   'Regenerate it with asePreProcess.sh --regenHdr.'], hdr_file);
        end
        warning('SliceThickness_mm missing from %s: z-shim gradient estimate skipped.', hdr_file);
        needZshim = false;
    end
    if needZshim && info.nZshimming == 0
        if opts.dBCorr, error('dBCorr requested but the header lists no z-shim frames.'); end
        warning('No z-shim frames in %s: z-shim gradient estimate skipped.', hdr_file);
        needZshim = false;
    end

    %% Handle averaging and sort ase data into single volume
    % For every frame of the combined data, record which echo and which raw-file frame
    % it came from (frameEcho, frameSrc); taus/TEs/mask are all derived from this map.
    repInds   = info.repInds;                          % nReps x nSweep raw-frame indices
    zShimInds = nSweep*nReps + (1:info.nZshimming);    % trailing z-shim frames (raw indices)
    nZshim    = numel(zShimInds);
    echoIdx   = (1:nEchos);
    if opts.averageReps
        if v; fprintf('Averaging %d reps of each of %d echos for fitting...\n',nReps,nEchos); end
        aseData_combined = combineAseReps(ase_files,info,"method",'avg');
        % layout: [e1 sweep | e2 sweep | ... | e1 zshim | e2 zshim | ...]
        frameEcho = [reshape(repmat(echoIdx,nSweep,1),[],1); reshape(repmat(echoIdx,nZshim,1),[],1)];
        frameSrc  = [repmat(repInds(1,:).',nEchos,1);      repmat(zShimInds(:),nEchos,1)];
        % an averaged sweep point is usable only if it is usable in every rep
        sweepUse = true(nSweep,nEchos);
        for r = 1:nReps
            sweepUse = sweepUse & (tauUsePts(repInds(r,:),:) ~= 0);
        end
        usePts_full = [sweepUse(:); reshape(tauUsePts(zShimInds,:) ~= 0,[],1)];
    else
        if v; fprintf('Stacking echos unaveraged for fitting (reps fit as extra samples)...\n');end
        aseData_combined = combineAseReps(ase_files,info,"method","stack");
        % layout: [all e1 frames | all e2 frames | ...]
        frameEcho = reshape(repmat(echoIdx,ndTEs,1),[],1);
        frameSrc  = repmat((1:ndTEs).',nEchos,1);
        usePts_full = tauUsePts(:) ~= 0;
    end
    nFrames = numel(frameSrc);
    templateHdr = niftiinfo(ase_files{1});

    %% CSF fraction from the quality masks (only when quality features are in use)
    qual = struct();
    anyQuality = q.enabled || q.repNoise || q.zshimEstimate || q.fitExtensions;
    if ~iscell(opts.qualityMaskFiles); opts.qualityMaskFiles = {opts.qualityMaskFiles}; end
    opts.qualityMaskFiles = opts.qualityMaskFiles(~cellfun(@isempty, opts.qualityMaskFiles));
    if anyQuality && ~isempty(opts.qualityMaskFiles)
        csf = zeros(templateHdr.ImageSize(1:3));
        for i = 1:numel(opts.qualityMaskFiles)
            fq = opts.qualityMaskFiles{i};
            isequalNiftiGeometry({templateHdr, niftiinfo(fq)}, "onMismatch","error");
            frac = min(max(double(niftiread(fq)),0),1);    % linear resampling -> tissue fraction
            csf = max(csf, 1 - frac);                       % several masks: worst case
        end
        qual.csfFrac = csf;
        if v; fprintf('Loaded %d quality mask(s); CSF fraction > %.2f in %.1f%% of voxels\n', ...
                numel(opts.qualityMaskFiles), q.csf_max, 100*mean(csf(:) > q.csf_max)); end
    end

    %% Repetition noise (quality diagnostics), always from the unaveraged data
    if needNoise
        if v; fprintf('Estimating repetition noise (sigma, SNR) per echo...\n'); end
        baseInds = find(taus_raw(repInds(1,:)) == 0);
        qual.sigma = NaN([templateHdr.ImageSize(1:3), nEchos]);
        qual.S0    = NaN([templateHdr.ImageSize(1:3), nEchos]);
        for e = 1:nEchos
            if opts.averageReps
                echoData = niftiread(ase_files{e});              % unaveraged frames of this echo
            else
                echoData = aseData_combined(:,:,:,(e-1)*ndTEs + (1:ndTEs));
            end
            [qual.sigma(:,:,:,e), qual.S0(:,:,:,e)] = estimateRepNoise(echoData, repInds, 'baseInds', baseInds);
        end
        clear echoData
        qual.SNR = qual.S0 ./ qual.sigma;
    end
    nt = size(aseData_combined,4);
    assert(nt==nFrames,"Combined ASE data has %d frames but the frame map expects %d.",nt,nFrames);
    %{
    At this point:
    taus_raw: ndTEs x 1, dTE per raw-file frame (acquisition order; includes every rep and z-shim)
    TEs: nEchos x 1
    zShim: ndTEs x nEchos, in Tms/m
    frameEcho, frameSrc: nFrames x 1, echo and raw frame index of each combined frame
    usePts_full: nFrames x 1 logical; false = frame excluded (z-shim or user mask)
    aseData_combined: X x Y x Z x nFrames
        averaged: nFrames = (nSweep + nZshim)*nEchos
        stacked:  nFrames = ndTEs*nEchos
    %}

    %% Combine Echos
    taus_full = taus_raw(frameSrc) + (TEs(frameEcho) - TEs(1));
    TEs_full  = TEs(frameEcho);
    taus_full = taus_full(:);  TEs_full = TEs_full(:);

    %% FITTING 
    gamma = 2.67502*1e8;            % rad/s/T
    dChi_0 = opts.dChi_0;             % sus. diff. between oxy and deoxy blood
    Hct = opts.Hct*opts.Hct_correction;               % corrected hematocrit
    StO2_max = opts.StO2_max;               % max expected StO2

    %sort out short/long regimes:
    ase_sorted = double(aseData_combined); 
    clear aseData_combined %free memory
    tauRange = opts.tauRange;
    minTau = min(taus_full,[],'omitmissing');
    maxTau = max(taus_full,[],'omitmissing');
    if isempty(tauRange)
        tauRange = [minTau,...
                    maxTau];
    end
    [~,minTauIdx] = min(abs(taus_full-tauRange(1)),[],"omitmissing");
    [~,maxTauIdx] = min(abs(taus_full-tauRange(2)),[],"omitmissing");
    if v; fprintf('Using times [%0.2fms:%0.2fms] out of %0.2fms\n', ...
            taus_full(minTauIdx),taus_full(maxTauIdx),maxTau);end
    

    t_crit = 1/(gamma*(4/3)*pi*dChi_0*Hct*(1-StO2_max)*B0)*1000;
    t_crit = round(1.5*t_crit,3);
    tcBuf = opts.tCrit_buffer(:).';
    assert(any(numel(tcBuf) == [1 2]) && all(tcBuf >= 0), 'tCrit_buffer must be a scalar or [below above], >= 0 (ms).');
    if isscalar(tcBuf), tcBuf = [tcBuf tcBuf]; end
    t_crit_lower = max(0,t_crit-tcBuf(1));
    t_crit_upper = t_crit+tcBuf(2);
    [~,critInd] = min(abs(taus_full-t_crit));
    [~,critInd_lower] = min(abs(taus_full-t_crit_lower));
    [~,critInd_upper] = min(abs(taus_full-t_crit_upper));
    absTau    = abs(taus_full);                                   % NaN for z-shim frames
    inRange   = taus_full >= tauRange(1) & taus_full <= tauRange(2);
    shortMask = inRange & absTau <= t_crit_lower;
    longMask  = inRange & absTau >  t_crit_upper;
    shortInds  = find(shortMask);   
    longInds = find(longMask);
    shortTimes = taus_full(shortInds);
    longTimes = taus_full(longInds);
    shortIndsPretty = compressIntArr(shortInds);
    longIndsPretty = compressIntArr(longInds);
    if v; fprintf(['Computed critical time (1.5x): %0.2fms (at idx %d)\n' ...
                   '\t- -%0.1f/+%0.1fms buffer => [%0.1fms-%0.1fms] = crit-window\n' ...
                   '\t=> - Short Regime: [%0.1f:%0.1f]ms (frames %s)\n' ...
                   '\t   - Long Regime:  [%0.1f:%0.1f]ms (frames %s)\n'], ...
                   t_crit,critInd, ...
                   tcBuf(1),tcBuf(2),taus_full(critInd_lower),taus_full(critInd_upper),...
                   min(shortTimes,[],'omitmissing'),max(shortTimes,[],'omitmissing'),shortIndsPretty, ...
                   min(longTimes,[],'omitmissing'),max(longTimes,[],'omitmissing'),longIndsPretty);
    end
    
    %% Z-shim: through-slice gradient estimate and (optional) correction
    if needZshim
        szV = size(ase_sorted,1:3);
        dz  = info.SliceThickness_mm/1000;                 % m, per echo
        qual.Gz         = NaN([szV, nEchos]);
        qual.zshimAmp   = NaN([szV, nEchos]);
        qual.zshimAtten = NaN([szV, nEchos]);
        qual.zshimEdge  = zeros([szV, nEchos]);
        tDeph = NaN(nEchos,1);
        for e = 1:nEchos
            idx = find(frameEcho == e & ismember(frameSrc, info.zShimInds));
            [~,ord] = sort(frameSrc(idx));  idx = idx(ord);   % same order as info.zMoments
            tD = abs(taus_full(idx));
            assert(all(abs(tD - tD(1)) < 1e-6), 'Echo %d z-shim frames have different dephasing times.', e);
            o = fitZshimGradient(ase_sorted(:,:,:,idx), info.zMoments(:,e), tD(1)/1000, dz(e), ...
                                 'GzRange', q.zshim_GzRange, 'GzStep', q.zshim_GzStep);
            qual.Gz(:,:,:,e) = o.Gz;          qual.zshimAmp(:,:,:,e)   = o.A;
            qual.zshimAtten(:,:,:,e) = o.atten; qual.zshimEdge(:,:,:,e) = double(o.edge);
            tDeph(e) = tD(1);
        end
        % summaries over voxels with real signal (echo-1 amplitude > 0.2*p99; background
        % voxels would otherwise dominate)
        A1 = qual.zshimAmp(:,:,:,1);
        sig = A1 > 0.2*prctile(A1(:),99);
        if v
            for e = 1:nEchos
                Ge = qual.Gz(:,:,:,e);  Ee = qual.zshimEdge(:,:,:,e);
                fprintf('Z-shim echo %d: t_deph=%.1fms, median Gz=%.1f uT/m, %.1f%% edge-flagged (voxels with signal)\n', ...
                    e, tDeph(e), median(Ge(sig),'omitnan'), 100*mean(Ee(sig)>0));
            end
        end
        % attenuation per voxel and frame (each echo's frames use that echo's Gz); voxels
        % without a Gz estimate (no signal) are left uncorrected (F = 1)
        usedLong = find(longMask & usePts_full);
        nLow = zeros(szV);
        % hard voxel x frame exclusion below Fcut; with F^2 weighting only frames where the
        % sinc model itself is unreliable are dropped, otherwise (legacy) below Fmin
        if q.zshim_weightF2, Fcut = q.zshim_Fcut; else, Fcut = q.zshim_Fmin; end
        for f = 1:nFrames
            if ~usePts_full(f), continue; end
            e = frameEcho(f);
            F = zshimAttenuation(qual.Gz(:,:,:,e), taus_full(f)/1000, dz(e));
            F(isnan(F)) = 1;
            low = F < q.zshim_Fmin;
            if ismember(f, usedLong), nLow = nLow + low; end
            if opts.dBCorr
                fr = ase_sorted(:,:,:,f) ./ F;
                fr(F < Fcut) = NaN;                        % voxel x frame exclusion
                ase_sorted(:,:,:,f) = fr;
            end
        end
        qual.zshimFracLow = nLow / max(numel(usedLong),1);
        qual.zshimFracLow(all(isnan(qual.Gz),4)) = NaN;
        if v && opts.dBCorr
            if q.zshim_weightF2, wTxt = sprintf('F^2 weights, NaN below F=%.2f', Fcut);
            else, wTxt = sprintf('unweighted, NaN below F=%.2f', Fcut); end
            fprintf(['Applied z-shim dB correction (%s); Fmin=%.2f: in voxels with signal, %.1f%% have >=1 ' ...
                     'long frame below Fmin, %.1f%% more than zshimFrac_max\n'], wTxt, q.zshim_Fmin, ...
                100*mean(qual.zshimFracLow(sig) > 0), 100*mean(qual.zshimFracLow(sig) > q.zshimFrac_max));
        end
    end

    %handle optional further restrictions on fitting points:
    %excluded frames -> NaN (keeps size; fitting ignores NaNs)
    frameMask = ones(1,1,1,nFrames);
    frameMask(~usePts_full) = NaN;
    ase_sorted = ase_sorted.*frameMask;   % broadcasts, no repmat

    %separate short and long regimes:

    longInds_good = longInds(~isnan(taus_full(longInds)));
    taus_long = taus_full(longInds_good);
    TEs_long = TEs_full(longInds_good);

    shortInds_good = shortInds(~isnan(taus_full(shortInds)));
    taus_short = taus_full(shortInds_good);
    TEs_short = TEs_full(shortInds_good);

    ase_short = ase_sorted(:,:,:,shortInds_good,:);
    ase_long = ase_sorted(:,:,:,longInds_good,:);

    % z-shim correction weights: dividing by F scales the log-space noise by 1/F, so each
    % corrected voxel-frame gets weight F^2 in the fits (empty = unweighted, the default)
    wLong = [];  wShort = [];
    if needZshim && opts.dBCorr && q.zshim_weightF2
        wLong  = zshimWeights(qual.Gz, taus_full(longInds_good),  frameEcho(longInds_good),  dz);
        wShort = zshimWeights(qual.Gz, taus_full(shortInds_good), frameEcho(shortInds_good), dz);
    end
        

    resolvedPaths = struct;
    physBds = struct("R2",opts.physBd_R2,...
                     "R2P",opts.physBd_R2P,...
                     "vCBV",opts.physBd_vCBV,...
                     "OEF",opts.physBd_OEF);
    wantStats = any(contains(opts.writeOpts,'stats')) && ~isempty(opts.segFiles);
    if any(contains(opts.writeOpts,{'maps','coeffs'})) || (q.enabled && wantStats)
        %If either maps or fit coeffs are asked for, do global fits.
        %The voxel quality mask also needs them, so they always run when it is applied to stats.
        fitRes = fitASERegimes(ase_long, ase_short, taus_long, TEs_long, taus_short, TEs_short, nEchos, ...
                    (4/3)*pi*gamma*dChi_0*Hct*B0, [], wLong, wShort);
        B_long = fitRes.BL;
        B_short = fitRes.BS;
        % rms_long = fit.rmsL;
        % rms_short = fit.rmsS;
        rawMaps = struct('R2',fitRes.R2,'R2P',fitRes.R2P,'vCBV',fitRes.vCBV,'OEF',fitRes.OEF);

        %% fit extensions (quality diagnostics): vCBV consistency and per-echo R2'
        if needExt
            if ~isempty(fitRes.k)
                qual.vCBV_quad = fitRes.vCBV_quad;
                % noise factor per combined frame: averaged sweep frames have sigma/sqrt(nReps)
                facFrame = ones(nFrames,1);
                if opts.averageReps, facFrame(1:nSweep*nEchos) = 1/sqrt(nReps); end
                sigLog = struct('sigma',qual.sigma, ...
                    'echoL',frameEcho(longInds_good), 'facL',facFrame(longInds_good), ...
                    'echoS',frameEcho(shortInds_good),'facS',facFrame(shortInds_good), ...
                    'wL',wLong,'wS',wShort);
                c = vcbvConsistency(ase_long, ase_short, fitRes.XL, fitRes.XS, fitRes.TES, fitRes, sigLog);
                qual.vCBV_SE = c.SE_int;  qual.vCBVquad_SE = c.SE_quad;  qual.vCBV_z = c.z;
            else
                warning('Short regime has no tau^2 column: vCBV_quad / vCBV_z not computed.');
            end
            if nEchos == 2
                pe = fitR2PerEcho(ase_long, taus_long, frameEcho(longInds_good), fitRes.rmsL, fitRes.nValidL, wLong);
                qual.R2P_e1 = pe.R2P_e1;  qual.R2P_e2 = pe.R2P_e2;
                qual.R2P_Fstat = pe.Fstat;  qual.R2P_pval = pe.pval;
            end
        end

        %% voxel quality mask (uses raw, unclamped maps)
        if q.enabled
            [qual.keep, qual.flags, qual.bitInfo, qual.qInfo] = computeVoxelQuality(rawMaps, qual, physBds, q);
            if v && q.useVein
                fprintf('Vein criterion (%s): threshold vCBV > %.4g (reference median %.4g over %d tissue voxels)\n', ...
                    q.vein_mode, qual.qInfo.veinThreshold, qual.qInfo.veinRef, qual.qInfo.veinRefN);
            end
            if v
                fprintf('Voxel quality: %d of %d voxels with a fit kept (whole FOV)\n', ...
                    nnz(qual.keep & isfinite(rawMaps.vCBV)), nnz(isfinite(rawMaps.vCBV)));
                disp(qual.bitInfo);
            end
        end

        C = exp(B_short(:,:,:,1)); %amplitude
        [imErr,~] = computeASEError(ase_sorted,taus_full,TEs_full,rawMaps.R2,rawMaps.R2P,rawMaps.vCBV,C,"normalize",opts.normErr);

        maps = applyPhysBounds(rawMaps, physBds, opts.physBdClipVal);
        
        if any(contains(opts.writeOpts,'coeffs'))
            %write fitting coefficients and RMS values:
            overwriteCoeffs = any(contains(opts.overwriteOpts,'coeffs'));
            cLong_path = sprintf("%s_CoeffsLong",opts.outPref);

            resolvedPaths.CoeffsLong = writeVec(templateHdr,B_long,cLong_path,'verbose',v,...
                                                 'overwrite',overwriteCoeffs,...
                                                 'vecName','Long_Coeffs',...
                                                 'vecDesc','Lin. Fit Coeffs.',...
                                                 'intent','LinFitParams');
            cShort_path = sprintf("%s_CoeffsShort",opts.outPref);
            resolvedPaths.CoeffsShort = writeVec(templateHdr,B_short,cShort_path,'verbose',v,...
                                                  'overwrite',overwriteCoeffs,...
                                                  'vecName','Short_Coeffs',...
                                                  'vecDesc','Lin. Fit Coeffs.',...
                                                  'intent','LinFitParams');
            % fit RMS (log space) per regime
            resolvedPaths.RMS_L = writeMap(templateHdr, fitRes.rmsL, sprintf("%s_RMS_L",opts.outPref), ...
                                           'verbose',v, 'overwrite',overwriteCoeffs, 'mapName','RMS_L');
            resolvedPaths.RMS_S = writeMap(templateHdr, fitRes.rmsS, sprintf("%s_RMS_S",opts.outPref), ...
                                           'verbose',v, 'overwrite',overwriteCoeffs, 'mapName','RMS_S');
        end
        %add error map to list of maps to write
        maps.imErr = imErr;
        if any(contains(opts.writeOpts,'maps'))
            %write maps:
            overwriteMaps = any(contains(opts.overwriteOpts,'maps'));
            outNames = fieldnames(maps);
            for i = 1:numel(outNames)
                nm = outNames{i};
                map = maps.(nm);
                if isempty(map), continue; end
                pth_Bd = writeMap(templateHdr, map, sprintf("%s_%s", opts.outPref, nm), ...
                         'verbose',v, 'overwrite',overwriteMaps, 'mapName',nm);
                resolvedPaths.(nm) = pth_Bd;
                if opts.writeRawMaps && isfield(rawMaps,nm) && ~isempty(rawMaps.(nm))
                    pth_Raw = writeMap(templateHdr, rawMaps.(nm), sprintf("%s_%s_raw", opts.outPref, nm), ...
                         'verbose',v, 'overwrite',overwriteMaps, 'mapName',nm);
                    resolvedPaths.(sprintf("%s_raw",nm)) = pth_Raw;
                end
            end
        end
    end
    
    %% write quality diagnostic maps
    if any(contains(opts.writeOpts,'maps')) && ~isempty(fieldnames(qual))
        overwriteMaps = any(contains(opts.overwriteOpts,'maps'));
        resolvedPaths = writeQualityMaps(qual, templateHdr, opts.outPref, nEchos, overwriteMaps, v, resolvedPaths);
    end

    if any(contains(opts.writeOpts,'stats')) && ~isempty(opts.segFiles)
        if v; fprintf('Organizing regions...\n'); end
        orig_state = warning;
        if ~opts.warnSegIssues, warning('off','all'); end
        [segVol,LUT] = combineSegVols(opts.segFiles,opts.segLUTs,"verbose",v,"keepOverlaps",1);
        warning(orig_state);
        segInfo = niftiinfo(opts.segFiles{1});
        M = segToAseMembership(segVol, LUT, segInfo, templateHdr,"maxSampling",[5,5,5],"minFrac",opts.minOverlapFrac);
        nVoxReg = full(sum(M,2));            % effective voxel count (weighted)
        % remove the size-equality assert and the buildMembership call
        nReg = height(LUT);

        % ---- voxel quality: hard exclusion (and optional noise weighting) ----
        useQ = q.enabled && isfield(qual,'keep');
        if useQ
            M0 = M;                                            % occupancy before exclusion
            M  = M .* double(qual.keep(:).');
            nVoxKept = full(sum(M,2));
            if q.noiseWeighting
                % weight ~ 1/var(log S) = S^2/sigma^2 (echo-1 SNR^2); scaled so a region's
                % total weight stays comparable. Off by default until compared.
                w = qual.SNR(:,:,:,1).^2;  w(~isfinite(w)) = 0;
                M = M .* (w(:).' / median(w(qual.keep & w > 0)));
            end
            lowKept = nVoxKept < q.minVoxKept;
            if v; fprintf('Quality mask: %.1f%% of ROI occupancy excluded; %d region(s) below minVoxKept=%g\n', ...
                    100*(1 - sum(nVoxKept)/sum(nVoxReg)), nnz(lowKept & nVoxReg > 0), q.minVoxKept); end
        end
        
    
        % ---- z-shim: one fixed weight per voxel for the region means ----
        % Same weight in every frame, so every frame of a region curve averages the same voxel
        % mix (no curve-shape bias); weight = F^2 at the voxel's worst used frame, and 0 if any
        % used frame was dropped (F < Fcut), which keeps the voxel set identical across frames.
        if ~isempty(wLong)
            vw = min(min(wLong,[],[4 5]), min(wShort,[],[4 5]));
            vw(vw < q.zshim_Fcut^2) = 0;
            occTot  = full(sum(M,'all'));
            occZero = full(sum(M(:, vw(:) == 0),'all'));
            M = M .* double(vw(:).');
            if v; fprintf(['Z-shim ROI voxel weights (F^2 at worst used frame): total ROI weight %.1f%% of ' ...
                    'occupancy; %.1f%% of occupancy at weight 0 (a frame below Fcut)\n'], ...
                    100*full(sum(M,'all'))/occTot, 100*occZero/occTot); 
            end
        end

        if v; fprintf('Averaging signal and fitting %d regions...\n',nReg); end
        meanLong  = regionMeans(ase_long,  M);
        meanShort = regionMeans(ase_short, M);
    
        oefScale = (4/3)*pi*gamma*dChi_0*Hct*B0;
        wLongR = [];  wShortR = [];
        if ~isempty(wLong)
            wLongR  = regionWeights(ase_long,  wLong,  M);
            wShortR = regionWeights(ase_short, wShort, M);
        end
        roi = fitASERegimes(meanLong, meanShort, taus_long, TEs_long, taus_short, TEs_short, nEchos, oefScale, [], ...
                            wLongR, wShortR);
    
        if isempty(roi.R2), roi.R2 = NaN(nReg,1); end       % single echo: R2 unavailable
        roiRaw = struct('R2',roi.R2,'R2P',roi.R2P,'vCBV',roi.vCBV,'OEF',roi.OEF);
        roiBd  = applyPhysBounds(roiRaw, physBds, opts.physBdClipVal);
        nonEmpty = nVoxReg > 0;                              % regions with no voxels -> NaN
        if useQ, nonEmpty = nonEmpty & ~lowKept; end         % too few kept voxels -> NaN
        for f = fieldnames(roiBd)'
            roiBd.(f{1})(~nonEmpty) = NaN;  roiRaw.(f{1})(~nonEmpty) = NaN;
        end
    
        roiData = table(roiBd.OEF(:), roiBd.R2P(:), roiBd.vCBV(:), roiBd.R2(:), ...
                        reshape(roi.BL,nReg,[]), roi.rmsL(:), roi.nValidL(:), ...
                        reshape(roi.BS,nReg,[]), roi.rmsS(:), roi.nValidS(:), nVoxReg, ...
                        'VariableNames',{'OEF','R2P','vCBV','R2','B_L','RMS_L','nFrames_L', ...
                                         'B_S','RMS_S','nFrames_S','nVox'}, ...
                        'RowNames',LUT.name);
        if opts.writeRawMaps
            for f = {'OEF','R2P','vCBV','R2'}
                roiData.([f{1} '_raw']) = roiRaw.(f{1})(:);
            end
        end
        if useQ
            roiData = addQualityColumns(roiData, M0, nVoxReg, nVoxKept, qual, rawMaps);
        end
    
        overwriteStats = any(contains(opts.overwriteOpts,'stats'));
        ROI_Tbl_Name = findGoodPath(sprintf("%s_ROI_Stats.csv",opts.outPref),overwriteStats);
        writetable(roiData, ROI_Tbl_Name, 'WriteRowNames', true, 'Delimiter', ',');
        resolvedPaths.ROI_Stats = ROI_Tbl_Name;
        if v; fprintf('Wrote ROI stats to %s\n',ROI_Tbl_Name); end
    end

    %%write manifest%%
    runTimestamp = string(datetime('now',"Format","yyyyMMdd_HHmmss"));
    manifest = struct;
    manifest.runTimestamp = runTimestamp;
    manifest.B0 = B0;
    manifest.hdr_File = hdr_file;
    manifest.sorted_ASE_Files = ase_files;
    manifest.taus_raw = taus_raw;
    manifest.averageReps = opts.averageReps;
    manifest.taus_full = taus_full;
    manifest.TEs_full = TEs_full;
    manifest.frameEcho = frameEcho;
    manifest.frameSrc = frameSrc;
    manifest.usePts_full = usePts_full;
    manifest.shortInds_good = shortInds_good;
    manifest.longInds_good = longInds_good;
    manifest.t_crit = t_crit;
    manifest.fitOpts = opts;
    manifest.paths = resolvedPaths;
    manifestPath = sprintf("%s_manifest.json", opts.outPref);
    fid = fopen(manifestPath, 'w');
    fprintf(fid, '%s', jsonencode(manifest,'PrettyPrint',true));
    fclose(fid);
    if v; fprintf('Wrote manifest to %s\n', manifestPath); end
    if v;fprintf('DONE!\n');end
end



function roiData = addQualityColumns(roiData, M0, nVoxReg, nVoxKept, qual, rawMaps)
% ROI-table columns describing the voxel exclusion: kept occupancy, excluded fraction (total
% and per criterion bit; bits overlap, so they need not sum to the total), and medians of the
% raw voxel maps over the kept voxels of each region.
    roiData.nVox_kept    = nVoxKept;
    roiData.fracExcluded = 1 - nVoxKept ./ nVoxReg;
    bi = qual.bitInfo;
    for i = 1:height(bi)
        if ~bi.Enabled(i), continue; end
        hit = double(bitand(qual.flags(:), uint16(2^bi.Bit(i))) > 0).';
        roiData.(sprintf('fracExcl_b%d_%s', bi.Bit(i), bi.Name{i})) = full(sum(M0 .* hit, 2)) ./ nVoxReg;
    end
    keepIdx = find(qual.keep(:));
    Mk = M0(:, keepIdx) > 0;                                % kept members of each region
    for f = {'OEF','R2P','vCBV','R2'}
        x = rawMaps.(f{1})(keepIdx);
        med = NaN(size(Mk,1),1);
        for r = 1:size(Mk,1)
            med(r) = median(x(Mk(r,:)), 'omitnan');
        end
        roiData.([f{1} '_med']) = med;
    end
end

function resolvedPaths = writeQualityMaps(qual, templateHdr, outPref, nEchos, overwrite, v, resolvedPaths)
% Write per-echo quality diagnostics (4th dim = echo) as <outPref>_<name>_e<k>,
% and single maps as <outPref>_<name>. Every file gets a manifest entry.
    if isfield(qual,'keep')
        qual.qualityMask  = double(qual.keep);     % 1 = kept
        qual.qualityFlags = double(qual.flags);    % uint16 bitmask, exact in single precision
    end
    perEcho = {'sigma','SNR','Gz','zshimAmp','zshimAtten','zshimEdge'};
    for i = 1:numel(perEcho)
        nm = perEcho{i};
        if ~isfield(qual,nm) || isempty(qual.(nm)), continue; end
        for e = 1:nEchos
            key = sprintf('%s_e%d', nm, e);
            resolvedPaths.(key) = writeMap(templateHdr, qual.(nm)(:,:,:,e), sprintf("%s_%s", outPref, key), ...
                'verbose',v, 'overwrite',overwrite, 'mapName',key);
        end
    end
    single = {'qualityMask','qualityFlags','csfFrac','zshimFracLow','vCBV_quad','vCBV_SE','vCBVquad_SE','vCBV_z','R2P_e1','R2P_e2','R2P_Fstat'};
    for i = 1:numel(single)
        nm = single{i};
        if ~isfield(qual,nm) || isempty(qual.(nm)), continue; end
        resolvedPaths.(nm) = writeMap(templateHdr, qual.(nm), sprintf("%s_%s", outPref, nm), ...
            'verbose',v, 'overwrite',overwrite, 'mapName',nm);
    end
end

function maps = applyPhysBounds(maps, physBds, clipVal)
% Apply [lo hi] bounds to each field. NaN stays NaN in both modes.
    nms = fieldnames(physBds);
    for i = 1:numel(nms)
        nm = nms{i};
        if ~isfield(maps,nm) || isempty(maps.(nm)) || isempty(physBds.(nm)), continue; end
        lo = physBds.(nm)(1);  hi = physBds.(nm)(2);
        M = maps.(nm);
        bad = isnan(M);
        if isempty(clipVal)
            M = min(max(M,lo),hi);
        else
            M(M < lo | M > hi) = clipVal;
        end
        M(bad) = NaN;                       % max(NaN,lo) returns lo, so restore NaNs
        maps.(nm) = M;
    end
end

function out = fitASERegimes(aseLong, aseShort, tausLong, TEsLong, tausShort, TEsShort, nEchos, oefScale, mask, wLong, wShort)
% Shared by the voxelwise and ROI paths. ase*: X x Y x Z x N x R (or nReg x 1 x 1 x N x R).
% taus/TEs in ms (column vectors). oefScale = (4/3)*pi*gamma*dChi_0*Hct*B0.
% wLong/wShort: optional fit weights, same size as aseLong/aseShort (empty = OLS).
    if nargin < 10, wLong = []; end
    if nargin < 11, wShort = []; end
    tauL = abs(tausLong(:))/1000;   TEL = TEsLong(:)/1000;
    tauS = abs(tausShort(:))/1000;  TES = TEsShort(:)/1000;
    %Long regime:
    % fitting: ln(S) = ln(C) + vCBV − R2*TE − R2P*|dTE|
    if nEchos == 1
        %unable to compute R2 so we just regress against dTE:
        XLdes = [ones(numel(tauL),1), -tauL];
        [BL,rmsL,nL] = fitMultiRegressor(aseLong, XLdes, mask, [], wLong);
        % BL(:,:,:,1) = ln(C) + vCBV − R2*TE
        R2  = [];
        R2P = BL(:,:,:,2);
        modShort = aseShort;                 % constant TE: R2*TE cancels in vCBV
    else
        %here we can (and need) to get R2 so we regress against TE AND dTE
        XLdes = [ones(numel(tauL),1), -TEL, -tauL];
        [BL,rmsL,nL] = fitMultiRegressor(aseLong, XLdes, mask, [], wLong);
        % B_long(:,:,:,1) = ln(C) + vCBV
        R2  = BL(:,:,:,2);
        R2P = BL(:,:,:,3);
        %   adding R2*TE(1) from long fit gives:
        %       ln(S)+R2*TE(1) = ln(C) - k*vCBV*dTE²
        modShort = aseShort .* exp(R2 .* reshape(TES,1,1,1,[]));   % ln(S) + R2*TE
    end
    %Short Regime:
    % fitting: ln(S) = ln(C) − R2·TE − k*vCBV*dTE²
    tau2 = tauS.^2;
    if max(tau2) - min(tau2) > 0
        Xs = [ones(numel(tauS),1), tau2];
    else
        Xs = ones(numel(tauS),1);            % all short frames at tau = 0: intercept only
    end
    [BS,rmsS,nS] = fitMultiRegressor(modShort, Xs, mask, [], wShort);

    vCBV = BL(:,:,:,1) - BS(:,:,:,1);
    OEF  = R2P ./ (vCBV * oefScale);
    OEF(vCBV <= 0 | R2P <= 0) = NaN;         % avoids negative/negative = spurious positive OEF

    % quadratic (short-regime curvature) vCBV: ln S = c - 0.3*vCBV*(dw*tau)^2 and R2' = vCBV*dw
    %   => k = -0.3*R2'^2/vCBV  =>  vCBV_quad = -0.3*R2'^2/k   (only defined for k < 0)
    if size(Xs,2) == 2
        k = BS(:,:,:,2);
        vCBV_quad = -0.3*R2P.^2 ./ k;
        vCBV_quad(~(k < 0)) = NaN;
    else
        k = [];  vCBV_quad = [];
    end

    out = struct('BL',BL,'rmsL',rmsL,'nValidL',nL, 'BS',BS,'rmsS',rmsS,'nValidS',nS, ...
                 'R2',R2,'R2P',R2P,'vCBV',vCBV,'OEF',OEF,'k',k,'vCBV_quad',vCBV_quad, ...
                 'XL',XLdes,'XS',Xs,'TES',TES);
end

function out = fitR2PerEcho(aseLong, tausLong, echoLong, rmsJoint, nJoint, wLong)
% Separate long-regime fits per echo (ln S = c_e - R2'_e*tau) and an F-test of the joint
% model (shared R2', TE term) against the split model (4 parameters).
%   F = (RSS_joint - RSS_split) / (RSS_split/(n-4)),  df = (1, n-4)
    tau = abs(tausLong(:))/1000;
    sz = size(aseLong,1:3);
    R2P = NaN([sz 2]);  rss = zeros(sz);  n = zeros(sz);
    for e = 1:2
        sel = echoLong(:) == e;
        if isempty(wLong), we = []; else, we = wLong(:,:,:,sel); end
        [B,r,ne] = fitMultiRegressor(aseLong(:,:,:,sel), [ones(nnz(sel),1), -tau(sel)], [], [], we);
        R2P(:,:,:,e) = B(:,:,:,2);
        rss = rss + r.^2.*ne;  n = n + ne;
    end
    rssJ = rmsJoint.^2 .* nJoint;
    F = (rssJ - rss) ./ (rss ./ (n - 4));
    F(n <= 4 | ~isfinite(F)) = NaN;
    F(F < 0) = 0;                           % nested models: negative only by round-off
    p = 1 - fcdf(F, 1, n - 4);
    out = struct('R2P_e1',R2P(:,:,:,1),'R2P_e2',R2P(:,:,:,2),'Fstat',F,'pval',p);
end


function [meanSig, nValid] = regionMeans(data, M)
% Mean over each region's voxels, per frame, ignoring NaN/Inf/<=0 values.
% data: X x Y x Z x N x R  ->  meanSig, nValid: nReg x 1 x 1 x N x R
    sz = size(data);  sz(end+1:5) = 1;
    nReg = size(M,1);  nCol = sz(4)*sz(5);
    D = reshape(data, [], nCol);                  % no copy
    sumSig = zeros(nReg,nCol);  nValid = zeros(nReg,nCol);
    for c = 1:16:nCol                             % column blocks keep memory low
        cols = c:min(c+15,nCol);
        X  = D(:,cols);
        ok = isfinite(X) & X > 0;
        X(~ok) = 0;
        sumSig(:,cols) = M*X;
        nValid(:,cols) = M*double(ok);
    end
    meanSig = reshape(sumSig ./ nValid, nReg,1,1,sz(4),sz(5));   % 0/0 -> NaN
    nValid  = reshape(nValid,           nReg,1,1,sz(4),sz(5));
end

function w = zshimWeights(Gz, taus, echo, dz)
% F^2 fit weights per voxel and frame for the z-shim-corrected data (same F as the correction:
% each frame uses its own echo's Gz and t = |tau|; voxels without a Gz estimate get 1).
% Gz: X x Y x Z x nEchos (uT/m); taus in ms; dz per echo (m).  w: X x Y x Z x N (single)
    sz = size(Gz,1:3);  N = numel(taus);
    w = ones([sz N], 'single');
    for f = 1:N
        F = zshimAttenuation(Gz(:,:,:,echo(f)), taus(f)/1000, dz(echo(f)));
        F(isnan(F)) = 1;
        w(:,:,:,f) = F.^2;
    end
end

function wR = regionWeights(data, w, M)
% Weight of each region-mean frame for the ROI fit. The mean of corrected signals has log-space
% variance proportional to sum(m^2/F^2)/sum(m)^2, against sum(m^2)/sum(m)^2 without attenuation,
% so its relative weight is sum(m^2)/sum(m^2/w) over the voxels valid in that frame (w = F^2).
% data, w: X x Y x Z x N x R  ->  wR: nReg x 1 x 1 x N x R
    sz = size(data);  sz(end+1:5) = 1;
    nReg = size(M,1);  nCol = sz(4)*sz(5);
    D  = reshape(data, [], nCol);
    Wv = reshape(w, [], nCol);
    M2 = M.^2;
    num = zeros(nReg,nCol);  den = zeros(nReg,nCol);
    for c = 1:16:nCol
        cols = c:min(c+15,nCol);
        ok = isfinite(D(:,cols)) & D(:,cols) > 0;
        iw = double(ok) ./ double(Wv(:,cols));
        iw(~ok) = 0;
        num(:,cols) = M2*double(ok);
        den(:,cols) = M2*iw;
    end
    wR = reshape(num ./ den, nReg,1,1,sz(4),sz(5));    % 0/0 -> NaN (frame unused anyway)
end
