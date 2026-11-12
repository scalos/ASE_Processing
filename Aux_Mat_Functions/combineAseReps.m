function aseData = combineAseReps(aseFiles,hdrInfo, opts)
% combineAseReps  Combine per-echo ASE files into one 4D volume.
%   aseFiles : cell (nEchos x 1), one processed 4D file per echo, ordered by echo
%   hdrInfo  : info struct (4th output of readASEHeader)
%   opts.method:
%     'avg'   -> X x Y x Z x (nSweep*nEchos + nZshim*nEchos), single
%                [e1 sweep avg | e2 sweep avg | ... | e1 zshim | e2 zshim | ...]
%                Each sweep frame is averaged across reps (same position in each rep);
%                z-shim frames are copied as-is.
%     'stack' -> X x Y x Z x (ndTE*nEchos), native class
%                [all e1 frames | all e2 frames | ...]
%   Frame order within each block is acquisition order (as listed in the header).
%
%   Each file's frames are expected as: nReps sweeps of nSweep frames, then nZshim frames.
    arguments
        aseFiles cell % nEchos
        hdrInfo struct
        opts.method {mustBeMember(opts.method,{'avg','stack'})} = 'avg';
    end

    nEchos  = hdrInfo.nEchos;
    nReps   = hdrInfo.nReps;
    nSweep  = hdrInfo.nSweep;
    nZshim  = hdrInfo.nZshimming;
    repInds = hdrInfo.repInds;   % nReps x nSweep frame indices
    ndTE    = nSweep*nReps + nZshim;

    assert(numel(aseFiles)==nEchos, ...
        'Expected one ase file per echo (%d). Got %d.', nEchos, numel(aseFiles));
    isequalNiftiGeometry(aseFiles(:),"onMismatch","error");

    % frame counts from headers, so we fail before loading any data
    for eIdx = 1:nEchos
        fInfo = niftiinfo(aseFiles{eIdx});
        nt = 1;
        if numel(fInfo.ImageSize) >= 4; nt = fInfo.ImageSize(4); end
        assert(nt==ndTE, ...
            ['%s has %d frames but header implies %d ' ...
             '(%d reps x %d sweep + %d z-shim).'], ...
            aseFiles{eIdx}, nt, ndTE, nReps, nSweep, nZshim);
    end
    sz = fInfo.ImageSize(1:3);
    zShimInds = nSweep*nReps + (1:nZshim);

    switch opts.method
        case 'avg'
            aseData = zeros([sz, nSweep*nEchos + nZshim*nEchos], 'single');
            zOff = nSweep*nEchos;   % z-shim blocks follow all sweep blocks
            for eIdx = 1:nEchos
                echoData = single(niftiread(aseFiles{eIdx}));
                sweepSum = zeros([sz, nSweep], 'single');
                for rIdx = 1:nReps
                    sweepSum = sweepSum + echoData(:,:,:,repInds(rIdx,:));
                end
                aseData(:,:,:,nSweep*(eIdx-1) + (1:nSweep)) = sweepSum./nReps;
                aseData(:,:,:,zOff + nZshim*(eIdx-1) + (1:nZshim)) = echoData(:,:,:,zShimInds);
            end

        case 'stack'
            echoData = cell(1,nEchos);
            for eIdx = 1:nEchos
                echoData{eIdx} = niftiread(aseFiles{eIdx});
            end
            aseData = cat(4, echoData{:});
    end
end
