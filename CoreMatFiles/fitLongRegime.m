function cMat = fitLongRegime(ydata, xdata)
    % Returns coefficient map where cMat(:,:,:,1) is y-int,
    % cMat(:,:,:,2) is the slope, and cMat(:,:,:,3) is the RMS residual (log space)
    [Np,Nf,Ns,Nfr] = size(ydata);
    cMat = zeros(Np,Nf,Ns, 3);
    X = [ones(numel(xdata),1), xdata(:)];   
    for sl = 1:size(ydata,3)
        series = squeeze(ydata(:,:,sl,:));  
        data_ = reshape(series, Np*Nf, Nfr)';
        data_(data_ <= 0) = NaN;
        logData = log(data_);
        coeffs = NaN(Np*Nf, 3);
        for v = 1:(Np*Nf)
            y = logData(:,v);
            mask = ~isnan(y);
            % Require at least 2 points
            if nnz(mask) >= 2
                beta = X(mask,:) \ y(mask);
                resid = X(mask,:)*beta - y(mask);
                rmsVal = sqrt(mean(resid.^2));
                coeffs(v,:) = [beta.', rmsVal];
            end
        end
        cMat(:,:,sl,:) = reshape(coeffs, Np, Nf, 3);
    end
end