function cMat = FitShortRegime(ydata, xdata)
    % Returns coefficient map where cMat(:,:,:,1) is C (additive),
    % cMat(:,:,:,2) is A (multiplicative), fitting Y = A*X^2 + C,
    % and cMat(:,:,:,3) is the RMS residual (log space)
    [Np,Nf,Ns,Nfr] = size(ydata);
    ydata(ydata <= 0) = NaN;
    ydata = log(ydata);
    cMat = zeros(Np,Nf,Ns, 3);
    % Design matrix: column of ones for the constant term, then x^2
    X = [ones(numel(xdata),1), xdata(:).^2];
    for sl = 1:size(ydata,3)
        series = squeeze(ydata(:,:,sl,:));
        data_ = reshape(series, Np*Nf, Nfr)';   % Nfr x (Np*Nf)
        coeffs = NaN(Np*Nf, 3);
        for v = 1:(Np*Nf)
            y = data_(:,v);
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