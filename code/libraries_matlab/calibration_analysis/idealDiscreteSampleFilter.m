function filterProfile = idealDiscreteSampleFilter(sourceFreqsHz,dTsignal)
    nCycles = 100;
    dTsource = 0.001; % seconds

    for ii = 1:length(sourceFreqsHz)
        % Define the signal length
        sourceDurSecs = nCycles/sourceFreqsHz(ii);
        sourceTime = 0:dTsource:sourceDurSecs-dTsource;
        % Create a source modulation
        source = sin(sourceTime/(1/sourceFreqsHz(ii))*2*pi);
        % Downsample the source to create the signal
        signalTime = 0:dTsignal:sourceDurSecs-dTsignal;
        signal = interp1(sourceTime,source,signalTime,'linear');
        % Set up the regression matrix
        X = [];
        X(:,1) = sin(  sourceTime./(1/sourceFreqsHz(ii)).*2*pi );
        X(:,2) = cos(  sourceTime./(1/sourceFreqsHz(ii)).*2*pi );
        % Perform the fit
        y = interp1(signalTime,signal,sourceTime,'nearest','extrap')';
        b = X\y;
        filterProfile(ii)  = norm(b);
    end
end