function channelFluxDensity = estimateIMX219FluxFromMinispect(radianceModel, radianceModelS)
% ESTIMATEIMX219FLUXFROMMINISPECT Computes the estimated radiant flux density 
% in the R, G, and B channels of the IMX219 camera by internally simulating 
% the minispect values from the radianceModel and reconstructing the mean 
% spectral radiance, accounting for Bayer pattern channel pixel fractions.
%
%   channelFluxDensity = estimateIMX219FluxFromMinispect(radianceModel, radianceModelS)

    % 1. Internally compute simulated minispect values from the radiance model
    [minispectValues, spectralRadiance, spectralRadianceS] = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS);

    % 2. Reconstruct the mean spectral radiance vector from the minispect values
    %[spectralRadiance, spectralRadianceS] = estimateRadianceSpectrumFromMinispect(minispectValues);

    % 3. Load the IMX219 sensitivity functions using the correct structure
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(dataFileName, 'T');
    wlsSensor = T.wls;
    channelNames = {'red', 'green', 'blue'};

    % Extract and max-normalize sensitivities
    imxSensitivities = table2array(T(:, channelNames))'; % [3 x nWls]
    imxSensitivities = imxSensitivities ./ max(imxSensitivities, [], 2);

    % 4. Get wavelengths and bin width from reconstructed spectral radiance
    targetWls = SToWls(spectralRadianceS);
    deltaWl = spectralRadianceS(2); % Wavelength band sampling interval (e.g., 2 nm)

    % 5. Resample IMX219 sensitivities to match spectralRadiance wavelength sampling
    nChannels = length(channelNames);
    A_imx = nan(nChannels, length(targetWls));
    for i = 1:nChannels
        A_imx(i, :) = interp1(wlsSensor, imxSensitivities(i, :), targetWls, 'linear', 0);
    end

    % 6. Load radiometric correction map from derived folder matching pipeline structure
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radiometricCorrectionRGB.mat');
    load(paramFileName, 'radiometricCorrectionMap');
    meanCorrectionRGB = mean(radiometricCorrectionMap(:), 'omitnan');

    % 7. Compute integrated radiance (dot product with sensitivity and wavelength step)
    integratedFlux = (A_imx * spectralRadiance(:)) * deltaWl;
    
    % Apply radiometric correction scaling factor
    adjustedFlux = integratedFlux .* meanCorrectionRGB;

    % 8. Account for the fraction of the pixel array covered by each channel (RGGB Bayer pattern)
    % Red: 1/4 (0.25), Green: 2/4 (0.5), Blue: 1/4 (0.25)
    pixelFractions = [0.25; 0.5; 0.25];
    finalFlux = adjustedFlux .* pixelFractions;

    % 9. Package into a structured output
    channelFluxDensity = struct(...
        'Red', finalFlux(1), ...
        'Green', finalFlux(2), ...
        'Blue', finalFlux(3) ...
    );
end