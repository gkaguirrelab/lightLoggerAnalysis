function channelEnergy = calculateIMX219ChannelRadiantEnergyViaMS(radianceModel, radianceModelS)
% COMPUTEUNIFORMFIELDCHANNELENERGY Computes the IMX219 channel radiant energies
% estimated via the ASM7341 minispect simulation and Tikhonov radiance reconstruction.

persistent T bayerIdx channelNames endIdx imxWls deltaSteradians
if isempty(T)
    projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');
    
    % 1. Load IMX219 spectral sensitivities and setup Bayer indices
    dataFileName = fullfile(projectRoot, 'data', 'IMX219_spectralSensitivity.mat');
    load(dataFileName, 'T');
    endIdx = find(T.wls == 780);
    imxWls = T.wls(1:endIdx);
    channelNames = {'red', 'green', 'blue'};
    bayerPattern = "BGGR";
    
    rows = 480;
    columns = 640;
    [bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(zeros(rows, columns), bayerPattern);
    
    % 2. Load precomputed per-pixel solid angles
    deltaPath = fullfile(projectRoot, 'derived', 'deltaSteradians.mat');
    load(deltaPath, 'deltaSteradians');
end

% Get wavelengths from radianceModelS
modelWls = SToWls(radianceModelS);
nmStep = modelWls(2) - modelWls(1);

% 1. Model expected ASM7341 sensor values from the radiance model[cite: 4]
[minispectValues, ~, ~] = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS);

% 2. Estimate mean environmental radiance spectrum from minispect values[cite: 5]
[estimatedSpectrumRaw, estimatedS] = estimateRadianceSpectrumFromMinispect(minispectValues);
estWls = SToWls(estimatedS);

% Interpolate estimated spectrum onto the IMX219 wavelength grid
estimatedSpectrum = interp1(estWls, estimatedSpectrumRaw, imxWls, 'linear', 0);

% 3. Predict IMX219 channel radiant energies (Watts) from the estimated spectrum
channelEnergy = struct();
for cc = 1:3
    thisSensitivity = T.(channelNames{cc})(1:endIdx);
    thisSensitivityNormed = thisSensitivity ./ max(thisSensitivity);
    
    channelRadianceEst = (thisSensitivityNormed' * estimatedSpectrum) * nmStep;
    channelSolidAngleSum = sum(deltaSteradians(bayerIdx{cc}));
    channelEnergy.(channelNames{cc}) = channelRadianceEst * channelSolidAngleSum;
end

end