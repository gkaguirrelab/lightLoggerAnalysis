function [I_raw, radianceMap] = synthesisPipeline(radianceModel, radianceModelS, AGCSettings)
% synthesisPipeline performs the inverse of reconstructionPipeline.
% Given a radiance model, its sampling structure, and camera AGCSettings, 
% this function returns the original raw sensor image as a uint8 array.

% Declare persistent variables for derived parameters and correction maps
persistent clippingExponent darkSignal ...
    correctionMap radiometricCorrectionMap ...
    agcToRadianceP ...
    meanCorrectionFielding meanCorrectionRGB Smax ...
    azimuthMap elevationMap T channelNames bayerPattern

% Load non-linear clipping exponent
if isempty(clippingExponent)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName, 'clippingExponent');
end

% Load dark signal and compute Smax
if isempty(darkSignal)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'darkSignal.mat');
    load(paramFileName, 'darkSignal');
    Smax = 2^8 - 1 - darkSignal;
end

% Load flat fielding correction map and compute its mean scaling factor
if isempty(correctionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'flatFieldingFunction.mat');
    load(paramFileName, 'correctionMap');
    % Calculate the mean using the harmonic-like formulation
    meanCorrectionFielding = 1 / mean(1 ./ correctionMap(:), 'omitnan');
end

% Load RGB radiometric correction map and compute its mean scaling factor
if isempty(radiometricCorrectionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radiometricCorrectionRGB.mat');
    load(paramFileName, 'radiometricCorrectionMap');
    meanCorrectionRGB = 1 / mean(1 ./ radiometricCorrectionMap(:), 'omitnan');
end

% Load camera score to effective integrated radiance mapping parameters
if isempty(agcToRadianceP)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToIntegratedRadiance.mat');
    load(paramFileName, 'agcToRadianceP');
end

% Load camera intrinsics and compute azimuth/elevation maps
if isempty(azimuthMap)
    projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');
    intrinsicsPath = fullfile(projectRoot, 'derived', 'arducamB0392cameraIntrinsics.mat');
    
    intrinsicsData = load(intrinsicsPath, 'arducamB0392cameraIntrinsics');
    fisheyeIntrinsics = intrinsicsData.arducamB0392cameraIntrinsics.results.Intrinsics;
    
    % Define the 640x480 sensor grid
    rows = 480;
    columns = 640;
    [xCoordinates, yCoordinates] = meshgrid(1:columns, 1:rows);
    sensorPoints = [xCoordinates(:), yCoordinates(:)];
    
    % Convert pixel centers to visual angles using the calibrated intrinsics
    visualAngles = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    visualAngles = reshape(visualAngles, rows, columns, 2);
    
    azimuthMap = deg2rad(visualAngles(:, :, 1));
    elevationMap = deg2rad(visualAngles(:, :, 2));
end

% Load IMX219 spectral sensitivities and setup channel info
if isempty(T)
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(dataFileName, 'T');
    channelNames = {'red', 'green', 'blue'};
    bayerPattern = "BGGR";
end

% Compute Spatially-Varying Channel Radiance Map from radianceModel
spectralRadianceMap = radianceModel(azimuthMap, elevationMap);
numWls = size(spectralRadianceMap, 1);
spectralRadianceFlat = reshape(spectralRadianceMap, numWls, []);

rows = size(azimuthMap, 1);
columns = size(azimuthMap, 2);
radianceMap = zeros(rows, columns);

% Get Bayer channel indices
[bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(radianceMap, bayerPattern);

% Get model wavelengths from radianceModelS
modelWls = SToWls(radianceModelS);

% Project spectral radiance onto each max-normalized Bayer channel sensitivity function
for cc = 1:3
    % Interpolate tabular sensitivity onto model wavelengths and zero-pad outside bounds
    thisSensitivity = interp1(T.wls, T.(channelNames{cc}), modelWls, 'linear', 0);
    thisSensitivityNormed = thisSensitivity ./ max(thisSensitivity);
    
    % Compute dot product across wavelengths and scale by wavelength bin width
    channelValsFlat = (thisSensitivityNormed' * spectralRadianceFlat) * radianceModelS(2);
    
    % Place into the corresponding Bayer locations in radianceMap
    tempMap = zeros(rows, columns);
    tempMap(bayerIdx{cc}) = channelValsFlat(bayerIdx{cc});
    radianceMap = radianceMap + tempMap;
end

% Calculate the linearized set point for the current AGC settings
n = clippingExponent;
setPoint = 127;
setPoint = (setPoint / AGCSettings.Dgain) - darkSignal;
linearizedSetPoint = setPoint ./ (1 - (setPoint ./ Smax).^n).^(1./n);
linearizedSetPoint = linearizedSetPoint * meanCorrectionFielding * meanCorrectionRGB;

% Determine mean integrated radiance implied by the given AGCSettings
thisCameraScore = (AGCSettings.exposure * AGCSettings.Again) / linearizedSetPoint;
meanIntegratedRadiance = 10.^polyval(agcToRadianceP, log10(thisCameraScore));

% Inverse of Radiance Conversion
I_sensor_corrected = (radianceMap / meanIntegratedRadiance) * linearizedSetPoint;

% Inverse of RGB Radiometric Correction
I_flat = I_sensor_corrected ./ radiometricCorrectionMap;

% Inverse of Flat Fielding Correction
yLinear = I_flat ./ correctionMap;

% Inverse of Sensor Linearization (handling the max(0, ...) logic from the forward pipeline)
yPrime = zeros(size(yLinear));
posIdx = yLinear >= 0;
negIdx = yLinear < 0;

% For positive values, apply the algebraic inverse of the asymptotic gain
yPrime(posIdx) = (yLinear(posIdx) .* Smax) ./ (Smax.^n + yLinear(posIdx).^n).^(1./n);

% For negative values, the forward pipeline applied an asymptotic gain of 1
yPrime(negIdx) = yLinear(negIdx);

% Restore dark signal
y = yPrime + darkSignal;

% Clip, round, and convert back to uint8 raw sensor counts
I_raw = uint8(round(max(0, min(255, y))));

end