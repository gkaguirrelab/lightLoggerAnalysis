function [I_raw, radianceMap] = synthesisPipeline(radianceModel, radianceModelS, AGCSettings)
% synthesisPipeline performs the inverse of reconstructionPipeline.
% Given a radiance model, its sampling structure, and camera AGCSettings, 
% this function returns the original raw sensor image as a uint8 array.

% Declare persistent variables for derived parameters and correction maps
persistent clippingExponent linearizedSetPoint darkSignal ...
    correctionMap radiometricCorrectionMap ...
    effectiveRadiance cameraScore ...
    meanCorrectionFielding meanCorrectionRGB Smax ...
    azimuthMap elevationMap T channelNames bayerPattern

% Load non-linear clipping exponent and linearized set point
if isempty(clippingExponent)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName, 'clippingExponent', 'linearizedSetPoint');
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

% Load flat fielding correction map and compute its mean
if isempty(correctionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'flatFieldingFunction.mat');
    load(paramFileName, 'correctionMap');
    meanCorrectionFielding = mean(correctionMap(:), 'omitnan');
end

% Load RGB radiometric correction map and compute its mean
if isempty(radiometricCorrectionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radiometricCorrectionRGB.mat');
    load(paramFileName, 'radiometricCorrectionMap');
    meanCorrectionRGB = mean(radiometricCorrectionMap(:), 'omitnan');
end

% Load camera score to effective integrated radiance mapping parameters
if isempty(effectiveRadiance)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToEffectiveRadiance.mat');
    load(paramFileName, 'effectiveRadiance', 'cameraScore');
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

% 0. Compute Spatially-Varying Channel Radiance Map from radianceModel
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

% Calculate the effective set point accounting for spatial scaling
effectiveSetPoint = linearizedSetPoint * meanCorrectionFielding * meanCorrectionRGB;

% Determine mean effective radiance implied by the given AGCSettings
thisCameraScore = AGCSettings.exposure * AGCSettings.Again * AGCSettings.Dgain;
logThisEffectiveRadiance = interp1(log10(cameraScore), log10(effectiveRadiance), log10(thisCameraScore), 'linear');
thisEffectiveRadiance = 10.^logThisEffectiveRadiance;

% Calculate the Bayer-weighted mean effective radiance (1 Red, 2 Green, 1 Blue)
meanEffectiveRadiance = (thisEffectiveRadiance(1) + 2*thisEffectiveRadiance(2) + thisEffectiveRadiance(3)) / 4;

% 1. Inverse of Radiance Conversion (Stage 6 -> Stage 5)
I_sensor_corrected = (radianceMap / meanEffectiveRadiance) * effectiveSetPoint;

% 2. Inverse of RGB Radiometric Correction (Stage 5 -> Stage 4)
I_flat = I_sensor_corrected ./ radiometricCorrectionMap;

% 3. Inverse of Flat Fielding Correction (Stage 4 -> Stage 3/2)
yLinear = I_flat ./ correctionMap;

% 4. Inverse of Sensor Linearization (Stage 2 -> Stage 1 double)
% Note: The inverse of Stage 3 (imputing ceiling/floor pixels) is skipped 
% because absolute pixel loss at saturation boundaries cannot be deterministically inverted.
n = clippingExponent;
yPrime = (yLinear * Smax) ./ (Smax.^n + yLinear.^n).^(1./n);
y = yPrime + darkSignal;

% 5. Clip, round, and convert back to uint8 raw sensor counts (Stage 1 -> uint8)
I_raw = uint8(round(max(0, min(255, y))));

end