function I_raw = synthesisPipeline(radianceMap, AGCSettings)
% synthesisPipeline performs the inverse of reconstructionPipeline.
% Given a radiance map and camera AGCSettings, this function returns the
% original raw sensor image as a uint8 array.

% Declare persistent variables for derived parameters and correction maps
persistent clippingExponent linearizedSetPoint darkSignal ...
    correctionMap radiometricCorrectionMap ...
    avgSceneRadiance cameraScore ...
    meanCorrectionFielding meanCorrectionRGB Smax

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

% Load camera score to average radiance mapping parameters
if isempty(avgSceneRadiance)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToAverageRadiance.mat');
    load(paramFileName, 'avgSceneRadiance', 'cameraScore');
end

% Calculate the effective set point accounting for spatial scaling
effectiveSetPoint = linearizedSetPoint * meanCorrectionFielding * meanCorrectionRGB;

% Determine mean environmental radiance implied by the given AGCSettings
thisCameraScore = AGCSettings.exposure * AGCSettings.Again * AGCSettings.Dgain;
logMeanSceneRadiance = interp1(log10(cameraScore), log10(avgSceneRadiance), log10(thisCameraScore), 'linear');
meanSceneRadiance = 10.^logMeanSceneRadiance;

% 1. Inverse of Radiance Conversion (Stage 5 -> Stage 4)
I_sensor_corrected = (radianceMap / meanSceneRadiance) * effectiveSetPoint;

% 2. Inverse of RGB Radiometric Correction (Stage 4 -> Stage 3)
I_flat = I_sensor_corrected ./ radiometricCorrectionMap;

% 3. Inverse of Flat Fielding Correction (Stage 3 -> Stage 2)
yLinear = I_flat ./ correctionMap;

% 4. Inverse of Sensor Linearization (Stage 2 -> Stage 1 double)
n = clippingExponent;
yPrime = (yLinear * Smax) ./ (Smax.^n + yLinear.^n).^(1./n);
y = yPrime + darkSignal;

% 5. Clip, round, and convert back to uint8 raw sensor counts (Stage 1 -> uint8)
I_raw = uint8(round(max(0, min(255, y))));

end