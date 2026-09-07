function [radianceMap, imageStages] = reconstructionPipeline(I, AGCSettings)
% Performs the complete forward reconstruction pipeline in a single
% integrated function, using persistent variables to load all necessary
% derived parameters upfront.

% Declare persistent variables for all derived parameters and maps
persistent clippingExponent linearizedSetPoint darkSignal ...
    correctionMap radiometricCorrectionMap ...
    avgSceneRadiance cameraScore ...
    meanCorrectionFielding meanCorrectionRGB Smax

% 1. Load non-linear clipping exponent and linearized set point
if isempty(clippingExponent)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName, 'clippingExponent', 'linearizedSetPoint');
end

% 2. Load dark signal and compute Smax for linearization
if isempty(darkSignal)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'darkSignal.mat');
    load(paramFileName, 'darkSignal');
    Smax = 2^8 - 1 - darkSignal;
end

% 3. Load flat fielding correction map and compute its mean scaling factor
if isempty(correctionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'flatFieldingFunction.mat');
    load(paramFileName, 'correctionMap');
    meanCorrectionFielding = mean(correctionMap(:), 'omitnan');
end

% 4. Load RGB radiometric correction map and compute its mean scaling factor
if isempty(radiometricCorrectionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radiometricCorrectionRGB.mat');
    load(paramFileName, 'radiometricCorrectionMap');
    meanCorrectionRGB = mean(radiometricCorrectionMap(:), 'omitnan');
end

% 5. Load camera score to average radiance mapping parameters
if isempty(avgSceneRadiance)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToAverageRadiance.mat');
    load(paramFileName, 'avgSceneRadiance', 'cameraScore');
end

% --- Processing Pipeline Stages ---

% Stage 1: Convert from uint8 to double float
imageStages{1} = double(I);

% Stage 2: Linearize sensor counts
y = imageStages{1};
y(y < darkSignal) = darkSignal;
yPrime = y - darkSignal;
n = clippingExponent;
imageStages{2} = yPrime ./ (1 - (yPrime ./ Smax).^n).^(1./n);

% Stage 3: Flat fielding correction
imageStages{3} = imageStages{2} .* correctionMap;

% Stage 4: Radiometric correction RGB
imageStages{4} = imageStages{3} .* radiometricCorrectionMap;

% Stage 5: Convert to absolute radiance units
effectiveSetPoint = linearizedSetPoint * meanCorrectionFielding * meanCorrectionRGB;
thisCameraScore = AGCSettings.exposure * AGCSettings.Again * AGCSettings.Dgain;
logMeanSceneRadiance = interp1(log10(cameraScore), log10(avgSceneRadiance), log10(thisCameraScore), 'linear');
meanSceneRadiance = 10.^logMeanSceneRadiance;
imageStages{5} = (imageStages{4} / effectiveSetPoint) * meanSceneRadiance;

% Return the final stage as the radiance map
radianceMap = imageStages{5};

end