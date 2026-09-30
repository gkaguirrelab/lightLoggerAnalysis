function [radianceMap, imageStages] = reconstructionPipeline(I, AGCSettings)
% Performs the complete forward reconstruction pipeline in a single
% integrated function, using persistent variables to load all necessary
% derived parameters upfront.

% Declare persistent variables for all derived parameters and maps
persistent clippingExponent darkSignal ...
    correctionMap radiometricCorrectionMap ...
    integratedRadiance cameraScore ...
    meanCorrectionFielding meanCorrectionRGB Smax

% Load non-linear clipping exponent and linearized set point
if isempty(clippingExponent)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName, 'clippingExponent');
end

% Load dark signal and compute Smax for linearization
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
    % Calculate the mean of the fielding correction map
    meanCorrectionFielding = mean(correctionMap(:), 'omitnan');
end

% Load RGB radiometric correction map and compute its mean scaling factor
if isempty(radiometricCorrectionMap)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radiometricCorrectionRGB.mat');
    load(paramFileName, 'radiometricCorrectionMap');
    meanCorrectionRGB = mean(radiometricCorrectionMap(:), 'omitnan');
end

% Load camera score to effective integrated radiance mapping parameters
if isempty(integratedRadiance)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToIntegratedRadiance.mat');
    load(paramFileName,'agcToRadianceP');
end


% Set the saturation threshold based upon the clippingExponent value

% Define the maximum allowable noise amplification (derivative) A value of
% 3.0 to 5.0 is typically a safe boundary for Bayesian conditioning
maxAllowedDerivative = 4.0; 

% Dynamically calculate the saturation threshold based on the derivative
exponentRatio = clippingExponent / (clippingExponent + 1);
yPrimeThresh = Smax * (1 - maxAllowedDerivative^(-exponentRatio))^(1/clippingExponent);

% Add dark signal to map back to absolute raw sensor counts
saturationThreshold = floor(yPrimeThresh + darkSignal);


% --- Processing Pipeline Stages ---

% Stage 1: Convert from uint8 to double float
imageStages{1} = double(I);

% Stage 2: Linearize sensor counts; set to Inf any values above the
% saturation threshold in the raw image
y = imageStages{1};
y(y < darkSignal) = darkSignal;
yPrime = y - darkSignal;
n = clippingExponent;
asymptoticGain = 1 ./ (1 - (yPrime ./ Smax).^n).^(1./n);
linearized = yPrime .* asymptoticGain;
linearized(y >= saturationThreshold) = Inf;
imageStages{2} = linearized;

% Apply digital gain
imageStages{2} = imageStages{2} * AGCSettings.Dgain;

% Stage 3: Impute values for ceiling and floor pixels
imageStages{3} = imputePixelValues(imageStages{2});

% Stage 4: Flat fielding correction
imageStages{4} = imageStages{3} .* correctionMap;

% Stage 5: Equalize RGB channels
imageStages{5} = imageStages{4} .* radiometricCorrectionMap;

% Obtain a linearized set point. This is the sensor value (after
% linearization) that corresponds to the set point that the AGC attempts to
% obtain for the mean of the entire image. To do so, we take the initial
% set point, undo Dgain effects, linearize, and account for the mean
% fielding and RGB corrections
setPoint = 127;
setPoint = (setPoint / AGCSettings.Dgain) - darkSignal;
linearizedSetPoint = setPoint ./ (1 - (setPoint ./ Smax).^n).^(1./n);
linearizedSetPoint = linearizedSetPoint * meanCorrectionFielding * meanCorrectionRGB;

% Stage 6: Convert to absolute radiance units
thisCameraScore = AGCSettings.exposure * AGCSettings.Again * AGCSettings.Dgain;

% Obtain the mean integrated radiance implied by this camera score
meanIntegratedRadiance = 10.^polyval(agcToRadianceP,log10(thisCameraScore))

% Scale the radiometrically balanced image to absolute radiance
imageStages{6} = (imageStages{5} / linearizedSetPoint) * meanIntegratedRadiance;

% Return the final stage as the radiance map
radianceMap = imageStages{6};

end