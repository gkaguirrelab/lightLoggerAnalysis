function I_radiance = convertSensorValuesToRadiance(I_sensor,AGCSettings)
% This function converts an input image from sensor values to absolute
% radiance values. The input image must have already undergone several
% stages of correction: linearization, flat fielding, and radiometric
% correction of the channels.

% The "linearizedSetPoint" is the sensor value that corresponds to the set
% point of the gain control of the camera in linearized values.
persistent linearizedSetPoint
if isempty(linearizedSetPoint)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName,'linearizedSetPoint');
end

% Load the correction map used for flat fielding to determine the 
% scale factor that altered the mean of the image array.
persistent meanCorrection
if isempty(meanCorrection)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'flatFieldingFunction.mat');
    load(paramFileName,'correctionMap');
    % Use omitnan in case the spatial map includes masked boundaries
    meanCorrection = mean(correctionMap(:), 'omitnan');
end

% Adjust the expected set point to account for the overall spatial 
% scaling caused by the flat fielding normalization step.
effectiveSetPoint = linearizedSetPoint * meanCorrection;

% Load the mapping of camera AGCSettings to mean environmental radiance
persistent avgSceneRadiance cameraScore
if isempty(avgSceneRadiance)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraScoreToAverageRadiance.mat');
    load(paramFileName,'avgSceneRadiance','cameraScore');
end

% Obtain the mean environmental radiance implied by the AGC settings
thisCameraScore = AGCSettings.exposure * AGCSettings.Again * AGCSettings.Dgain;

% Perform linear interpolation in log10 space to accurately map the 
% inverse exponential relationship between AGC settings and radiance
logMeanSceneRadiance = interp1(log10(cameraScore), log10(avgSceneRadiance), log10(thisCameraScore), 'linear');
meanSceneRadiance = 10.^logMeanSceneRadiance;

% Perform the conversion using the geometrically adjusted set point
I_radiance = (I_sensor / effectiveSetPoint) * meanSceneRadiance;

end