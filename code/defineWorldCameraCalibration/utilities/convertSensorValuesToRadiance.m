function I_radiance = convertSensorValuesToRadiance(I_sensor,AGCSettings)
% This function converts an input image from sensor values to absolute
% radiance values. The input image must have already undergone several
% stages of correction: linearization, flat fielding, and radiometric
% correction of the channels.

% The "linearizedSetPoint" is the sensor value that corresponds to the set
% point of the gain control of the camera in linearized values. That is,
% the camera AGC attempts to keep the mean of the image at a sensor value
% of 127. The linearization step maps 127 --> 57.
persistent linearizedSetPoint
if isempty(linearizedSetPoint)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'nonLinearClippingExponent.mat');
    load(paramFileName,'linearizedSetPoint');
end

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
meanSceneRadiance = interp1(cameraScore, avgSceneRadiance, thisCameraScore, 'linear');

% Perform the conversion
I_radiance = (I_sensor / linearizedSetPoint) * meanSceneRadiance;

end
