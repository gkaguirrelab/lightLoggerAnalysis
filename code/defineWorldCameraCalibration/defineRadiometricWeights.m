% This script derives the multiplicative adjustments that should be applied
% to the R, G, and B channels so that the sensor values reflect the
% radiometric power of the light source.
%
% The PR670 was used to measure the SPD of a cloudy sky, at the same time
% that the IMX219 camera was used to collect images of the sky. We
% separately obtained the tabular spectral sensitivity functions of the R,
% G, and B channels of the camera chip.
%
% These measurements are processed below to yield the needed corrections.

% Housekeeping
clear
close all

% Load the SPD of the cloudy sky
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'cloudySkySPD_37degSolarElevation.mat');
load(dataFileName,'radiance','wls');

% Load the channel spectral sensitivity functions. This is a table with the
% first column providing the wavelength support.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');

% Load the IMX219 image(s) of the cloudy sky
% FIND THE DATA. WHAT WE WANT IS THE RAW IMAGES FROM THIS MEASUREMENT

% Process the data to extract the R, G, and B sensor values from within the
% region of the image that contains the cloudy sky. 

% CURRENTLY, THESE ARE THE VALUES ZACH POSTED TO SLACK
sensorValues = [146.07803403 223.66871456 177.44678639];

% Linearize the raw sensor values. To do so, we load the the derived,
% nonLinearClippingExponent, and then call the linearization function
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
load(paramFileName,'clippingExponent');
sensorValuesLinear = linearizeIMX219SensorCounts(sensorValues,clippingExponent);

% Normalize the sensor values to the blue channel
sensorValuesLinearNormed = sensorValuesLinear / sensorValuesLinear(3);

% Obtain the relative predicted channel values based upon the cloudy sky
% SPD. We first need to match up the wavelength support of the two
% measurements. This is hard-coded at the moment, given that the data
% inputs are known fixed.
endIdx = find(T.wls == 780);
assert(wls(1) == T.wls(1));
assert(wls(end) == T.wls(endIdx));
assert(diff(wls(1:2)) == diff(T.wls(1:2)));

% Set up a figure
figure
plot(wls,radiance,'.k');
hold on
xlabel('wavelength [nm]');
ylabel('radiance [w/m^2/sr/nm]');

% Loop through the channels
channelNames = {'red','green','blue'};
for ii = 1:3
    thisChannelSensitivity = T.(channelNames{ii})(1:endIdx);
    thisChannelSensitivityNormed = thisChannelSensitivity ./ max(thisChannelSensitivity);
    predictedSensorValue(ii) = radiance' * thisChannelSensitivityNormed;
    plot(T.wls(1:endIdx),thisChannelSensitivityNormed*predictedSensorValue(ii)/100,['-' channelNames{ii}(1)]);
end

% Normalize the predicted sensor values to the blue channel
predictedSensorValueNormed = predictedSensorValue / predictedSensorValue(3);

% Calculate the radiometric correction that must be applied to the observed
% sensor values to have them match the predicted sensor values
radiometricCorrectionRGB = predictedSensorValueNormed ./  sensorValuesLinearNormed;

% We further adjust this triplet so that the mean sensor value (across RGB)
% is unchanged by this operation
k = 3/sum(radiometricCorrectionRGB);
radiometricCorrectionRGB = radiometricCorrectionRGB * k;

% Save the radiometric correction to the "derived" directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
readme = ['Created by defineRadiometricWeights.\n'...
    'radiometricCorrectionRGB -- multiple the (linearized) sensor values by these correction factors.\n'];
save(saveFileName,'readme','radiometricCorrectionRGB');
