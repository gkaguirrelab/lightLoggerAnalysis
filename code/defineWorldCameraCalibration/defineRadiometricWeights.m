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
    'radiometricCorrectionRGB',...
    'cloudySkySPD_37degSolarElevation.mat');
load(dataFileName,'radiance','wls');

% Load the channel spectral sensitivity functions. This is a table with the
% first column providing the wavelength support.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');

% Identify the location of the IMX219 image(s) of the cloudy sky.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'radiometricCorrectionRGB',...
    'rawFrames',...
    '*.tiff');
fileSet = dir(dataFileName);

% Define a crop region from  these images that includes just the sky.
cropRegionLRTB = [211,410,151,300];
cropMask = zeros(480,640);
cropMask(cropRegionLRTB(3):cropRegionLRTB(4),cropRegionLRTB(1):cropRegionLRTB(2))=1;

% Prepare to linearize the raw sensor values. To do so, we load the the
% derived, nonLinearClippingExponent, and then call the linearization
% function
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
load(paramFileName,'clippingExponent');

% Prepare to correct for the fielding function
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'flatFieldingFunction.mat');
load(paramFileName,'correctionMap');

% What is the Bayer pattern in these data?
bayerPattern = "BGGR";

% Loop through the images
pixelValsRGB = {[],[],[]};
for ii = 1:length(fileSet)

    % Load the image
    fileName = fullfile(fileSet(ii).folder,fileSet(ii).name);
    I = double(imread(fileName));    

    % Retain only the first channel, as this is a raw image and the other
    % channels are simply a copy of this one. CAN REMOVE THIS AFTER ZACH
    % SAVES THESE FILES AS GRAYSCALE
    I = I(:,:,1);

    % Linearize
    linearI = linearizeIMX219SensorCounts(I,clippingExponent);

    % Correct for the fielding function
    linearFlatI = linearI .* correctionMap;

    % Set values outside the crop zone to nan
    outside = ~cropMask;
    linearFlatI(outside)=nan;

    % Obtain the R, G, and B components
    imageVals = {};
    [imageVals{1}, imageVals{2}, imageVals{3}] = ...
        splitBayerFrame(linearFlatI, bayerPattern);

    % Store the non-nan R, G, and B pixel values in the growing array
    for cc=1:3
        theseVals = imageVals{cc};
        theseVals = theseVals(~isnan(theseVals));
        pixelValsRGB{cc} = [pixelValsRGB{cc}; theseVals(:)];
    end
end

% Take the mean of the pixel values within each channel (R,G,B) across
% images and pixels within the region of interest
sensorValues = cellfun(@(x) mean(x),pixelValsRGB);

% Normalize the sensor values to the blue channel
sensorValuesNormed = sensorValues / sensorValues(3);

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
radiometricCorrectionRGB = predictedSensorValueNormed ./  sensorValuesNormed;

% We further adjust this triplet so that the mean sensor value (across RGB)
% is unchanged by this operation
k = 3/sum(radiometricCorrectionRGB);
radiometricCorrectionRGB = radiometricCorrectionRGB * k;

% Report the correction to the console
fprintf('The radiometric correction tuple (RGB) is: [%2.2f, %2.2f, %2.2f]\n',radiometricCorrectionRGB);

% Save the radiometric correction to the "derived" directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
readme = ['Created by defineRadiometricWeights.\n'...
    'radiometricCorrectionRGB -- multiple the (linearized) sensor values by these correction factors.\n'];
save(saveFileName,'readme','radiometricCorrectionRGB');
