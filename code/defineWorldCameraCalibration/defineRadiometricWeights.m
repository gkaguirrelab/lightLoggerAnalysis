% This script derives the multiplicative adjustments that should be applied
% to the R, G, and B channels so that the sensor values reflect the
% radiometric power of the light source, correcting for the differing peak 
% transmission efficiencies of the Bayer filters.
%
% The PR670 was used to measure the SPD of a cloudy sky, at the same time
% that the IMX219 camera was used to collect images of the sky. We
% separately obtained the tabular spectral sensitivity functions of the R,
% G, and B channels of the camera chip.
%
% These measurements are processed below to yield the needed corrections.

% Housekeeping
clear

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

% Define a crop region from these images that includes just the sky.
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

    % Linearize
    linearI = linearizeIMX219SensorCounts(I,clippingExponent);

    % Correct for the fielding function
    linearFlatI = linearI .* correctionMap;

    % Set values outside the crop zone to nan
    outside = ~cropMask;
    linearFlatI(outside)=nan;

    % Obtain the R, G, and B components
    imageVals = {};
    [rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = ...
        returnBayerIndices(linearFlatI, bayerPattern);

% Store the valid, non-saturated R, G, and B pixel values in the growing array
    for cc=1:3
        theseVals = linearFlatI(rgbIdx{cc});
        
        % Filter out both NaNs (outside crop region) and Infs (saturated pixels)
        theseVals = theseVals(isfinite(theseVals));
        
        pixelValsRGB{cc} = [pixelValsRGB{cc}; theseVals(:)];
    end
end

% Take the mean of the pixel values within each channel (R,G,B) across
% images and pixels within the region of interest
sensorValues = cellfun(@(x) mean(x),pixelValsRGB);

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
title('SPD of cloudy sky and weighted channel sensitivity functions');

% Loop through the channels
channelNames = {'red','green','blue'};
predictedSensorValue = zeros(1,3);
for ii = 1:3
    thisChannelSensitivity = T.(channelNames{ii})(1:endIdx);
    % We use the max-normalized curves as requested
    thisChannelSensitivityNormed = thisChannelSensitivity ./ max(thisChannelSensitivity);
    predictedSensorValue(ii) = radiance' * thisChannelSensitivityNormed;
    plot(T.wls(1:endIdx),thisChannelSensitivityNormed*predictedSensorValue(ii)/100,['-' channelNames{ii}(1)]);
end

% Calculate raw mapping weights. These map the camera counts to the theoretical 
% max-normalized radiance. These weights embed the AGC state of the sky photo.
rawWeights = predictedSensorValue ./ sensorValues;

% Anchor the weights to the Green channel (index 2). This explicitly factors out 
% the sky photo's AGC scalar, leaving exactly the ratio of the Bayer filters' 
% peak quantum efficiencies relative to green.
radiometricCorrectionRGB = rawWeights / rawWeights(2);

% Construct a map to apply this correction
radiometricCorrectionMap = ones(size(I));
[bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(radiometricCorrectionMap, bayerPattern);
for cc = 1:3
    radiometricCorrectionMap(bayerIdx{cc}) = radiometricCorrectionRGB(cc);
end

% Report the correction to the console
fprintf('The relative radiometric correction tuple (RGB) is: [%2.4f, %2.4f, %2.4f]\n',radiometricCorrectionRGB);

% Save the radiometric correction to the "derived" directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
readme = ['Created by defineRadiometricWeights.\n'...
    'radiometricCorrectionRGB -- multiply the (linearized) sensor values by these correction factors.\n'...
    'radiometricCorrectionMap -- a map of these corrections that can be applied to an entire image.\n'];
save(saveFileName,'readme','radiometricCorrectionRGB','radiometricCorrectionMap');