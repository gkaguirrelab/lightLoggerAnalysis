% Dark noise.

% Housekeeping
clear
close all

% Identify the location of the dark noise measurements.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'darkNoise',...
    'AGCstate2',...
    '*.tiff');
fileSet = dir(dataFileName);

% What is the Bayer pattern in these data?
bayerPattern = "BGGR";

% Loop through the images
pixelValsRGB = {[],[],[]};
for ii = 1:length(fileSet)

    % Load the image
    fileName = fullfile(fileSet(ii).folder,fileSet(ii).name);
    I = double(imread(fileName));    

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
% 
% % Take the mean of the pixel values within each channel (R,G,B) across
% % images and pixels within the region of interest
% sensorValues = cellfun(@(x) mean(x),pixelValsRGB);
% 
% % Normalize the sensor values to the blue channel
% sensorValuesNormed = sensorValues / sensorValues(3);
% 
% % Obtain the relative predicted channel values based upon the cloudy sky
% % SPD. We first need to match up the wavelength support of the two
% % measurements. This is hard-coded at the moment, given that the data
% % inputs are known fixed.
% endIdx = find(T.wls == 780);
% assert(wls(1) == T.wls(1));
% assert(wls(end) == T.wls(endIdx));
% assert(diff(wls(1:2)) == diff(T.wls(1:2)));
% 
% % Set up a figure
% figure
% plot(wls,radiance,'.k');
% hold on
% xlabel('wavelength [nm]');
% ylabel('radiance [w/m^2/sr/nm]');
% 
% % Loop through the channels
% channelNames = {'red','green','blue'};
% for ii = 1:3
%     thisChannelSensitivity = T.(channelNames{ii})(1:endIdx);
%     thisChannelSensitivityNormed = thisChannelSensitivity ./ max(thisChannelSensitivity);
%     predictedSensorValue(ii) = radiance' * thisChannelSensitivityNormed;
%     plot(T.wls(1:endIdx),thisChannelSensitivityNormed*predictedSensorValue(ii)/100,['-' channelNames{ii}(1)]);
% end
% 
% % Normalize the predicted sensor values to the blue channel
% predictedSensorValueNormed = predictedSensorValue / predictedSensorValue(3);
% 
% % Calculate the radiometric correction that must be applied to the observed
% % sensor values to have them match the predicted sensor values
% radiometricCorrectionRGB = predictedSensorValueNormed ./  sensorValuesNormed;
% 
% % We further adjust this triplet so that the mean sensor value (across RGB)
% % is unchanged by this operation
% k = 3/sum(radiometricCorrectionRGB);
% radiometricCorrectionRGB = radiometricCorrectionRGB * k;
% 
% % Report the correction to the console
% fprintf('The radiometric correction tuple (RGB) is: [%2.2f, %2.2f, %2.2f]\n',radiometricCorrectionRGB);
% 
% % Save the radiometric correction to the "derived" directory
% saveFileName = fullfile(...
%     tbLocateProjectSilent('lightLoggerAnalysis'),...
%     'derived',...
%     'radiometricCorrectionRGB.mat');
% readme = ['Created by defineRadiometricWeights.\n'...
%     'radiometricCorrectionRGB -- multiple the (linearized) sensor values by these correction factors.\n'];
% save(saveFileName,'readme','radiometricCorrectionRGB');
