% Housekeeping
clear
%close all

% Pick a measurement
measurementName = 'outdoor_AGCandMS_01.mat';
%measurementName = 'planetarium_AGCandMS_01.mat';

% Load the data
fileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'exampleWorldCameraImages',...
    measurementName);
load(fileName,'worldFrame','AGCSettings','minispectValue')

% Reconstruct the radiance image
[radianceMap,imageStages] = reconstructionPipeline(worldFrame,AGCSettings);

% Display the reconstruction
plotReconstructionStages(imageStages)

% Obtain the demosaiced image and show this
radianceMapDemosaiced = demosaicRadianceMapRCD(imageStages{end});

figure
logImage = log10(radianceMapDemosaiced);
logImage = logImage-min(logImage(:));
logImage = logImage/max(logImage(:));
imagesc(logImage)
box off
axis off
axis equal
title('log10 radiance');
