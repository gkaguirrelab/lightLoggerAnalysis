% Housekeeping
clear
%close all

% Pick a measurement
measurementName = 'indoor_AGCandMS_01.mat';
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

% Demosaic the image
imageStages{end+1} = demosaicRadianceMapRCD(imageStages{end});

% Display
figure
plotReconstructionStages(imageStages)

figure
logImage = log10(imageStages{end});
logImage = logImage-min(logImage(:));
logImage = logImage/max(logImage(:));
imagesc(logImage)
title('log10 radiance');
