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

% Demosaic the image
imageStages{end+1} = demosaicRadianceMapRCD(imageStages{end});

% Display
figure
plotReconstructionStages(imageStages)

figure
imagesc(log10(imageStages{end}))
title('log10 radiance');
