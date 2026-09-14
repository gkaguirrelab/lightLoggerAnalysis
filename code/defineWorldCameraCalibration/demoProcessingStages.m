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

% Bilinear demosaic
imageStages{end+1} = demosaicRadianceMap(imageStages{end});

% Impute radiance values for saturated areas
imageStages{end+1} = imputePixelValueBayes(imageStages{end});

plotReconstructionStages(imageStages)

