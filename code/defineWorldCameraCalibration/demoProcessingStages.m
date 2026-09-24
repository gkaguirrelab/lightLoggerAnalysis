% Housekeeping
clear

% Pick a measurement
measurementName = 'outdoor_AGCandMS_04.mat';
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

% Obtain the demosaiced image
radianceMapDemosaiced = demosaicRadianceMapRCD(imageStages{end});

% Obtain the estimated spectral radiance
[spectralRadiance,miniSpectS] = estimateRadianceSpectrumFromMinispect(minispectValue.AS);

% Show the final, demosaiced image
figure
logImage = log10(radianceMapDemosaiced);
logImage = logImage-min(logImage(:));
logImage = logImage/max(logImage(:));
imagesc(logImage)
box off
axis off
axis equal
title('log10 radiance');
