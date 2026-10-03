% Housekeeping
clear

% Pick a measurement
measurementName = 'outdoor_AGCandMS_01.mat';
%measurementName = 'planetarium_AGCandMS_01.mat';
%measurementName = 'macbeth_AGCandMS_01.mat';

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
radianceMapDemosaiced = demosaicRadianceMap(imageStages{end});

% Generate the sensor -> cone mapping for this observer
fprintf('Generating the sensor -> cone mapping...');
coneMapVar = generateSensorToConeMapping('age',22);
fprintf('done.\n')

% Obtain the cone isomerization map
isomerizationMap = computeConeIsomerizationMap(radianceMapDemosaiced, coneMapVar);

% Show the final, demosaiced image
figure
tiledlayout(1,2,'Padding','tight','TileSpacing','compact')
nexttile
logImage = log10(radianceMapDemosaiced);
logImage = logImage-min(logImage(:));
logImage = logImage/max(logImage(:));
imagesc(logImage)
box off
axis off
axis equal
title('log10 radiance');

nexttile
logImage = log10(isomerizationMap);
logImage = logImage-min(logImage(:));
logImage = logImage/max(logImage(:));
imagesc(logImage)
box off
axis off
axis equal
title('log10 isomerization rate');
