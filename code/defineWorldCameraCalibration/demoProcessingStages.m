% Housekeeping
clear

% Pick a measurement
%measurementName = 'outdoor_AGCandMS_01.mat';
measurementName = 'planetarium_AGCandMS_01.mat';
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

% Obtain the demosaiced image
radianceMapDemosaiced = demosaicRadianceMap(imageStages{end});

% Generate the sensor -> cone mapping for this observer
fprintf('Generating the sensor -> cone mapping...');
coneMapVar = generateSensorToConeMapping('age',22);
fprintf('done.\n')

% Obtain the cone isomerization map
isomerizationMap = computeConeIsomerizationMap(radianceMapDemosaiced, coneMapVar);

% Obtain the Linear Opponent post receptoral map
postRecepMap = computePostRecepExcitation(isomerizationMap);

% Display the reconstruction
imageStages{end+1}=postRecepMap;
plotReconstructionStages(imageStages)

% Show the final, demosaiced images
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
surf(postRecepMap(:,:,1),'FaceColor','k','FaceAlpha',0.5,'EdgeColor','none');
hold on
surf(postRecepMap(:,:,2),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none');
surf(postRecepMap(:,:,3),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none');
box off
a = gca();
a.XDir ="reverse";
view([-160,20])
zlabel('log pooled isomerization rate [R^*/c/s]');
title('Post-receptoral channels');
