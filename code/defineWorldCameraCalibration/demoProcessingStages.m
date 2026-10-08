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

radianceMapDemosaiced(:,:,1) = 1.0;
radianceMapDemosaiced(:,:,2) = 1.1;
radianceMapDemosaiced(:,:,3) = 1.2;


% Obtain the cone isomerization map, which is in units of R*/deg^2/s. We
% pass the current fixation location (in terms of image pixel) and pupil
% size to be used in this calculation
observerAge = 25;
gazeX  = 320;
gazeY = 240;
pupilDiameterMm = 3.0;
[isomerizationMap,AzGrid,ElGrid] = computeFoveatedConeIsomerizationMap(radianceMapDemosaiced, ...
    'observerAge',observerAge,...
    'gazeX', gazeX, ...
    'gazeY', gazeY, ...
    'pupilDiameterMm', pupilDiameterMm);

% Obtain the Linear Opponent post-receptoral map
[postRecepMap,weights] = computePostRecepExcitation(isomerizationMap,'balanceNeutral',false);

% Display the reconstruction
plotReconstructionStages(imageStages)

% Show postRecepMap channels
faceColors = {'k','r','b'};
figure
for cc = 1:3
    surf(AzGrid,ElGrid,postRecepMap(:,:,cc),'FaceColor',faceColors{cc},'FaceAlpha',0.5);
    hold on
    zlabel('pooled opponent isomerization rate [R^*/deg^2/s]');
end
xlim([-100 100]);
ylim([-50 50]);
pbaspect([2 1 1])
title('Post-receptoral channels');
