% Housekeeping
clear

% Pick a measurement
%measurementName = 'planetarium_AGCandMS_01.mat';
%measurementName = 'outdoor_AGCandMS_02.mat';
%measurementName = 'indoor_AGCandMS_01.mat';
measurementName = 'macbeth_AGCandMS_02.mat';

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

% Obtain the mean environmental radiance from the  minispect and
% radianceMap
[meanSceneRadiance, S, fVal] = ...
    estimateJointRadianceSpectrum(minispectValue.AS, radianceMapDemosaiced);

% Obtain the foveated integrated cone radiance, further converted into
% post-receptoral channels.
gazeX = 320;
gazeY = 240;

[opponentConeRadianceMap,coneWeights, AzGrid,ElGrid] = computeFoveatedOpponentMap(radianceMapDemosaiced, ...
    'gazeX', gazeX, ...
    'gazeY', gazeY,...
    'balanceNeutral', true);

% Display the reconstruction
plotReconstructionStages(imageStages)

% Show postRecepMap channels
faceColors = {'k','r','b'};
figure
for cc = 1:3
    surf(AzGrid,ElGrid,opponentConeRadianceMap(:,:,cc),'FaceColor',faceColors{cc},'FaceAlpha',0.5);
    hold on
    zlabel('opponent cone radiance');
end
xlim([-100 100]);
ylim([-50 50]);
pbaspect([2 1 1])
title('Post-receptoral channels');

% Plot the mean scene radiance
figure
plot(SToWls(S),meanSceneRadiance,'r');
xlabel('wavelength [nm]');
ylabel('radiance [W/m^2/sr]');
title('Mean scene radiance')
