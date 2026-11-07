% Housekeeping
clear

% Pick a measurement
measurementName = 'planetarium_AGCandMS_01.mat';
%measurementName = 'outdoor_AGCandMS_02.mat';
%measurementName = 'indoor_AGCandMS_01.mat';
%measurementName = 'macbeth_AGCandMS_02.mat';

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

% Obtain the minispect scene illuminance. Load the CIE 1931 XYZ color
% matching functions. Spline the CMFs to the wavelength support of our
% radiance spectrum Multiplying by 683 scales the functions to the maximum
% luminous efficacy (lm/W).
load('T_xyz1931.mat', 'T_xyz1931', 'S_xyz1931');
T_xyz = SplineCmf(S_xyz1931, 683 * T_xyz1931, S);
minispectIlluminanceLux = T_xyz(2, :) * meanSceneRadiance * pi;

% Obtain the foveated integrated cone radiance, further converted into
% post-receptoral channels.
gazeX = 320;
gazeY = 240;

[opponentConeRadianceMap,coneWeights, AzGrid,ElGrid, cameraIlluminanceLux] = computeFoveatedOpponentMap(radianceMapDemosaiced, ...
    'gazeX', gazeX, ...
    'gazeY', gazeY,...
    'balanceNeutral', true);

% Report the illuminance values
fprintf('minispect illuminance: %2.1f lux; camera illuminance: %2.1f lux\n',minispectIlluminanceLux,cameraIlluminanceLux);

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
