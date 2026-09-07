% Simulation of radiance map reconstruction. We begin by creating a model
% of a continuous light field with variation in spectral content and
% radiance. This radianceModel is then used to generate a
% cameraRadianceMap, the expected AGCSettings, and the minispect values. We
% then synthesize the camera image that would be expected. Finally, we
% attempt to reconstruct the radiance map given this camera image.

clear
close all

% Create the radiance model
[radianceModel,radianceModelS] = generateSimulatedLightField();
%[radianceModel,radianceModelS] = generateUniformLightField();
visualizeRadianceModel(radianceModel,radianceModelS);

% Derive the minispect values for this radianceModel. We also save the mean
% spectral radiance of the 
[minispectValues, meanSpectralRadiance] = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS);
minispectData.AS = minispectValues;

% Obtain the estimate of the mean spectral radiance of the radianceModel
% from the minispect values so we can compare how well all this worked.
[spectralRadianceFromMinispect,spectralRadianceFromMinispectS] = estimateRadianceSpectrumFromMinispect(minispectValues);

% Find the camera score that would be predicted for this radianceModel.
% This non-linear search is necessary as the pixel saturation as a function
% of the camera settings has a non-linear interaction with the calculated
% camera score.
myMap = @(x) synthesisPipeline(radianceModel,radianceModelS,cameraScoreToAGCSettings(x));
myObj = @(x) norm(mean(mean(myMap(x)))-127);
cameraScore = fminsearch(myObj,8000);
AGCSettings = cameraScoreToAGCSettings(cameraScore);

% Synthesize the raw camera image for this radianceModel and camera score
[I, cameraRadianceMap] = synthesisPipeline(radianceModel,radianceModelS, AGCSettings);

% Report the channel energies for the cameraRadianceMap directly, and as
% seen by the minispect
channelEnergy = calculateIMX219ChannelRadiantEnergyDirect(cameraRadianceMap);
channelEnergy = calculateIMX219ChannelRadiantEnergyViaMS(radianceModel, radianceModelS);

% Reconstruct the radiance map from the image and AGCSettings
[~,imageStages] = reconstructionPipeline(I,AGCSettings);


% Impute the missing values in the cameraRadianceMapEstimated
imageStages{end+1} = imputePixelValues(imageStages{end},minispectData);
cameraRadianceMapEstimated = imageStages{end};

% Save the source radiance map in the end of the imageStages
imageStages{end+1} = cameraRadianceMap;

% Visualize how well we did reconstruction the mean spectral radiance of
% the light field
figure
plot(SToWls(radianceModelS),meanSpectralRadiance,'.k');
hold on
plot(SToWls(spectralRadianceFromMinispectS),spectralRadianceFromMinispect,'-r');
xlim([380 780]);
xlabel('wavelength [nm]');
ylabel('Radiance [W/m2/sr/nm]');
legend({'source','reconstruction'});
title('Reconstruction of hemifield mean spectral radiance')

% Visualize our ability to reconstruct the camera radiance map
figure
tiledlayout(1,4)
nexttile
surf(cameraRadianceMap, 'EdgeColor', 'none'); colorbar; title('Source'); zlim([0 30]);
nexttile
surf(imageStages{end-2}, 'EdgeColor', 'none'); colorbar; title('Reconstructed'); zlim([0 30]);
nexttile
surf(cameraRadianceMapEstimated, 'EdgeColor', 'none'); colorbar; title('Imputed'); zlim([0 30]);
nexttile
surf(cameraRadianceMapEstimated-cameraRadianceMap, 'EdgeColor', 'none'); colorbar; title('Error');

% Show the reconstruction stages
plotReconstructionStages(imageStages)
