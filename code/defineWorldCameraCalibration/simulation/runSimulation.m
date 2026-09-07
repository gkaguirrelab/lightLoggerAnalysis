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
visualizeRadianceModel(radianceModel,radianceModelS);

% Obtain the camera radiance map
cameraRadianceMap = projectRadianceModelToCameraRadiance(radianceModel,radianceModelS);

% Find the camera score that would be predicted for this radiance map. This
% non-linear search is necessary as the pixel saturation as a function of
% the camera settings has a non-linear interaction with the calculated
% camera score.
myMap = @(x) synthesisPipeline(cameraRadianceMap, cameraScoreToAGCSettings(x));
myObj = @(x) norm(mean(mean(myMap(x)))-127);
cameraScore = fminsearch(myObj,8000);
AGCSettings = cameraScoreToAGCSettings(cameraScore);

% Synthesize the raw camera image for this radianceModel and camera score
I = synthesisPipeline(cameraRadianceMap, AGCSettings);

% Reconstruct the radiance map from the image and AGCSettings
[~,imageStages] = reconstructionPipeline(I,AGCSettings);

% Derive the minispect for this radianceModel
[minispectValues, meanSpectralRadiance] = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS);
minispectData.AS = minispectValues;

% Impute the missing values in the cameraRadianceMapEstimated
imageStages{end+1} = imputePixelValues(imageStages{end},minispectData);
cameraRadianceMapEstimated = imageStages{end};

% Save the source radiance map in the end of the imageStages
imageStages{end+1} = cameraRadianceMap;

% Obtain the estimate of the mean spectral radiance of the radianceModel
% from the minispect values so we can compare how well all this worked.
[spectralRadianceFromMinispect,spectralRadianceFromMinispectS] = estimateRadianceSpectrumFromMinispect(minispectValues);

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
tiledlayout(1,3)
nexttile
imagesc(cameraRadianceMap); colorbar; title('Source radiance map');
nexttile
imagesc(cameraRadianceMapEstimated); colorbar; title('Reconstructed radiance map');
nexttile
imagesc(cameraRadianceMap-cameraRadianceMapEstimated); colorbar; title('error map');

% Show the reconstruction stages
plotReconstructionStages(imageStages)
