function I = imputeSaturatedPixelValues(I,minispectData)
% This function imputes absolute values for saturated pixels

persistent deltaSteradians
if isempty(deltaSteradians)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'deltaSteradians.mat');
    load(paramFileName,'deltaSteradians');
end

persistent cameraT cameraWls
if isempty(cameraT)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(paramFileName,'T');
    cameraT = table2array(T(:,["red" "green" "blue"]))';
    cameraWls = T.wls;
end

% Derive an estimate of the environmental SPD from the minispect
[miniSpectSPD,miniSpectS] = estimateRadianceSpectrumFromMiniSpect(minispectData.AS(1:9));
miniSpectSPD = SplineRaw(SToWls(miniSpectS),miniSpectSPD,cameraWls);
miniSpectS = WlsToS(cameraWls);

% Obtain an estimate of the mean scene radiance experienced by the camera
miniSpectAvgSceneRadiance = mean(miniSpectSPD'*(cameraT ./ max(cameraT,[],2))');

% Obtain the visual angle area (in steradians) of the saturated pixels, and
% the total visual surface area of the camrea
infIdx = isinf(I);
cameraSaturatedSteradians = sum(deltaSteradians(infIdx));
cameraTotalSteradians = sum(deltaSteradians(:));

% Obtain the total radiance as seen by the unsaturated pixels, weighted by
% the visual angle area of each of these pixels
unsatPartition = sum(I(~infIdx) .* deltaSteradians(~infIdx));

% Calculate the mean radiance that should be present in the saturated pixels
cameraSaturatedAvgRadiance = ((miniSpectAvgSceneRadiance * cameraTotalSteradians) - unsatPartition) / cameraSaturatedSteradians;

% Assign the imputed saturated radiance to the inf values in the image
I(find(infIdx)) = cameraSaturatedAvgRadiance;

end
