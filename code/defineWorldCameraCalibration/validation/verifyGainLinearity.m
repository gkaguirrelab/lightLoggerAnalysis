% verifyGainLinearity.m
% This script tests whether the camera's pre-gain image values scale 
% inversely and linearly with Digital Gain (Dgain) by comparing Dataset 2 
% (Dgain ~2.00) and Dataset 5 (Dgain ~5.73), ignoring saturated pixels.

clear
close all

% Validation set definitions matching validateCameraByColorChecker
valSetOptions = {'indoor','indoor','indoor','outdoor','outdoor','outdoor'};
valNumOptions = {'1','2','4','1','2','3'};
plotOrder = [4,5,1,2,6,3];

% Target datasets based on DgainStore output:
% vv = 2 -> Dgain ~2.00
% vv = 5 -> Dgain ~5.73
targetIndices = [2, 5]; 
datasetLabels = {'Dgain ~2.00 (Dataset 2)', 'Dgain ~5.73 (Dataset 5)'};

extractedMeans = zeros(1, 2);
dGainValues = zeros(1, 2);

for ii = 1:2
    vv = targetIndices(ii);
    valIdx = plotOrder(vv);

    % Load the raw camera frame and AGC settings
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'macbethColorCheck', ...
        valSetOptions{valIdx}, valNumOptions{valIdx}, 'lightLogger', 'close_AGCandMS_01.mat');
    load(dataFileName, 'worldFrame', 'AGCSettings');

    dGainValues(ii) = AGCSettings.Dgain;

    % Run through reconstructionPipeline to access intermediate stages
    [~, imageStages] = reconstructionPipeline(worldFrame, AGCSettings);

    % Stage 2 contains linearized values, with saturated pixels set to Inf.
    % Compute mean of finite, non-NaN values to avoid Inf contamination.
    validPixels = imageStages{2}(isfinite(imageStages{2}));
    extractedMeans(ii) = mean(validPixels);
    
    fprintf('%s -> Dgain: %.4f, Stage 2 Finite Mean: %.4f\n', ...
        datasetLabels{ii}, dGainValues(ii), extractedMeans(ii));
end

% Calculate theoretical vs empirical ratios
theoreticalRatio = dGainValues(2) / dGainValues(1);
empiricalRatio = extractedMeans(1) / extractedMeans(2);

fprintf('\n--- Gain Linearity Verification Results ---\n');
fprintf('Theoretical Dgain Ratio (Dgain_5 / Dgain_2): %.4f\n', theoreticalRatio);
fprintf('Empirical Pre-Gain Signal Ratio (Mean_2 / Mean_5): %.4f\n', empiricalRatio);
fprintf('Linearity Discrepancy Factor: %.4f\n', empiricalRatio / theoreticalRatio);
fprintf('-------------------------------------------\n');