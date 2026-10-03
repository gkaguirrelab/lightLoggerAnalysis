% validateCameraBySphereImages.m
%
% This script validates the IMX219 camera radiometric reconstruction
% pipeline by passing the flat calibration images (obtained in the
% integrating sphere) back through the pipeline. It compares the
% reconstructed mean integrated radiance and the RGB channel fractions to
% the ground-truth predictions derived from the simultaneous PR670 spectral
% measurements.

% Housekeeping
clear
close all

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing)
commonS = [380, 1, 352];
commonWls = SToWls(commonS);

% Load the IMX219 sensitivity functions
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');
wlsSensor = T.wls;
channelNames = {'red','green','blue'};

% Spline the sensor sensitivities to the common 1 nm wavelength domain
% and scale so maximum value is unity
sensorSensitivities = zeros(length(commonWls), 3);
for cc = 1:length(channelNames)
    sens = SplineRaw(wlsSensor, T.(channelNames{cc}), commonWls);
    sensorSensitivities(:, cc) = sens ./ max(sens);
end

% Define the ND levels used in calibration
ndfLevels = [0 1 2 3];

% Preallocate arrays to hold predicted and measured radiance values
predictedRadiance = zeros(length(ndfLevels), 3);
measuredRadiance = zeros(length(ndfLevels), 3);
radianceMaps = cell(length(ndfLevels), 1);

% Loop over each ND level to process PR670 and camera data
for ii = 1:length(ndfLevels)
    
    fprintf('Processing NDF %d...\n', ndfLevels(ii));

    % ----------------------------------------------------
    % 1. PREDICTED RADIANCE (via PR670)
    % ----------------------------------------------------
    pr670FileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'PR670',...
        sprintf('AGCSettingsMeasure%dNDF.mat', ndfLevels(ii)));
    load(pr670FileName, 'measurement', 'S');
    
    % Average radiance and spline to 1 nm spacing
    spdSource_raw = mean(measurement, 1);
    spdSource = SplineSpd(SToWls(S), spdSource_raw', commonWls)';
    
    % Calculate the predicted integrated radiance using the dot product
    for cc = 1:3
        predictedRadiance(ii, cc) = spdSource * sensorSensitivities(:, cc);
    end
    
    % ----------------------------------------------------
    % 2. MEASURED RADIANCE (via Camera Pipeline)
    % ----------------------------------------------------
    cameraFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'worldCamera',...
        sprintf('NDF%d', ndfLevels(ii)),...
        sprintf('NDF%d_AGCandMS_01.mat', ndfLevels(ii)));
    load(cameraFileName, 'AGCSettings', 'worldFrame');
    
    % Convert the raw camera frame to a radiance map using the pipeline
    [radianceMap, ~] = reconstructionPipeline(worldFrame, AGCSettings);
    radianceMaps{ii} = radianceMap;
    
    % Demosaic the image to obtain RGB triplets for each pixel
    radianceMapRGB = demosaicRadianceMap(radianceMap);
    
    % Calculate the mean radiance for the R, G, and B channels across the image
    measuredRadiance(ii, 1) = mean(radianceMapRGB(:, :, 1), 'all', 'omitnan');
    measuredRadiance(ii, 2) = mean(radianceMapRGB(:, :, 2), 'all', 'omitnan');
    measuredRadiance(ii, 3) = mean(radianceMapRGB(:, :, 3), 'all', 'omitnan');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------
% FIGURE 1: Radiance Images by NDF Level
% ----------------------------------------------------
figure('Name', 'Radiance Images by ND Level', 'WindowStyle', 'Docked', ...
       'Position', [50, 500, 1200, 300]);
tiledlayout(1, length(ndfLevels), "TileSpacing", "compact", "Padding", "tight");

for ii = 1:length(ndfLevels)
    nexttile;
    logImage = log10(radianceMaps{ii});
    logImage = logImage - min(logImage(:));
    logImage = logImage / max(logImage(:));
    imagesc(logImage);
    colormap gray;
    box off;
    axis off;
    axis equal;
    title(sprintf('NDF %d', ndfLevels(ii)));
end

% ----------------------------------------------------
% FIGURE 2: Calibration Validation
% ----------------------------------------------------
figure('Name', 'Calibration Validation', 'WindowStyle', 'Docked', ...
       'Position', [50, 50, 1000, 400]);
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight");

% Calculate the mean across the R, G, and B channels
predMean = mean(predictedRadiance, 2);
measMean = mean(measuredRadiance, 2);

% Define shades of gray for each ND level (lightest to darkest)
grayColors = [0.9 0.9 0.9; 0.7 0.7 0.7; 0.4 0.4 0.4; 0.1 0.1 0.1];

% --- TILE 1: Mean Integrated Radiance ---
nexttile; hold on;

% Plot each ND level point and capture the handle for the legend
hScat = zeros(length(ndfLevels), 1);
for ii = 1:length(ndfLevels)
    hScat(ii) = scatter(predMean(ii), measMean(ii), 120, grayColors(ii,:), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
end

% Determine axis limits dynamically
maxVal = max([predMean(:); measMean(:)]) * 1.2;
minVal = min([predMean(:); measMean(:)]) * 0.8;

% 1:1 identity line
plot([minVal, maxVal], [minVal, maxVal], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

axis square;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([minVal, maxVal]);
ylim([minVal, maxVal]);
xlabel('Predicted via PR670 [W/m2/sr]');
ylabel('Measured via IMX219 Pipeline [W/m2/sr]');
title('Mean Integrated Radiance');
grid on; box on;

% Create a legend using the NDF levels
legendStr = arrayfun(@(x) sprintf('NDF %d', x), ndfLevels, 'UniformOutput', false);
legend(hScat, legendStr, 'Location', 'northwest');
hold off;

% --- TILE 2: R:G:B Channel Fractions ---
nexttile; hold on;

% Calculate total radiance per NDF level to extract fractional channel ratios
predSum = sum(predictedRadiance, 2);
measSum = sum(measuredRadiance, 2);

predFracR = predictedRadiance(:, 1) ./ predSum;
predFracG = predictedRadiance(:, 2) ./ predSum;
predFracB = predictedRadiance(:, 3) ./ predSum;

measFracR = measuredRadiance(:, 1) ./ measSum;
measFracG = measuredRadiance(:, 2) ./ measSum;
measFracB = measuredRadiance(:, 3) ./ measSum;

% Plot each fractional channel color contribution
scatter(predFracR, measFracR, 85, 'r', 'filled', 'MarkerEdgeColor', 'k');
scatter(predFracG, measFracG, 85, 'g', 'filled', 'MarkerEdgeColor', 'k');
scatter(predFracB, measFracB, 85, 'b', 'filled', 'MarkerEdgeColor', 'k');

% 1:1 identity line
plot([0, 0.75], [0, 0.75], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

axis square;
xlim([0, 0.75]);
ylim([0, 0.75]);
set(gca, 'XTick', [0 0.25 0.5 0.75], 'YTick', [0 0.25 0.5 0.75]);
xlabel('Predicted Channel Fraction');
ylabel('Measured Channel Fraction');
title('Relative Channel Fractions (R:G:B Ratio)');
grid on; box on;
hold off;

fprintf('\nValidation complete.\n');