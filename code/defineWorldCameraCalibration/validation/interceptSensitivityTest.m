% interceptSensitivityTest.m
% This script tests whether the calibration anchor (intercept) is stable 
% across all light levels or if it is being systematically dragged down 
% by non-linearities in the dimmest (high Dgain) integrating sphere scenes.

% Housekeeping
clear
close all

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing)
commonS = [380, 1, 352];
commonWls = SToWls(commonS);

% Load the AGC settings for each ND level
agcData.ndf = [0 1 3];
for ii = 1:length(agcData.ndf)
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'worldCamera',...
        sprintf('NDF%d',agcData.ndf(ii)),...
        sprintf('NDF%d_AGCandMS_01.mat',agcData.ndf(ii)));
    load(dataFileName,'AGCSettings','worldFrame');
    agcData.AGain(ii) = AGCSettings.Again;
    agcData.DGain(ii) = AGCSettings.Dgain;
    agcData.Exposure(ii) = AGCSettings.exposure;
end

% Load parameters required for set point calculation
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'darkSignal.mat'), 'darkSignal');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'nonLinearClippingExponent.mat'), 'clippingExponent');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'flatFieldingFunction.mat'), 'correctionMap');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'radiometricCorrectionRGB.mat'), 'radiometricCorrectionMap');

Smax = 2^8 - 1 - darkSignal;
n = clippingExponent;
meanCorrectionFielding = 1 / mean(1 ./ correctionMap(:), 'omitnan');
meanCorrectionRGB = 1 / mean(1 ./ radiometricCorrectionMap(:), 'omitnan');

% Derive an "effective camera score" that accounts for dark signal
cameraScore = zeros(1, length(agcData.ndf));
for ii = 1:length(agcData.ndf)
    setPoint = 127;
    setPoint = (setPoint / agcData.DGain(ii)) - darkSignal;
    linSetPoint = setPoint / (1 - (setPoint / Smax)^n)^(1/n);
    linSetPoint = linSetPoint * meanCorrectionFielding * meanCorrectionRGB;

    % The effective sensitivity score is the hardware gain divided by the targeted linear signal
    cameraScore(ii) = (agcData.AGain(ii) * agcData.Exposure(ii)) / linSetPoint;
end

% Load the IMX219 sensitivity functions.
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

% Preallocate an Nx3 array for the sensor-weighted integrated radiance
integratedRadiance = zeros(length(agcData.ndf), 3);

% Next, load radiance spectrum associated with each ND level
for ii = 1:length(agcData.ndf)
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'PR670',...
        sprintf('AGCSettingsMeasure%dNDF.mat',agcData.ndf(ii)));
    load(dataFileName,'measurement','S');
    
    % This is the average radiance in units of Watts/m2/sr/[S(2)*nm]
    spdSource_raw = mean(measurement,1);
    
    % Resample the source SPD to 1 nm spacing.
    spdSource = SplineSpd(SToWls(S), spdSource_raw', commonWls)';

    % Loop over the channels to calculate the sensor-weighted effective radiance.
    for cc = 1:length(channelNames)
        integratedRadiance(ii,cc) = spdSource * sensorSensitivities(:, cc);
    end
end

% Get the mean integrated radiance as a function of camera score
bayerWeightedRadiance = (integratedRadiance(:,1) + 2*integratedRadiance(:,2) + integratedRadiance(:,3)) / 4;

% Force the theoretical physical slope of -1
fixedSlope = -1.0;

% =========================================================================
% ANALYTICALLY CALCULATE THE INDIVIDUAL INTERCEPTS
% =========================================================================
% Instead of taking the mean(), we retain the individual intercept for each NDF
individualIntercepts = log10(bayerWeightedRadiance(:)) - (fixedSlope * log10(cameraScore(:)));
digitalGains = agcData.DGain(:);

% =========================================================================
% PLOT RESULTS
% =========================================================================
figure('Name', 'Intercept Sensitivity Test', 'WindowStyle', 'Docked', ...
       'Position', [100, 100, 600, 500]);
    
scatter(digitalGains, individualIntercepts, 100, 'filled', 'MarkerEdgeColor', 'k');
hold on;

% Add text labels for NDF level next to each point
for ii = 1:length(digitalGains)
    text(digitalGains(ii) + 0.1, individualIntercepts(ii), sprintf('NDF %d', agcData.ndf(ii)), ...
        'FontSize', 10, 'VerticalAlignment', 'bottom');
end

% Add a horizontal line representing the mean intercept (what your pipeline currently uses)
meanIntercept = mean(individualIntercepts);
yline(meanIntercept, 'r--', 'Mean Intercept (Current Calibration Anchor)', ...
    'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

% Formatting
xlabel('Digital Gain (Dgain)');
ylabel('Calculated Log_{10} Intercept');
title('Calibration Intercept Stability vs. Digital Gain');
grid on; box on;

% Expand X-axis limits to accommodate text labels
xlim([0, max(digitalGains) * 1.2]);

% Print to console
fprintf('\n--- Intercept Sensitivity Results ---\n');
for ii = 1:length(digitalGains)
    fprintf('NDF %d (Dgain %.2f): %.4f\n', agcData.ndf(ii), digitalGains(ii), individualIntercepts(ii));
end
fprintf('-------------------------------------\n');
fprintf('Mean Intercept (Global): %.4f\n\n', meanIntercept);