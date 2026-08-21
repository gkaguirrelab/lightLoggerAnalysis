% This script defines the relationship between camera exposure and the
% median dark signal of the IMX219 camera sensor. For each fixed AGC state,
% the median frame is obtained across the set of dark images, followed by
% the median pixel value of that median frame. A linear model is then fit to
% the median sensor values at the unique exposure settings.

% Housekeeping
clear
close all

% Define the set of sensor conditions under which the measurements were
% made. The labels are listed in natural state-number order.
stateLabels = { ...
    'AGCstate1_AGain-1_DGain-1_E-39', ...
    'AGCstate2_AGain-1_DGain-1_E-4180', ...
    'AGCstate3_AGain-1_DGain-1_E-8290', ...
    'AGCstate4_AGain-5.5652174949645996_DGain-1_E-8290', ...
    'AGCstate5_AGain-10.666666984558105_DGain-1_E-8290'};

% Parse the analog gain, digital gain, and exposure from each state label
AGainValues = zeros(size(stateLabels));
DGainValues = zeros(size(stateLabels));
exposureValues = zeros(size(stateLabels));
for ss = 1:length(stateLabels)
    settingTokens = regexp(stateLabels{ss},...
        'AGain-(?<AGain>[\d.]+)_DGain-(?<DGain>[\d.]+)_E-(?<Exposure>[\d.]+)$',...
        'names');
    AGainValues(ss) = str2double(settingTokens.AGain);
    DGainValues(ss) = str2double(settingTokens.DGain);
    exposureValues(ss) = str2double(settingTokens.Exposure);
end

% Loop over the states
medianFrames = cell(size(stateLabels));
medianPixelValues = zeros(size(stateLabels));
for ss = 1:length(stateLabels)

    % Identify the location of the dark measurements
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'darkNoise',...
        stateLabels{ss},...
        '*.tiff');
    fileSet = dir(dataFileName);

    % Load the images
    frameSet = [];
    for ii = 1:length(fileSet)
        fileName = fullfile(fileSet(ii).folder,fileSet(ii).name);
        frameSet(:,:,ii) = double(imread(fileName));
    end

    % Obtain the median frame
    medianFrames{ss} = median(frameSet,3);

    % Obtain the median pixel value of the median frame
    medianPixelValues(ss) = median(medianFrames{ss}(:));

end % states

% Select the states in which exposure varies while analog and digital gain
% are both fixed at one
exposureStateMask = AGainValues == 1 & DGainValues == 1;
exposureStateIdx = find(exposureStateMask);

% Obtain the unique exposures and their corresponding median pixel values,
% preserving natural state-number order
[uniqueExposureValues,uniqueIdx] = unique(...
    exposureValues(exposureStateMask),'stable');
uniqueExposureIdx = exposureStateIdx(uniqueIdx);
uniqueMedianPixelValues = medianPixelValues(uniqueExposureIdx);

% Fit a linear model to median sensor value as a function of exposure
exposureFitCoefficients = polyfit(...
    uniqueExposureValues,uniqueMedianPixelValues,1);
darkSignalFit = struct(...
    'intercept',exposureFitCoefficients(2),...
    'slope',exposureFitCoefficients(1));

% Evaluate the fit over the measured exposure range
xExposureFit = linspace(...
    min(uniqueExposureValues),max(uniqueExposureValues),100);
yExposureFit = polyval(exposureFitCoefficients,xExposureFit);

% Select the states in which analog gain varies while exposure is fixed at
% its maximum value and digital gain is fixed at one
maxExposure = max(exposureValues);
AGainStateMask = exposureValues == maxExposure & DGainValues == 1;
AGainStateIdx = find(AGainStateMask);

% Obtain the unique analog gains and their corresponding median pixel
% values, again preserving natural state-number order
[uniqueAGainValues,uniqueAGainIdx] = unique(...
    AGainValues(AGainStateMask),'stable');
uniqueAGainIdx = AGainStateIdx(uniqueAGainIdx);
uniqueAGainMedianPixelValues = medianPixelValues(uniqueAGainIdx);

% Fit a linear model to median sensor value as a function of analog gain
AGainFitCoefficients = polyfit(...
    uniqueAGainValues,uniqueAGainMedianPixelValues,1);
darkSignalAGainFit = struct(...
    'intercept',AGainFitCoefficients(2),...
    'slope',AGainFitCoefficients(1));

% Evaluate the fit over the measured analog-gain range
xAGainFit = linspace(min(uniqueAGainValues),max(uniqueAGainValues),100);
yAGainFit = polyval(AGainFitCoefficients,xAGainFit);

% Create a figure
figure('Position',[1 1 1000 400]);
tiledlayout(1,2)

% Plot the measurements and fitted line as a function of exposure
nexttile
plot(uniqueExposureValues,uniqueMedianPixelValues,'xk',...
    'LineStyle','none','MarkerSize',10,'LineWidth',2);
hold on
plot(xExposureFit,yExposureFit,'-r','LineWidth',2);
xlabel('Exposure');
ylabel('Median Sensor Value');
title('Camera Dark Signal as a Function of Exposure');
legend({'Measurements','Linear fit'},'Location','best');
box off

% Plot the measurements and fitted line as a function of analog gain
nexttile
plot(uniqueAGainValues,uniqueAGainMedianPixelValues,'xk',...
    'LineStyle','none','MarkerSize',10,'LineWidth',2);
hold on
plot(xAGainFit,yAGainFit,'-r','LineWidth',2);
xlabel('Analog Gain');
ylabel('Median Sensor Value');
title('Camera Dark Signal as a Function of Analog Gain');
legend({'Measurements','Linear fit'},'Location','best');
box off

% Save the linear fits in the "derived" directory for use during camera
% preprocessing
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'darkSignalFit.mat');
readme = ['Created by defineDarkSignal.\n'...
    'darkSignalFit -- a struct containing the intercept and slope of the linear relation\n'...
    'median sensor value = intercept + slope * exposure.\n'...
    'darkSignalAGainFit -- a struct containing the intercept and slope of the linear relation\n'...
    'median sensor value = intercept + slope * analog gain.\n'];
save(saveFileName,'readme','darkSignalFit','darkSignalAGainFit');
