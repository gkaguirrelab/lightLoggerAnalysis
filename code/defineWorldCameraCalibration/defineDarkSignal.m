% This script defines the relationship between camera exposure and the
% median dark signal of the IMX219 camera sensor. Video data were collected
% while the camera was in a state of total darkness (lens cap on, wrapped
% in black cloth, in a dark room). These data were collected for a range of
% AGC settings to explore the relationship between the dark signal value
% and camera sensitivity (principally exposure time). As we find that the
% dark signal is the same across all recording states, we simply save the
% single dark single value into the derived parameters directory.

% Housekeeping
clear

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
        'darkSignal',...
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

% Test the assertion that the medianPixelValues contain a single unique
% value across camera AGC settings
darkSignal = unique(medianPixelValues);
assert(isscalar(darkSignal));

% Save the darkSignal in the "derived" directory for use during camera
% preprocessing
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'darkSignal.mat');
README = ['Created by defineDarkSignal.\n'...
    'darkSignal -- the scalar dark signal value observed in the camera data.\n'];
derivedVariables = struct();
derivedVariables.darkSignal = darkSignal;
saveDerivedFile(saveFileName, README, derivedVariables);
