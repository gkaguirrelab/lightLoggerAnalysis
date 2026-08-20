% Dark signal.

% Housekeeping
clear
close all

% Define the set of sensor conditions under which the measurements were
% made
stateLabels = {'AGCstate0','AGCstate1','AGCstate2'};

% What is the Bayer pattern in these data?
bayerPattern = "BGGR";

% Loop over the states
for ss = 1:length(stateLabels)

% Identify the location of the dark measurements.
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
    frameSet(:,:,1) = double(imread(fileName));    
end

% Obtain the mean frame
meanFrame{ss} = mean(frameSet,3);

% Obtain the median pixel value
pixelValues(ss) = median(frameSet(:));

end % states


    
