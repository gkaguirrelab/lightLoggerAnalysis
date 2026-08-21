%% defineFlatFieldingFunction.m

% This is an analysis of camera images collected using the Fels Planetarium
% dome. A full description of the manner of data collection is described in
% the readme.md within this directory. Raw frames were selected from each of
% four camera-rotation periods. This routine loads and linearizes the raw
% frames, averages them, separates out the R, G, and B channels, fits a
% flattened 2D Gaussian function, and then uses this fit to generate and
% save a derived correction function that imposes a flat intensity field.

% Housekeeping
clear; close all; clc

% Load the parameters needed to linearize the raw sensor values
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
load(paramFileName,'clippingExponent');

% Identify the raw frames from the Fels Planetarium recording
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'flatFieldingFunction',...
    'rawFrames',...
    '**',...
    '*.tiff');
fileSet = dir(dataFileName);

% What is the Bayer pattern in these data?
bayerPattern = "BGGR";

% What is the channel order?
channelOrder = {'r','g','b'};

% Call local function to obtain the average across frames for each of the
% channels
avgImageByChannel = makeBayerChannelAverages(...
    fileSet,bayerPattern,clippingExponent);

% Loop over the channels and fit a flattened Gaussian
results = struct();
for ii = 1:numel(avgImageByChannel)
    [pFit, modelFit, residual, X, Y] = fitFlattenedGaussian(avgImageByChannel{ii});
    results(ii).source = avgImageByChannel{ii};
    results(ii).pFit = pFit;
    results(ii).modelFit = modelFit;
    results(ii).residual = residual;
end

% Create separate correction maps for the channels
corrR = max(results(1).modelFit(:)) ./ results(1).modelFit;
corrG = max(results(2).modelFit(:)) ./ results(2).modelFit;
corrB = max(results(3).modelFit(:)) ./ results(3).modelFit;

% Create an integrated correction map
[Hsmall, Wsmall] = size(modelFit);
Hfull = Hsmall * 2;
Wfull = Wsmall * 2;
correctionMap = nan(Hfull, Wfull);

% Must match the Bayer pattern used in fitting: BGGR
correctionMap(1:2:end, 1:2:end) = corrB;
correctionMap(1:2:end, 2:2:end) = corrG;
correctionMap(2:2:end, 1:2:end) = corrG;
correctionMap(2:2:end, 2:2:end) = corrR;

% Plot the correction map
figure
imagesc(correctionMap)
axis image
colorbar
title('Correction for spatial fielding function')
xlabel('X Pixel')
ylabel('Y Pixel')

% Plot the image and model fit through the horizontal center for the three
% channels, normalized to the max of the model fit
figure
hold on
for ii=1:length(results)
    maxVal = max(results(ii).modelFit(:));
    yVals = results(ii).source(round(size(results(ii).source,1)/2),:)/maxVal;
    plot(yVals,['.',channelOrder{ii}]);
    yVals = results(ii).modelFit(round(size(results(ii).source,1)/2),:)/maxVal;
    plot(yVals,['-',channelOrder{ii}],'LineWidth',2);
end
title('Normalized image intensity along horiontal center')
ylabel('image intensity relative to max')
xlabel('horizontal pixel position')

% Save the derived fielding functions
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'flatFieldingFunction.mat');
readme = ['Created by defineFlatFieldingFunction.\n'...
    'correctionMap -- the multiplicative correction at each pixel.\n',...
    'corrR, G, B -- the correction map separated for each channel.\n'];
save(saveFileName,'readme',...
    'correctionMap', ...
    'corrR', 'corrG', 'corrB')




%% LOCAL FUNCTIONS

% Load and linearize the raw frames, and obtain the average across frames
% for each Bayer channel
function avgImageByChannel = makeBayerChannelAverages(...
    fileSet,bayerPattern,clippingExponent)

nFrames = length(fileSet);

for k = 1:nFrames

    % Load and linearize the raw frame
    fileName = fullfile(fileSet(k).folder,fileSet(k).name);
    rawFrame = double(imread(fileName));
    linearFrame = linearizeIMX219SensorCounts(...
        rawFrame,clippingExponent);

    % Separate the Bayer channels
    [R, G, B] = splitBayerFrame(linearFrame, bayerPattern);

    if k == 1
        sum_R = zeros(size(R));
        sum_G = zeros(size(G));
        sum_B = zeros(size(B));
    end

    sum_R = sum_R + R;
    sum_G = sum_G + G;
    sum_B = sum_B + B;

end

avg_R = sum_R ./ nFrames;
avg_G = sum_G ./ nFrames;
avg_B = sum_B ./ nFrames;

% Assemble to return
avgImageByChannel = {avg_R,avg_G,avg_B};

end



% A flattened Gaussian that fits our spatial intensity variation well
function [pFit, hotspot_fit, residual, X, Y] = fitFlattenedGaussian(I)

[H, W] = size(I);
[X, Y] = meshgrid(1:W, 1:H);

baseline0 = min(I(:));
amp0 = max(I(:)) - min(I(:));

[~, maxIdx] = max(I(:));
[y0_guess, x0_guess] = ind2sub(size(I), maxIdx);

sigma0 = min(H,W) / 3;
k0 = 1;

flatGauss2D = @(p, x, y) ...
    p(1) + p(2) .* tanh( ...
    p(7) .* exp( ...
    -((x - p(3)).^2 ./ (2*p(5)^2) + ...
    (y - p(4)).^2 ./ (2*p(6)^2)) ) );

modelFun = @(p, xy) flatGauss2D(p, xy(:,1), xy(:,2));

xy = [X(:), Y(:)];
z = I(:);

p0 = [baseline0, amp0, x0_guess, y0_guess, sigma0, sigma0, k0];

lb = [0, 0, 1, 1, 10, 10, 0.01];
ub = [Inf, Inf, W, H, W, H, 100];

opts = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'MaxFunctionEvaluations', 5000);

pFit = lsqcurvefit(modelFun, p0, xy, z, lb, ub, opts);

hotspot_fit = reshape(modelFun(pFit, xy), H, W);
residual = I - hotspot_fit;

end
