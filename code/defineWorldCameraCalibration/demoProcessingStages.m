
% Housekeeping
clear

% Identify a raw TIFF
tiffFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'exampleWorldCameraImages',...
    'planetarium_1.tiff');

% Define a saturation threshold value
satThresh = 250;

% Load the image and convert to double
imageStages{1} = double(imread(tiffFileName));

% Linearize
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
load(paramFileName,'clippingExponent');
imageStages{2} = linearizeIMX219SensorCounts(...
    imageStages{1},clippingExponent);

% Threshold
I = imageStages{2};
I(I>satThresh) = Inf;
imageStages{3} = I;

% Flat fielding function
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'flatFieldingFunction.mat');
load(paramFileName,'correctionMap');
imageStages{4} = imageStages{3} .* correctionMap;

% Radiometric correction
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
load(paramFileName,'radiometricCorrectionRGB');
bayerPattern = "BGGR";
I = [];
[I(:,:,1),I(:,:,2),I(:,:,3)] = splitBayerFrame(imageStages{4}, bayerPattern);
I = I .* reshape(radiometricCorrectionRGB, 1, 1, 3);
imageStages{5} = combineBayerFrame(I(:,:,1),I(:,:,2),I(:,:,3), bayerPattern);

% Plot

stages = {'raw','threshold','linearized','flattened','radiometric correction'};
channelColor = {'r','g','b'};

figure
tiledlayout(2,5,"TileSpacing","tight");
for ss = 1:5
    nexttile(ss)

    % Get the image and some basic stats
    I = imageStages{ss};
    isinfI = isinf(I);
    nInf = sum(isinfI(:));
    maxVal = max(I(~isinfI));
    showI = I;
    showI(isinfI) = maxVal;

    % Show the image
    imshow(I,[]);
    title(stages{ss});

    % Show a histogram
    nexttile(ss+5)
    I(isinfI) = nan;
    [rgb(:,:,1),rgb(:,:,2),rgb(:,:,3)] = splitBayerFrame(I, bayerPattern);
    edges = 0:0.01:1;
    for cc = 1:3

        % Grab this channel
        vec = rgb(:,:,cc);

        % Record min and max prior to normalizing, rounded to nearest whole number
        minVals(cc) = round(min(vec(:)));
        maxVals(cc) = round(max(vec(:)));

        vec = vec(:)/maxVal;
        N = histcounts(vec,edges);
        plot(edges(1:end-1),N,['-' channelColor{cc}]);
        hold on
    end

    % Add min and max text to the upper left corner
    textStr = {sprintf('min = [%d, %d, %d]', minVals(1), minVals(2), minVals(3)), ...
        sprintf('max = [%d, %d, %d]', maxVals(1), maxVals(2), maxVals(3)),...
        sprintf('sat = %d', nInf)};
    text(0.05, 0.95, textStr, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);

    axis off
    box off
end
