
% Housekeeping
clear

% Pick a measurement
measurementName = 'planetarium_AGCandMS_01.mat';

% Identify a raw TIFF
fileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'exampleWorldCameraImages',...
    measurementName);

% Load the image and convert to double
load(fileName,'worldFrame','AGCSettings','minispectValue')
imageStages{1} = double(worldFrame);

% Linearize
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
load(paramFileName,'clippingExponent');
imageStages{2} = linearizeIMX219SensorCounts(...
    imageStages{1},clippingExponent);

% Flat fielding function
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'flatFieldingFunction.mat');
load(paramFileName,'correctionMap');
imageStages{3} = imageStages{2} .* correctionMap;

% Radiometric correction
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
load(paramFileName,'radiometricCorrectionMap');
imageStages{4} = imageStages{3} .* radiometricCorrectionMap;

% Convert to radiance units
imageStages{5} = convertSensorValuesToRadiance(imageStages{4},AGCSettings);

% Impute radiance values for saturated areas
imageStages{6} = imputeSaturatedPixelValues(imageStages{5},minispectValue);

% Plot
stages = {'raw','linearized','flattened','radiometric correction','absolute radiance','impute saturated'};
nStages = length(stages);
channelColor = {'r','g','b'};

figure
tiledlayout(2,nStages,"TileSpacing","tight");
for ss = 1:nStages

    nexttile(ss)

    % Get the image and some basic stats
    I = imageStages{ss};
    isinfI = isinf(I);
    nInf = sum(isinfI(:));
    maxVal = max(I(~isinfI));
    showI = I;
    showI(isinfI) = maxVal;

    % Show the image
    imshow(I,[0 prctile(I(~isinf(I)), 85)]);
    title([sprintf('%d. ',ss),stages{ss}]);

    % Show a histogram
    nexttile(ss+nStages)
    I(isinfI) = nan;
    [rgbIdx{1},rgbIdx{2},rgbIdx{3}] = returnBayerIndices(I, 'BGGR');
    edges = (0:255)/255;
    for cc = 1:3

        % Grab this channel
        vec = I(rgbIdx{cc});

        % Record min and max prior to normalizing, rounded to nearest whole number
        minVals(cc) = round(min(vec(:)));
        maxVals(cc) = round(max(vec(:)));

        % Set the nan values to the maxVal
        vec(isnan(vec))=maxVal;

        vec = vec(:)/maxVal;
        N = histcounts(vec,edges);

        N=N./length(rgbIdx{cc});

        plot(edges(1:end-1)*100,N*100,['.-' channelColor{cc}]);
        hold on
    end

    % Add min and max text to the upper left corner
    textStr = {sprintf('min = [%d, %d, %d]', minVals(1), minVals(2), minVals(3)), ...
        sprintf('max = [%d, %d, %d]', maxVals(1), maxVals(2), maxVals(3)),...
        sprintf('ceil, floor = [%d, %d]', nInf,sum(I(:)==0))};
    text(0.25, 0.95, textStr, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);

    if ss ==1
        ylabel('Percentage of pixels');
        xlabel('Percentage max value');
        a=gca();
        a.TickDir="out";
    else
        axis off
    end
    box off
end
