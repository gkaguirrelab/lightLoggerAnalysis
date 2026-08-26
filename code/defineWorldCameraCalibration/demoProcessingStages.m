
% Housekeeping
clear

% Identify a raw TIFF
tiffFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'exampleWorldCameraImages',...
    'indoor_1.tiff');

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

% Set any pixel that is saturated to inf
I = imageStages{2};
I(I>250)=Inf;
imageStages{3}=I;

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
load(paramFileName,'radiometricCorrectionMap');
imageStages{5} = imageStages{4} .* radiometricCorrectionMap;

% Plot
stages = {'raw','linearized','threshold','flattened','radiometric correction'};
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
    imshow(I,[]);
    title(stages{ss});

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
    text(0.65, 0.95, textStr, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);

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
