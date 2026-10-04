% defineAGCToIntegratedRadianceViaMacbeth.m
% The purpose of this script is to define the relationship between the
% custom AGC settings we use to control the sensitivity of the IMX219
% camera and the true effective integrated radiance of the environment
% as seen by the R, G, and B channels independently.
%
% This version uses a set of Macbeth color checker validation measurements
% to optimize the agcToRadianceP [slope, intercept] mapping parameters by
% minimizing the sum of squared errors between predicted and measured 
% radiance in log10 space.

% Housekeeping
clear
close all

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing)
commonS = [380, 1, 352];
commonWls = SToWls(commonS);

% Define the datasets and plot order
valSetOptions = {'indoor','indoor','indoor','outdoor','outdoor','outdoor'};
valNumOptions = {'1','2','3','1','2','3'};
plotOrder = [4,5,1,2,6,3]; % Order of decreasing irradiance

% Set the rows and columns of the Macbeth chart
nRows = 4;
nColumns = 6;

% Corner locations for the "close" camera image for each validation set
cornerSets{1} = [200.8987 153.3239; 406.0814 154.3870; 405.0183 284.0880; 197.7093 287.2774];
cornerSets{2} = [150.9319 166.0814; 457.1113 176.7126; 484.7525 386.1478; 99.9020 398.9053];
cornerSets{3} = [217.9086 160.7658; 449.6694 179.9020; 440.1013 335.1179; 197.7093 315.9817];
cornerSets{4} = [198.7724 194.7857; 418.8389 191.5963; 430.5332 339.3704; 186.0150 344.6860];
cornerSets{5} = [233.8555 234.1213; 430.5332 211.7957; 447.5432 342.5598; 242.3605 370.2010];
cornerSets{6} = [134.9850 128.8721; 468.8056  83.1578; 485.8156 287.2774; 184.9518 351.0648];

% Define standard sRGB reference colors for the 24 Macbeth patches
macbethRGB = zeros(nRows, nColumns, 3);
macbethRGB(1,:,:) = [115,82,68; 194,150,130; 98,122,157; 87,108,67; 133,128,177; 103,189,170];
macbethRGB(2,:,:) = [214,126,44; 80,91,166; 193,90,99; 94,60,108; 157,188,64; 224,163,46]; 
macbethRGB(3,:,:) = [56,61,150; 70,148,73; 175,54,60; 231,199,31; 187,86,149; 8,133,161];  
macbethRGB(4,:,:) = [243,243,242; 200,200,200; 160,160,160; 122,122,121; 85,85,85; 52,52,52]; 
macbethRGB = macbethRGB / 255; 

% Load parameters required for set point calculation
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'darkSignal.mat'), 'darkSignal');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'nonLinearClippingExponent.mat'), 'clippingExponent');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'flatFieldingFunction.mat'), 'correctionMap');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'radiometricCorrectionRGB.mat'), 'radiometricCorrectionMap');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'arducamB0392cameraIntrinsics.mat'), 'arducamB0392cameraIntrinsics');

Smax = 2^8 - 1 - darkSignal;
n = clippingExponent;
meanCorrectionFielding = 1 / mean(1 ./ correctionMap(:), 'omitnan');
meanCorrectionRGB = 1 / mean(1 ./ radiometricCorrectionMap(:), 'omitnan');

% Preallocate storage for optimization
allPredicted = [];
allUnscaledMeasured = [];
allCameraScores = [];

% Storage for plotting loops
predictedRGBRadianceStore = cell(length(plotOrder), 1);
unscaledMeasuredRGBRadianceStore = cell(length(plotOrder), 1);
cameraScoresStore = zeros(length(plotOrder), 1);
AgainStore = zeros(length(plotOrder), 1);
DgainStore = zeros(length(plotOrder), 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA EXTRACTION (NO OPTIMIZATION YET)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Extracting patch data across validation sets...\n');

% Load the IMX219 sensitivity functions
dataFileName = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');
wlsSensor = T.wls;
T_common = SplineRaw(wlsSensor, [T.red, T.green, T.blue], commonWls);
sensorSensitivities = T_common ./ max(T_common, [], 1);

for vv = 1:length(plotOrder)
    valIdx = plotOrder(vv);

    % Get the list of spectral radiance measurements of checks
    dirName = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'macbethColorCheck', ...
        valSetOptions{valIdx}, valNumOptions{valIdx}, 'PR670', '*.mat');
    fileList = dir(dirName);

    spectralRadiance = cell(nColumns, nRows);
    for ii = 1:length(fileList)
        fileName = fullfile(fileList(ii).folder,fileList(ii).name);
        load(fileName,'measurement','S')
        patchLocation = regexp(fileList(ii).name, '_R(\d+)C(\d+)\.mat$', 'tokens', 'once');
        r = str2double(patchLocation{1});
        c = str2double(patchLocation{2});
        mySPD = mean(measurement,1);
        spectralRadiance{c, r} = SplineSpd(SToWls(S), mySPD', commonWls);
    end

    [spectralReflectance,spectralReflectanceS] = loadMacbethReflectance();
    for c = 1:nColumns
        for r = 1:nRows
            spectralReflectance{c, r} = SplineRaw(SToWls(spectralReflectanceS), spectralReflectance{c, r}(:), commonWls);
        end
    end

    % Estimate Illuminant
    illuminantEstimates = [];
    measCols = []; measRows = [];
    for c = 1:nColumns
        for r = 1:nRows
            if ~isempty(spectralRadiance{c, r})
                illuminantEstimates(:, end+1) = spectralRadiance{c, r} ./ spectralReflectance{c, r};
                measCols(end+1, 1) = c; measRows(end+1, 1) = r;
            end
        end
    end

    X_meas = [ones(length(measCols), 1), measCols, measRows];
    illuminantBetas = X_meas \ illuminantEstimates';
    illuminantBetas(:,end) = illuminantBetas(:,end-1); % handle nans

    % Predicted Radiance
    predictedRGBRadiance = zeros(nRows, nColumns, 3);
    for c = 1:nColumns
        for r = 1:nRows
            localIlluminant = ([1, c, r] * illuminantBetas)';
            if isempty(spectralRadiance{c, r})
                spectralRadiance{c, r} = localIlluminant .* spectralReflectance{c, r}(:);
            end
            predictedRGBRadiance(r, c, :) = (spectralRadiance{c, r}' * sensorSensitivities);
        end
    end
    predictedRGBRadianceStore{vv} = predictedRGBRadiance;

    % Load IMX219 Data
    dataFileName = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'macbethColorCheck', ...
        valSetOptions{valIdx}, valNumOptions{valIdx}, 'lightLogger', 'close_AGCandMS_01.mat');
    load(dataFileName,'worldFrame','AGCSettings');
    AgainStore(vv) = AGCSettings.Again;
    DgainStore(vv) = AGCSettings.Dgain;

    % Obtain unscaled radiance directly from the pipeline stages to bypass
    % the internally loaded agcToRadianceP file. Stage 5 is radiometrically equalized.
    [~, imageStages] = reconstructionPipeline(worldFrame, AGCSettings);
    
    setPoint = 127;
    setPoint = (setPoint / AGCSettings.Dgain) - darkSignal;
    linSetPoint = setPoint / (1 - (setPoint / Smax)^n)^(1/n);
    linSetPoint = linSetPoint * meanCorrectionFielding * meanCorrectionRGB;
    
    % This isolates the linear camera response without polynomial scaling
    radianceMapUnscaled = imageStages{5} / linSetPoint;
    radianceMapUnscaled = demosaicRadianceMap(radianceMapUnscaled);

    % Derive the effective camera score
    thisCameraScore = (AGCSettings.exposure * AGCSettings.Again) / linSetPoint;
    cameraScoresStore(vv) = thisCameraScore;
    
    % Extract checking pixels
    rawPixelIndices = extractCheckerPixels(worldFrame*AGCSettings.Dgain, arducamB0392cameraIntrinsics.results.Intrinsics, cornerSets{valIdx});
    measuredRGBRadianceUnscaled = zeros(nRows, nColumns, 3);
    for r = 1:nRows
        for c = 1:nColumns
            idx = rawPixelIndices{r, c};
            R_channel = radianceMapUnscaled(:, :, 1);
            G_channel = radianceMapUnscaled(:, :, 2);
            B_channel = radianceMapUnscaled(:, :, 3);
            measuredRGBRadianceUnscaled(r, c, 1) = mean(R_channel(idx), 'omitnan');
            measuredRGBRadianceUnscaled(r, c, 2) = mean(G_channel(idx), 'omitnan');
            measuredRGBRadianceUnscaled(r, c, 3) = mean(B_channel(idx), 'omitnan');
        end
    end
    unscaledMeasuredRGBRadianceStore{vv} = measuredRGBRadianceUnscaled;

    % Aggregate for optimization
    allPredicted = [allPredicted; predictedRGBRadiance(:)];
    allUnscaledMeasured = [allUnscaledMeasured; measuredRGBRadianceUnscaled(:)];
    allCameraScores = [allCameraScores; repmat(thisCameraScore, numel(predictedRGBRadiance), 1)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMIZE agcToRadianceP PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Optimizing agcToRadianceP...\n');

% Objective function: minimize sum of squared errors in log10 space
% predicted = measured_unscaled * 10^polyval(p, log10(cameraScore))
objFun = @(p) sum((log10(allPredicted) - log10(allUnscaledMeasured .* 10.^polyval(p, log10(allCameraScores)))).^2);

% Run fminsearch (start with theoretical -1.0 slope and empirical guess)
initialGuess = [-1.0, 2.0];
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
agcToRadianceP = fminsearch(objFun, initialGuess, options);

fprintf('==================================================\n');
fprintf('Optimized agcToRadianceP: [%.4f, %.4f]\n', agcToRadianceP(1), agcToRadianceP(2));
fprintf('==================================================\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AGGREGATE PLOT: ALL VALIDATION SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'All Validation Sets - Integrated Radiance', 'WindowStyle', 'Docked', 'Position', [50, 50, 600, 600]);
hold on;
plot([-4 0],[-4 0], 'k--', 'LineWidth', 1.5);

markers = {'o', 's', '^', 'd', 'v', 'p'};
for vv = 1:length(plotOrder)
    
    % Scale the unscaled measured values with our newly optimized polynomial
    optimalMeanIntegratedRadiance = 10.^polyval(agcToRadianceP, log10(cameraScoresStore(vv)));
    mMean = mean(unscaledMeasuredRGBRadianceStore{vv} * optimalMeanIntegratedRadiance, 3);
    pMean = mean(predictedRGBRadianceStore{vv}, 3);
    
    for r = 1:nRows
        for c = 1:nColumns
            patchColor = squeeze(macbethRGB(r, c, :))';
            scatter(log10(pMean(r, c)), log10(mMean(r, c)), 85, patchColor, 'filled', ...
                'Marker', markers{vv}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end
end

axis square;
xlim([-4 0]); ylim([-4 0]);
xlabel('via PR670 [log_{10} W/m^2/sr]'); ylabel('via IMX219 [log_{10} W/m^2/sr]');
title('Integrated Radiance vs. Camera AGC Sensitivity (Optimized)');
a = gca(); a.XTick = a.YTick; a.TickDir = 'out';
grid on; box on;

dummyPlots = gobjects(length(plotOrder), 1);
legendLabels = cell(length(plotOrder), 1);
for vv = 1:length(plotOrder)
    valIdx = plotOrder(vv);
    dummyPlots(vv) = scatter(NaN, NaN, 85, [0.5 0.5 0.5], 'filled', 'Marker', markers{vv}, 'MarkerEdgeColor', 'k');

    switch valSetOptions{valIdx}
        case 'indoor'
            setting = 'in';
        case 'outdoor'
            setting = 'out';
    end
    legendLabels{vv} = sprintf('%s, A: %2.2f, D: %2.2f', setting, AgainStore(vv), DgainStore(vv));
end
legend(dummyPlots, legendLabels, 'Location','southeast');
hold off; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE DERIVED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the values that relate camera score to channel-specific effective radiance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'cameraScoreToIntegratedRadiance.mat');
readme = ['Created by defineAGCToIntegratedRadiance.\n'...
    'A linear function (in log10 space) maps AGC values to integrated radiance.\n',...
    'Optimized utilizing the Macbeth Color Checker validation routines.\n',...
    'agcToRadianceP -- the slope and intercept.\n'];
save(saveFileName,'readme','agcToRadianceP');
fprintf('Saved derived values to: %s\n', saveFileName);