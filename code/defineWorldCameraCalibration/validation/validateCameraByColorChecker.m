% Validation with Macbeth Color Checker. We obtained an image of the
% Macbeth color checker with the IMX219 camera. We also measured the
% spectral radiance of 5 of the 24 patches using the PR670. We have
% available that tabular reflectance spectra of all 24 color squares. Usimg
% these data, we test the ability of the camera reconstruction pipeline to
% produce integrated radiance values for the R, G, and B channels of the
% camera that match what we predict based upon the direct spectral
% radiance measurements. To do so we:
%
% - Estimate the illuminant present in the scene
% - Combine the illuminant with the tabular spectral reflectance to obtain
%   the estimated spectral radiance for all 24 checks
% - Using the tabular spectral sensitivities of the IMX219 camera and the
%   estimated spectral radiances, obtain the predicted integrated radiance
%   for the three channels (RGB) for each of the 24 checks.
% - Convert the raw IMX219 image into a map of integrated spectral radiance
% - Identify the pixels within the IMX219 image that correspond to each of
%   the 24 Macbeth patches and obtain the mean integrated radiance values
%   for each of the color channels for each of the color patches
% - Compare these reconstructed integrated radiance values with the
%   predicted values.
%


%%%%%%%%%%%%%%%%%%%%
%% SETUP AND LOADING
%%%%%%%%%%%%%%%%%%%%


% Housekeeping. We clear all to make sure we have fresh persistent vals
clear all

% Indoor or outdoor data set?
valSetOptions = {'indoor','outdoor'};
valIdx = 1;

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing) for
% this analysis
commonS = [380, 1, 352]; 

% Set the rows and columns of the Macbeth chart
nRows = 4;
nColumns = 6;

% Using the "extractCheckerPixels" in the GUI mode, I defined the corner
% locations for the "close" camera image for each validation set.
cornerSets{1} = [
    128.6063  103.3571
    530.4668  109.7359
    530.4668  364.8854
    122.2276  376.5797];

cornerSets{2} = [
   73.3239  108.6728
  503.8887  103.3571
  542.1611  376.5797
   48.8721  419.1047];

% Get the list of spectral radiance measurements of checks
dirName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'macbethColorCheck',...
    valSetOptions{valIdx},...
    'PR670',...
    '*.mat');
fileList = dir(dirName);

% Load the spectral measurements and resample to 1 nm. Note that SplineSpd
% automatically handles the change in power due to the resampling from the
% 2 nm sampling originally (with units of W/m2/sr/wavelength band) to the 1
% nm sampling (with units which are now W/m2/sr/nm).
for ii = 1:length(fileList)
    fileName = fullfile(fileList(ii).folder,fileList(ii).name);
    load(fileName,'measurement','S')
    myIndex = int32(sscanf(fileList(ii).name, 'Index-%d'));
    [c, r] = ind2sub([nColumns nRows], myIndex);
    mySPD = mean(measurement,1);
    spectralRadiance{c, r} = SplineSpd(SToWls(S), mySPD', SToWls(commonS));    
end

% Get the table of reflectance spectra of the Macbeth color checker and
% then spline the reflectance to the commonS
[spectralReflectance,spectralReflectanceS] = loadMacbethReflectance();
for c = 1:nColumns
    for r = 1:nRows
            spectralReflectance{c, r} = SplineRaw(SToWls(spectralReflectanceS), spectralReflectance{c, r}(:), SToWls(commonS));
    end
end

% Load the world camera channel spectral sensitivity functions. This is a
% table with the first column providing the wavelength support.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');

% Extract wavelength support and sensor sensitivities from the table T
% Assuming the first column is Wavelength and the next three are R, G, B
sensorWls = T{:, 1};

% Resample spectral sensitivities to the common wavelength domain (380-730 nm)
T_common = SplineRaw(sensorWls, [T.red, T.green, T.blue], SToWls(commonS));

% Scale the camera spectral sensitivity functions so their maximum value is unity
sensorSensitivities = T_common ./ max(T_common, [], 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE THE ILLUMINANT
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize arrays to hold the estimated illuminant spectrum and its grid coordinates
illuminantEstimates = [];
measCols = [];
measRows = [];

% Loop over the columns and rows of the color checker
for c = 1:nColumns
    for r = 1:nRows
        if ~isempty(spectralRadiance{c, r})
            % Estimate effective illuminant and record its coordinate position
            illuminantEstimates(:, end+1) = spectralRadiance{c, r} ./ spectralReflectance{c, r};
            measCols(end+1, 1) = c;
            measRows(end+1, 1) = r;
        end
    end
end

% Fit a planar spatial model to the illuminant for each wavelength
% Model: Illuminant(lambda) = beta0(lambda) + beta1(lambda)*col + beta2(lambda)*row
% X_meas design matrix is [5 patches x 3 predictors (intercept, col, row)]
X_meas = [ones(length(measCols), 1), measCols, measRows];

% Solve for the coefficients across all wavelengths simultaneously
illuminantBetas = X_meas \ illuminantEstimates';

% We end up with some nans in the last entry of illuminantBetas from this
% regression. Not sure why. Replace these with the nearest value
illuminantBetas(:,end) = illuminantBetas(:,end-1);

% Evaluate if the spatial model is satisfactory by calculating R-squared
% model predictions for the 5 measured locations. Report the value.
predictedMeasIlluminants = X_meas * illuminantBetas;
ssTotal = sum((illuminantEstimates' - mean(illuminantEstimates', 1)).^2, 1);
ssResid = sum((illuminantEstimates' - predictedMeasIlluminants).^2, 1);
rSquared = 1 - (ssResid ./ ssTotal);
fprintf('Mean spatial model R-squared across all wavelengths: %.4f\n', mean(rSquared, 'omitnan'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE SPECTRAL RADIANCE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize storage arrays for the predicted and measured RGB radiance
predictedRGBRadiance = zeros(nRows, nColumns, 3);
measuredRGBRadiance = zeros(nRows, nColumns, 3);

% Loop over the columns and rows of the color checker
for c = 1:nColumns
    for r = 1:nRows

        % Obtain the modeled illuminant for this specific [col, row]
        % position
        localIlluminant = ([1, c, r] * illuminantBetas)'; 

        % Construct the spectral radiance if we did not measure it
        if isempty(spectralRadiance{c, r})
            spectralRadiance{c, r} = localIlluminant .* spectralReflectance{c, r}(:);
        end

        % Calculate the predicted RGB integrated radiance via dot product
        % of the source radiance and scaled world camera sensitivities
        predictedRGBRadiance(r, c, :) = (spectralRadiance{c, r}' * sensorSensitivities);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEASURED SPECTRAL RADIANCE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load the world camera lens intrinsics.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'arducamB0392cameraIntrinsics.mat');
load(dataFileName,'arducamB0392cameraIntrinsics');

% Load the data for the "close" camera acquisition of the color checker 
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'macbethColorCheck',...
    valSetOptions{valIdx},...
    'lightLogger',...
    'close_AGCandMS_01.mat');
load(dataFileName,'worldFrame','AGCSettings');

% Convert the raw worldFrame to a radiance map
radianceMap = reconstructionPipeline(worldFrame, AGCSettings);

% Demosaic the image so that we have an RGB triplet of radiance values for
% each pixel
radianceMap = demosaicRadianceMapRCD(radianceMap);

% Obtain the pixel indices within the world image for each check.
rawPixelIndices = extractCheckerPixels(worldFrame,arducamB0392cameraIntrinsics.results.Intrinsics,cornerSets{valIdx});

% Now loop through the rows and columns of the checker chart and obtain the
% measured RGB radiance values
for r = 1:nRows
    for c = 1:nColumns
        
        % Extract linear indices for the 75% central region of this check
        idx = rawPixelIndices{r, c};
        
        % Extract the individual channels from the demosaiced radiance map
        R_channel = radianceMap(:, :, 1);
        G_channel = radianceMap(:, :, 2);
        B_channel = radianceMap(:, :, 3);
        
        % Calculate the mean radiance for this check
        measuredRGBRadiance(r, c, 1) = mean(R_channel(idx), 'omitnan');
        measuredRGBRadiance(r, c, 2) = mean(G_channel(idx), 'omitnan');
        measuredRGBRadiance(r, c, 3) = mean(B_channel(idx), 'omitnan');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT PREDICTED VS MEASURED INTEGRATED RADIANCE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Predicted vs Measured Integrated Radiance', 'Position', [150, 150, 800, 600]);
hold on;

% Flatten the 4x6 matrices into 24x1 vectors for each channel
predR = reshape(predictedRGBRadiance(:, :, 1), [], 1);
predG = reshape(predictedRGBRadiance(:, :, 2), [], 1);
predB = reshape(predictedRGBRadiance(:, :, 3), [], 1);

measR = reshape(measuredRGBRadiance(:, :, 1), [], 1);
measG = reshape(measuredRGBRadiance(:, :, 2), [], 1);
measB = reshape(measuredRGBRadiance(:, :, 3), [], 1);

% Create an alpha map (0.1 for estimated, 0.5 for measured) and flatten it
alphaMap = 0.75 * ones(nRows, nColumns);
for i = 1:length(measCols)
    alphaMap(measRows(i), measCols(i)) = 0.75;
end
alphaFlat = reshape(alphaMap, [], 1);

% Plot each channel with a distinct color and literal transparency values
sR = scatter(predR, measR, 75, 'r', 'filled', 'MarkerEdgeColor', 'none', 'DisplayName', 'Red Channel');
sR.AlphaData = alphaFlat;
sR.MarkerFaceAlpha = 'flat';
sR.AlphaDataMapping = 'none';

sG = scatter(predG, measG, 75, 'g', 'filled', 'MarkerEdgeColor', 'none', 'DisplayName', 'Green Channel');
sG.AlphaData = alphaFlat;
sG.MarkerFaceAlpha = 'flat';
sG.AlphaDataMapping = 'none';

sB = scatter(predB, measB, 75, 'b', 'filled', 'MarkerEdgeColor', 'none', 'DisplayName', 'Blue Channel');
sB.AlphaData = alphaFlat;
sB.MarkerFaceAlpha = 'flat';
sB.AlphaDataMapping = 'none';

% Determine the axis limits based on the data to create a proportional plot
maxVal = max([predR; predG; predB; measR; measG; measB]) * 1.05;
minVal = min([predR; predG; predB; measR; measG; measB; 0]);

% Plot the 1:1 identity line
plot([minVal, maxVal], [minVal, maxVal], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Agreement');

% Calculate overall R-squared for all channels combined
allPred = [predR; predG; predB];
allMeas = [measR; measG; measB];

% Remove any potential NaN values (e.g., if a patch mask entirely failed) before fitting
validIdx = ~isnan(allPred) & ~isnan(allMeas);
mdl = fitlm(allPred(validIdx), allMeas(validIdx));
rSq = mdl.Rsquared.Ordinary;

% Formatting
axis equal;
xlim([minVal, maxVal]);
ylim([minVal, maxVal]);
xlabel('Predicted integrated radiance [W/m2/sr]');
ylabel('Measured integrated radiance [W/m2/sr]');
title(sprintf('Camera Validation: Predicted vs. Measured RGB Radiance\nOverall R^2 = %.4f', rSq));
legend('Location', 'northwest');
a = gca();
a.XTick = a.YTick;
grid on;
box on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT PREDICTED VS MEASURED CHANNEL RATIOS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the total radiance across all channels for each patch
predSum = predR + predG + predB;
measSum = measR + measG + measB;

% Normalize each channel to obtain the relative ratio (chromaticity)
predRatioR = predR ./ predSum;
predRatioG = predG ./ predSum;
predRatioB = predB ./ predSum;

measRatioR = measR ./ measSum;
measRatioG = measG ./ measSum;
measRatioB = measB ./ measSum;

% Create a new figure for the ratio agreement plot
figure('Name', 'Predicted vs Measured RGB Ratios', 'Position', [150, 150, 800, 600]);
hold on;

% Plot each channel ratio with a distinct color
scatter(predRatioR, measRatioR, 75, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Red Ratio');
scatter(predRatioG, measRatioG, 75, 'g', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Green Ratio');
scatter(predRatioB, measRatioB, 75, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Blue Ratio');

% Plot the 1:1 identity line for ratios (which strictly bound between 0 and 1)
plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Agreement');

% Calculate overall R-squared for the relative ratios
allPredRatios = [predRatioR; predRatioG; predRatioB];
allMeasRatios = [measRatioR; measRatioG; measRatioB];

% Remove any potential NaN values before fitting
validRatioIdx = ~isnan(allPredRatios) & ~isnan(allMeasRatios);
mdlRatio = fitlm(allPredRatios(validRatioIdx), allMeasRatios(validRatioIdx));
rSqRatio = mdlRatio.Rsquared.Ordinary;

% Formatting
axis equal;
xlim([0, 1]);
ylim([0, 1]);
xlabel('Predicted Channel Ratio (Channel / (R+G+B))');
ylabel('Measured Channel Ratio (Channel / (R+G+B))');
title(sprintf('Camera Validation: Predicted vs. Measured RGB Ratios\nOverall R^2 = %.4f', rSqRatio));
legend('Location', 'northwest');
grid on;
box on;
hold off;