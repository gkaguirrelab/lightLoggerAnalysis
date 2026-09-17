% Validation

% Housekeeping
clear all

% Get the list of spectral radiometric measurements of checks
dirName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'macbethColorCheck',...
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
    [col, row] = ind2sub([6 4], myIndex);
    newS = S; newS(2) = 1; newS(3) = S(3)*2;
    mySPD = mean(measurement,1);
    spectralRadiance{col, row} = SplineSpd(SToWls(S), mySPD', SToWls(newS));    
end
spectralRadianceS = newS;

% Get the table of reflectance spectra of the macbeth color checker
[spectralReflectance,spectralReflectanceS] = loadMacbethReflectance();

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing).
commonS = [380, 1, 352]; 

% Initialize arrays to hold the estimated illuminant spectrum and its grid coordinates
illuminantEstimates = [];
measCols = [];
measRows = [];

% 1. Estimate the illuminant at the measured locations
for c = 1:6
    for r = 1:4
        if ~isempty(spectralRadiance{c, r})

            % Resample both radiance and reflectance to the common wavelength domain
            rad_common = SplineSpd(SToWls(spectralRadianceS), spectralRadiance{c, r}(:), SToWls(commonS));
            ref_common = SplineRaw(SToWls(spectralReflectanceS), spectralReflectance{c, r}(:), SToWls(commonS));

            % Estimate effective illuminant and record its coordinate position
            illuminantEstimates(:, end+1) = rad_common ./ ref_common;
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
% Model predictions for the 5 measured locations
predictedMeasIlluminants = X_meas * illuminantBetas;
ssTotal = sum((illuminantEstimates' - mean(illuminantEstimates', 1)).^2, 1);
ssResid = sum((illuminantEstimates' - predictedMeasIlluminants).^2, 1);
rSquared = 1 - (ssResid ./ ssTotal);
fprintf('Mean spatial model R-squared across all wavelengths: %.4f\n', mean(rSquared, 'omitnan'));

% 2. Predict spectral radiance for unmeasured patches using the spatial model
for c = 1:6
    for r = 1:4
        % Predict the illuminant for this specific [col, row] position
        localIlluminant = ([1, c, r] * illuminantBetas)'; 

        if isempty(spectralRadiance{c, r})
            ref_patch = spectralReflectance{c, r}(:);
            ref_common = SplineRaw(SToWls(spectralReflectanceS), ref_patch, SToWls(commonS));
            spectralRadiance{c, r} = localIlluminant .* ref_common;
        else
            rad_common = SplineSpd(SToWls(spectralRadianceS), spectralRadiance{c, r}(:), SToWls(commonS));
            spectralRadiance{c, r} = rad_common;
        end
    end
end

% Update the radiance S vector to reflect the new common wavelength domain
spectralRadianceS = commonS;

% Create a new figure sized for a 6x4 grid
figure('Name', 'Spatial Illuminant Model vs Measurements', 'Position', [100, 100, 1400, 800]);

% Get wavelength support for the x-axis
wls = SToWls(commonS);

% Determine global y-axis limits to ensure consistent scaling across all subplots
% Evaluate the model at all 24 positions to find the maximum predicted value
allCols = repmat(1:6, 1, 4)';
allRows = kron(1:4, ones(1, 6))';
allPredicted = [ones(24, 1), allCols, allRows] * illuminantBetas;
maxY = max(allPredicted(:)) * 1.1;

% Ensure the empirical measurements don't exceed the calculated maxY
if exist('illuminantEstimates', 'var') && ~isempty(illuminantEstimates)
    maxY = max([maxY, max(illuminantEstimates(:)) * 1.1]);
end

% Loop through the 4 rows and 6 columns of the Macbeth checker
for r = 1:4
    for c = 1:6
        % Calculate subplot index (1 to 24, moving row by row)
        plotIdx = (r - 1) * 6 + c;
        subplot(4, 6, plotIdx);
        hold on;

        % 1. Plot the modeled (estimated) illuminant
        % Calculate beta0 + beta1*c + beta2*r
        estIlluminant = ([1, c, r] * illuminantBetas)';
        plot(wls, estIlluminant, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Model');

        % 2. Plot the measured illuminant if available for this patch
        % Check if the current (c, r) coordinate exists in the measured data
        idx = find(measCols == c & measRows == r);
        if ~isempty(idx)
            measIlluminant = illuminantEstimates(:, idx);
            plot(wls, measIlluminant, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Measured');
        end

        % Formatting
        title(sprintf('Col %d, Row %d', c, r));
        xlim([min(wls), max(wls)]);
        ylim([0, maxY]);
        grid on;

        % Add legend only to the first subplot
        if plotIdx == 1
            legend('Location', 'best');
        end

        % Clean up axis labels for a cleaner grid layout
        if r == 4
            xlabel('Wavelength (nm)');
        else
            set(gca, 'XTickLabel', []);
        end

        if c == 1
            ylabel('Radiance');
        else
            set(gca, 'YTickLabel', []);
        end
    end
end

% Create a figure to compare measured and modeled spectral radiance
figure('Name', 'Measured vs Modeled Spectral Radiance', 'Position', [200, 200, 700, 500]);
hold on;

% Select the first available measured patch to plot
if ~isempty(measCols) && ~isempty(measRows)
    c = measCols(1);
    r = measRows(1);

    % Get wavelength support from the common domain
    wls = SToWls(commonS);

    % 1. Retrieve the Measured Radiance
    measRad = spectralRadiance{c, r};

    % 2. Calculate the Modeled Radiance
    localIlluminant = ([1, c, r] * illuminantBetas)';

    % Resample the patch reflectance to the common domain
    ref_patch = spectralReflectance{c, r}(:);
    ref_common = SplineRaw(SToWls(spectralReflectanceS), ref_patch, SToWls(commonS));

    % Calculate modeled radiance and convert to W/m2/sr/2nm-band
    modRad = localIlluminant .* ref_common;

    % Plot both spectra
    plot(wls, measRad, 'r-', 'LineWidth', 2, 'DisplayName', 'Measured Radiance');
    plot(wls, modRad, 'k--', 'LineWidth', 2, 'DisplayName', 'Modeled Radiance');

    % Formatting
    title(sprintf('Spectral Radiance Match for Patch (Col %d, Row %d)', c, r));
    xlabel('Wavelength (nm)');
    ylabel('Radiance (W/m^2/sr/nm');
    legend('Location', 'best');
    grid on;
    box on;
else
    disp('No measured patches available to plot.');
end
hold off;


% Load the world camera lens intrinsics.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'arducamB0392cameraIntrinsics.mat');
load(dataFileName,'arducamB0392cameraIntrinsics');

% Load the world camera channel spectral sensitivity functions. This is a
% table with the first column providing the wavelength support.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');

% Load the data for the "close" camera acquisition of the color checker 
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'macbethColorCheck',...
    'lightLogger',...
    'close_AGCandMS_01.mat');
load(dataFileName,'worldFrame','AGCSettings');

% Convert the raw worldFrame to a radiance map
radianceMap = reconstructionPipeline(worldFrame, AGCSettings);

% Demosaic the image so that we have an RGB triplet of radiance values for
% each pixel
radianceMap = demosaicRadianceMapRCD(radianceMap);

% Using the "extractCheckerPixels" in the GUI mode, I defined the corner
% locations for the "close" camera image. We now call this routine again
% with the defined corners to obtain the locations of pixels within the
% image corresponding to each of the checks
corners = [
    128.6063  103.3571
    530.4668  109.7359
    530.4668  364.8854
    122.2276  376.5797];
rawPixelIndices = extractCheckerPixels(worldFrame,arducamB0392cameraIntrinsics.results.Intrinsics,true,corners);

% Now loop through the rows and columns of the checker chart and obtain the
% estimated RGB radiance values given the estimated illuminant at each
% location, the spectral reflectance of each check, and the spectral
% sensitivities of the world camera sensors.

% Extract wavelength support and sensor sensitivities from the table T
% Assuming the first column is Wavelength and the next three are R, G, B
sensorWls = T{:, 1};
sensorSensitivities = [T.red, T.green, T.blue];

% Resample spectral sensitivities to the common wavelength domain (380-730 nm)
T_common = SplineRaw(sensorWls, sensorSensitivities, SToWls(commonS));

% Scale the camera spectral sensitivity functions so their maximum value is unity
T_scaled = T_common ./ max(T_common, [], 1);

% Initialize storage arrays for the predicted and measured RGB radiance
predictedRGBRadiance = zeros(4, 6, 3);
measuredRGBRadiance = zeros(4, 6, 3);

% Loop through the rows (1-4) and columns (1-6) of the color checker
for r = 1:4
    for c = 1:6
        % --- 1. Calculate Predicted RGB Radiance ---
        % Reconstruct the estimated local illuminant for this check position
        localIlluminant = ([1, c, r] * illuminantBetas)';
        
        % Resample this check's reflectance to the common wavelength domain
        ref_patch = spectralReflectance{c, r}(:);
        ref_common = SplineRaw(SToWls(spectralReflectanceS), ref_patch, SToWls(commonS));
        
        % Calculate the source spectral radiance (W/m2/sr/nm)
        sourceRadiance = localIlluminant .* ref_common;
        
        % Calculate the predicted RGB radiance via dot product of the source 
        % radiance and scaled sensitivities
        predictedRGBRadiance(r, c, :) = (sourceRadiance' * T_scaled);
        
        % --- 2. Calculate Measured RGB Radiance ---
        % Extract linear indices for the 75% central region of this check[cite: 6]
        idx = rawPixelIndices{r, c};
        
        % Extract the individual channels from the 3D demosaiced radiance map[cite: 5]
        R_channel = radianceMap(:, :, 1);
        G_channel = radianceMap(:, :, 2);
        B_channel = radianceMap(:, :, 3);
        
        % Calculate the mean radiance for this check, ignoring any Inf/NaN values[cite: 5]
        measuredRGBRadiance(r, c, 1) = mean(R_channel(idx), 'omitnan');
        measuredRGBRadiance(r, c, 2) = mean(G_channel(idx), 'omitnan');
        measuredRGBRadiance(r, c, 3) = mean(B_channel(idx), 'omitnan');
    end
end

% Create a new figure for the agreement plot
figure('Name', 'Predicted vs Measured RGB Radiance', 'Position', [150, 150, 800, 600]);
hold on;

% Flatten the 4x6 matrices into 24x1 vectors for each channel
predR = reshape(predictedRGBRadiance(:, :, 1), [], 1);
predG = reshape(predictedRGBRadiance(:, :, 2), [], 1);
predB = reshape(predictedRGBRadiance(:, :, 3), [], 1);

measR = reshape(measuredRGBRadiance(:, :, 1), [], 1);
measG = reshape(measuredRGBRadiance(:, :, 2), [], 1);
measB = reshape(measuredRGBRadiance(:, :, 3), [], 1);

% Plot each channel with a distinct color
scatter(predR, measR, 75, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Red Channel');
scatter(predG, measG, 75, 'g', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Green Channel');
scatter(predB, measB, 75, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Blue Channel');

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
xlabel('Predicted RGB Radiance');
ylabel('Measured RGB Radiance');
title(sprintf('Camera Validation: Predicted vs. Measured RGB Radiance\nOverall R^2 = %.4f', rSq));
legend('Location', 'northwest');
grid on;
box on;
hold off;

% --- Append to validateCameraByColorChecker_2.m ---

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