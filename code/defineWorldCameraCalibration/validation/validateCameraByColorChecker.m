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
close all

% Define the save directory on the Desktop and ensure it exists
desktopPath = fullfile(getenv('HOME'), 'Desktop');
saveDir = fullfile(desktopPath, 'Macbeth_Validation_Figures');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Indoor or outdoor data set?
valSetOptions = {'indoor','indoor','indoor','outdoor','outdoor','outdoor'};
valNumOptions = {'1','2','3','1','2','3'};

% Plot the measurements in order of decreasing irradiance
plotOrder = [4,5,1,2,6,3];

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing) for
% this analysis
commonS = [380, 1, 352];

% Set the rows and columns of the Macbeth chart
nRows = 4;
nColumns = 6;

% Using the "extractCheckerPixels" in the GUI mode, I defined the corner
% locations for the "close" camera image for each validation set.
cornerSets{1} = [
    200.8987  153.3239
    406.0814  154.3870
    405.0183  284.0880
    197.7093  287.2774
    ];

cornerSets{2} = [
    150.9319  166.0814
    457.1113  176.7126
    484.7525  386.1478
    99.9020  398.9053
    ];

cornerSets{3} = [
    217.9086  160.7658
    449.6694  179.9020
    440.1013  335.1179
    197.7093  315.9817
    ];

cornerSets{4} = [
    198.7724  194.7857
    418.8389  191.5963
    430.5332  339.3704
    186.0150  344.6860
    ];

cornerSets{5} = [
    233.8555  234.1213
    430.5332  211.7957
    447.5432  342.5598
    242.3605  370.2010
    ];

cornerSets{6} = [
    134.9850  128.8721
    468.8056   83.1578
    485.8156  287.2774
    184.9518  351.0648
    ];

% Define standard sRGB reference colors for the 24 Macbeth patches (4x6 layout)
macbethRGB = zeros(nRows, nColumns, 3);
macbethRGB(1,:,:) = [115,82,68; 194,150,130; 98,122,157; 87,108,67; 133,128,177; 103,189,170]; % Row 1
macbethRGB(2,:,:) = [214,126,44; 80,91,166; 193,90,99; 94,60,108; 157,188,64; 224,163,46];    % Row 2
macbethRGB(3,:,:) = [56,61,150; 70,148,73; 175,54,60; 231,199,31; 187,86,149; 8,133,161];     % Row 3
macbethRGB(4,:,:) = [243,243,242; 200,200,200; 160,160,160; 122,122,121; 85,85,85; 52,52,52]; % Row 4 (Neutrals)
macbethRGB = macbethRGB / 255; % Normalize to [0, 1] for MATLAB plots

% Initialize storage arrays for the aggregate plot across all validation sets
allPredMean = cell(length(plotOrder), 1);
allMeasMean = cell(length(plotOrder), 1);

%% Loop over the validation sets

for vv = 1:length(plotOrder)

    valIdx = plotOrder(vv);

    % Get the list of spectral radiance measurements of checks
    dirName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'macbethColorCheck',...
        valSetOptions{valIdx},...
        valNumOptions{valIdx},...
        'PR670',...
        '*.mat');
    fileList = dir(dirName);

    % Load the spectral measurements and resample to 1 nm. Note that SplineSpd
    % automatically handles the change in power due to the resampling from the
    % 2 nm sampling originally (with units of W/m2/sr/wavelength band) to the 1
    % nm sampling (with units which are now W/m2/sr/nm).
    spectralRadiance = [];
    for ii = 1:length(fileList)
        fileName = fullfile(fileList(ii).folder,fileList(ii).name);
        load(fileName,'measurement','S')
        patchLocation = regexp(fileList(ii).name, '_R(\d+)C(\d+)\.mat$', 'tokens', 'once');
        if isempty(patchLocation)
            error('Could not parse ColorChecker row and column from %s', fileList(ii).name);
        end
        r = str2double(patchLocation{1});
        c = str2double(patchLocation{2});
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
        valNumOptions{valIdx},...
        'lightLogger',...
        'close_AGCandMS_01.mat');
    load(dataFileName,'worldFrame','AGCSettings','minispectValue');

    % Convert the raw worldFrame to a radiance map
    [radianceMap, imageStages] = reconstructionPipeline(worldFrame, AGCSettings);

    % Demosaic the image so that we have an RGB triplet of radiance values for
    % each pixel
    radianceMap = demosaicRadianceMap(radianceMap);

    % Obtain the pixel indices within the world image for each check. We
    % scale the image by the Dgain so that it is visible if we wish to
    % display it.
    rawPixelIndices = extractCheckerPixels(worldFrame*AGCSettings.Dgain,arducamB0392cameraIntrinsics.results.Intrinsics,cornerSets{valIdx});

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

    % Obtain the spectral radiance estimated from the minispect
    [minispectSPD,minispectS,fVal,fitErrors] = estimateRadianceSpectrumFromMinispect(minispectValue.AS);
    minispectSPD = SplineSpd(SToWls(minispectS), minispectSPD, SToWls(commonS));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name', sprintf('Validation: %s', vv), ...
        'WindowStyle', 'Docked',...
        'Position', [50+(vv-1)*100, 50+(vv-1)*100, 1600, 400]);
    tiledlayout(1,4,"TileSpacing","compact","Padding","tight");

    % Show the camera image
    nexttile
    logImage = log10(radianceMap);
    logImage = logImage-min(logImage(:));
    logImage = logImage/max(logImage(:));
    imagesc(logImage)
    box off
    axis off
    axis equal

    % Show the illuminant
    nexttile
    plot(SToWls(commonS),mean(predictedMeasIlluminants),'-k','LineWidth',2);
    ylabel('radiance [W/m2/sr/nm]');
    xlabel('wavelength [nm]');
    title('Illuminant');
    axis square

    % Show the mean integrated radiance with true-color patches
    nexttile; hold on;

    predMean = mean(predictedRGBRadiance, 3);
    measMean = mean(measuredRGBRadiance, 3);
    
    % Store the current predictions and measurements for the final aggregate plot
    allPredMean{vv} = predMean;
    allMeasMean{vv} = measMean;

    for r = 1:nRows
        for c = 1:nColumns
            patchColor = squeeze(macbethRGB(r, c, :))';
            scatter(predMean(r, c), measMean(r, c), 85, patchColor, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end

    % Store the slope and intercept of a robust fit
    predMeasFit(vv,:) = robustfit(predMean(:),measMean(:));
    AgainStore(vv) = AGCSettings.Again;
    DgainStore(vv) = AGCSettings.Dgain;

    % Determine axis limits
    maxValMean = max([predMean(:); measMean(:)]) * 1.05;
    minValMean = min([predMean(:); measMean(:); 0]);
    plot([minValMean, maxValMean], [minValMean, maxValMean], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    axis square;
    xlim([minValMean, maxValMean]);
    ylim([minValMean, maxValMean]);
    xlabel('via PR670 [W/m2/sr]');
    ylabel('via IMX219 [W/m2/sr]');
    title('Integrated Radiance');
    a = gca();
    a.XTick = a.YTick;
    grid on;
    box on;
    hold off;

    % Show the R:G:B ratio agreement via channel fractions
    nexttile; hold on;

    % Calculate total radiance per patch to extract fractional channel ratios
    predSum = sum(predictedRGBRadiance, 3);
    measSum = sum(measuredRGBRadiance, 3);

    predFracR = reshape(predictedRGBRadiance(:,:,1) ./ predSum, [], 1);
    predFracG = reshape(predictedRGBRadiance(:,:,2) ./ predSum, [], 1);
    predFracB = reshape(predictedRGBRadiance(:,:,3) ./ predSum, [], 1);

    measFracR = reshape(measuredRGBRadiance(:,:,1) ./ measSum, [], 1);
    measFracG = reshape(measuredRGBRadiance(:,:,2) ./ measSum, [], 1);
    measFracB = reshape(measuredRGBRadiance(:,:,3) ./ measSum, [], 1);

    % Plot each fractional channel color contribution
    scatter(predFracR, measFracR, 85, 'r', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
    scatter(predFracG, measFracG, 85, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
    scatter(predFracB, measFracB, 85, 'b', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);

    % 1:1 identity line
    plot([0, 0.75], [0, 0.75], 'k--', 'LineWidth', 1.5);

    a = gca();
    a.XTick = [0 0.25 0.5 0.75];
    a.YTick = [0 0.25 0.5 0.75];

    axis square;
    xlim([0, 0.75]);
    ylim([0, 0.75]);
    xlabel('Predicted Channel Fraction');
    ylabel('Measured Channel Fraction');
    title('R:G:B Ratio');
    grid on;
    box on;
    hold off;

    % Draw the figure to ensure it is fully rendered before saving
    drawnow;

    % Construct a descriptive filename based on the environment and set number
    envStr = valSetOptions{valIdx};
    numStr = valNumOptions{valIdx};

    fileName = sprintf('Validation_%s_set%s_NDF_Dgain_%1.2f.pdf', envStr, numStr, AGCSettings.Dgain);
    fullSavePath = fullfile(saveDir, fileName);
    exportgraphics(gcf, fullSavePath, 'ContentType', 'vector');

    % Export the figure as a PNG
    %{
    fileName = sprintf('Validation_%s_set%s_NDF_Dgain_%1.2f.png', envStr, numStr, AGCSettings.Dgain);
    fullSavePath = fullfile(saveDir, fileName);
    exportgraphics(gcf, fullSavePath);
    %}

    fprintf('Saved figure to: %s\n', fullSavePath);

end % loop over valSetOptions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AGGREGATE PLOT: ALL VALIDATION SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'All Validation Sets - Integrated Radiance', ...
    'WindowStyle', 'Docked', ...
    'Position', [50, 50, 600, 600]);
hold on;

% 1:1 identity line
plot([-4 0],[-4 0], 'k--', 'LineWidth', 1.5);

% Define distinct marker symbols for the different validation sets
markers = {'o', 's', '^', 'd', 'v', 'p'};
maxValMeanAgg = 0;
minValMeanAgg = inf;

for vv = 1:length(plotOrder)
    pMean = allPredMean{vv};
    mMean = allMeasMean{vv};
    
    % Update absolute max and min for the axis limits
    maxValMeanAgg = max([maxValMeanAgg; pMean(:); mMean(:)]);
    minValMeanAgg = min([minValMeanAgg; pMean(:); mMean(:); 0]);

    for r = 1:nRows
        for c = 1:nColumns
            patchColor = squeeze(macbethRGB(r, c, :))';
            scatter(log10(pMean(r, c)), log10(mMean(r, c)), 85, patchColor, 'filled', ...
                'Marker', markers{vv}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end
end


axis square;
xlim([-4 0]);
ylim([-4 0]);
xlabel('via PR670 [log10 W/m2/sr]');
ylabel('via IMX219 [log10 W/m2/sr]');
title('Integrated radiance');
a = gca();
a.XTick = a.YTick;
a.TickDir = 'out';
grid on;
box on;

% Create custom handles to generate a legend for the marker types
dummyPlots = gobjects(length(plotOrder), 1);
legendLabels = cell(length(plotOrder), 1);
for vv = 1:length(plotOrder)
    valIdx = plotOrder(vv);
    % Plot invisible points mapped to the marker types
    dummyPlots(vv) = scatter(NaN, NaN, 85, [0.5 0.5 0.5], 'filled', ...
        'Marker', markers{vv}, 'MarkerEdgeColor', 'k');

    switch valSetOptions{valIdx}
        case 'indoor'
            setting = 'in';
        case 'outdoor'
            setting = 'out';
    end

    % Update legend label to include environment and gains, excluding the index number
    legendLabels{vv} = sprintf('%s, A: %2.2f, D: %2.2f', ...
        setting, AgainStore(vv), DgainStore(vv));
end
legend(dummyPlots, legendLabels, 'Location','southeast');


hold off;
drawnow;

% Save aggregate figure
fileNameAgg = 'Validation_All_Sets_Integrated_Radiance.pdf';
fullSavePathAgg = fullfile(saveDir, fileNameAgg);
exportgraphics(gcf, fullSavePathAgg, 'ContentType', 'vector');
fprintf('Saved aggregate figure to: %s\n', fullSavePathAgg);
