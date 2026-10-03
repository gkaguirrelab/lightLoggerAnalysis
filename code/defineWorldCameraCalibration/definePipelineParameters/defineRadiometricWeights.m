% This script derives the multiplicative adjustments that should be applied
% to the R, G, and B channels so that the sensor values reflect the
% integrated radiance of the light source, correcting for the differing
% peak transmission efficiencies of the Bayer filters.
%
% It uses Macbeth color checker images and PR670 measurements to find the
% optimal radiometric weights by matching the predicted and measured 
% channel fractions across five validation sets. 

%%%%%%%%%%%%%%%%%%%%
%% SETUP AND LOADING
%%%%%%%%%%%%%%%%%%%%

% Housekeeping. We clear all to make sure we have fresh persistent vals
clear all
close all

% What is the Bayer pattern in these data?
bayerPattern = "BGGR";

% Indoor or outdoor data set?
valSetOptions = {'indoor','indoor','outdoor','outdoor','outdoor'};
valNumOptions = {'1','2','1','2','3'};

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
    198.7724  194.7857
    418.8389  191.5963
    430.5332  339.3704
    186.0150  344.6860
    ];

cornerSets{4} = [
    233.8555  234.1213
    430.5332  211.7957
    447.5432  342.5598
    242.3605  370.2010
    ];

cornerSets{5} = [
    134.9850  128.8721
    468.8056   83.1578
    485.8156  287.2774
    184.9518  351.0648
    ];

% Initialize an array to hold the adjustments for each validation set
allAdjustments = zeros(length(valSetOptions), 3);
imgSize = [480, 640]; % Default to be updated by the actual camera frames

%% Loop over the validation sets to calculate weight adjustments
for valIdx = 1:length(valSetOptions)

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

    % Load the spectral measurements and resample to 1 nm.
    spectralRadiance = cell(nColumns, nRows);
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

    % Load the world camera channel spectral sensitivity functions. 
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(dataFileName,'T');

    % Extract wavelength support and sensor sensitivities from the table T
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

    for c = 1:nColumns
        for r = 1:nRows
            if ~isempty(spectralRadiance{c, r})
                illuminantEstimates(:, end+1) = spectralRadiance{c, r} ./ spectralReflectance{c, r};
                measCols(end+1, 1) = c;
                measRows(end+1, 1) = r;
            end
        end
    end

    % Fit a planar spatial model to the illuminant for each wavelength
    X_meas = [ones(length(measCols), 1), measCols, measRows];
    illuminantBetas = X_meas \ illuminantEstimates';

    % Replace any nans in the last entry of illuminantBetas
    illuminantBetas(:,end) = illuminantBetas(:,end-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESTIMATE SPECTRAL RADIANCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    predictedRGBRadiance = zeros(nRows, nColumns, 3);
    measuredRGBRadiance = zeros(nRows, nColumns, 3);

    for c = 1:nColumns
        for r = 1:nRows
            localIlluminant = ([1, c, r] * illuminantBetas)';
            if isempty(spectralRadiance{c, r})
                spectralRadiance{c, r} = localIlluminant .* spectralReflectance{c, r}(:);
            end
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
    load(dataFileName,'worldFrame','AGCSettings');
    
    % Update global imgSize dynamically
    imgSize = size(worldFrame);

    % Convert the raw worldFrame to a radiance map
    radianceMap = reconstructionPipeline(worldFrame, AGCSettings);

    % Demosaic the image
    radianceMap = demosaicRadianceMap(radianceMap);

    % Obtain the pixel indices within the world image for each check.
    rawPixelIndices = extractCheckerPixels(worldFrame*AGCSettings.Dgain,arducamB0392cameraIntrinsics.results.Intrinsics,cornerSets{valIdx},true);

    % Obtain the measured RGB radiance values
    for r = 1:nRows
        for c = 1:nColumns
            idx = rawPixelIndices{r, c};
            measuredRGBRadiance(r, c, 1) = mean(radianceMap(idx + (0*imgSize(1)*imgSize(2))), 'omitnan');
            measuredRGBRadiance(r, c, 2) = mean(radianceMap(idx + (1*imgSize(1)*imgSize(2))), 'omitnan');
            measuredRGBRadiance(r, c, 3) = mean(radianceMap(idx + (2*imgSize(1)*imgSize(2))), 'omitnan');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE CHANNEL FRACTION ADJUSTMENTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate total radiance per patch to extract fractional channel ratios
    predSum = sum(predictedRGBRadiance, 3);
    measSum = sum(measuredRGBRadiance, 3);

    predFracR = reshape(predictedRGBRadiance(:,:,1) ./ predSum, [], 1);
    predFracG = reshape(predictedRGBRadiance(:,:,2) ./ predSum, [], 1);
    predFracB = reshape(predictedRGBRadiance(:,:,3) ./ predSum, [], 1);

    measFracR = reshape(measuredRGBRadiance(:,:,1) ./ measSum, [], 1);
    measFracG = reshape(measuredRGBRadiance(:,:,2) ./ measSum, [], 1);
    measFracB = reshape(measuredRGBRadiance(:,:,3) ./ measSum, [], 1);

    % Find the median adjustment for this specific validation set
    adjR = median(predFracR ./ measFracR, 'omitnan');
    adjG = median(predFracG ./ measFracG, 'omitnan');
    adjB = median(predFracB ./ measFracB, 'omitnan');

    allAdjustments(valIdx, :) = [adjR, adjG, adjB];

end % loop over valSetOptions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATE AND SAVE RADIOMETRIC WEIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain the overall adjustment by taking the median across all validation sets
overallAdjustment = median(allAdjustments, 1, 'omitnan');

% Load the existing radiometric correction to update it
paramFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
load(paramFileName, 'radiometricCorrectionRGB');

% Apply the correction fraction to the existing weights
newRadiometricCorrectionRGB = radiometricCorrectionRGB .* overallAdjustment;

% Adjust this triplet so that the mean sensor value (across RGB) is
% unchanged by this operation. Need to account for the twice as numerous G
% pixels in the BGGR pattern.
k = 4 / (newRadiometricCorrectionRGB(1) + 2*newRadiometricCorrectionRGB(2) + newRadiometricCorrectionRGB(3));
radiometricCorrectionRGB = newRadiometricCorrectionRGB * k;

% Construct a map to apply this correction
radiometricCorrectionMap = ones(imgSize(1), imgSize(2));
[bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(radiometricCorrectionMap, bayerPattern);
for cc = 1:3
    radiometricCorrectionMap(bayerIdx{cc}) = radiometricCorrectionRGB(cc);
end

% Report the correction to the console
fprintf('\nThe absolute integrated radiance calibration scalar tuple (RGB) is: [%2.4f, %2.4f, %2.4f]\n\n', radiometricCorrectionRGB);

% Save the radiometric correction to the "derived" directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radiometricCorrectionRGB.mat');
readme = ['Created by defineRadiometricWeights (via Macbeth ColorChecker data).\n'...
    'radiometricCorrectionRGB -- multiply the (linearized) sensor values by these absolute calibration factors.\n'...
    'radiometricCorrectionMap -- a map of these absolute corrections that can be applied to an entire image.\n'];
save(saveFileName,'readme','radiometricCorrectionRGB','radiometricCorrectionMap');