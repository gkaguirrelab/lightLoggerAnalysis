% calculatePatchSNR.m
% This script extracts the raw, pre-gain pixel counts for the 24 Macbeth 
% color patches across the six validation datasets. It calculates the true 
% photo-signal (raw - darkSignal) and the resulting Signal-to-Noise Ratio (SNR) 
% to evaluate the impact of the 8-bit quantization noise floor at high digital gains.

% Housekeeping
clear
close all

% Validation set definitions
valSetOptions = {'indoor','indoor','indoor','outdoor','outdoor','outdoor'};
valNumOptions = {'1','2','4','1','2','3'};
plotOrder = [4,5,1,2,6,3];

nRows = 4;
nColumns = 6;

% Corner sets from validation script
cornerSets{1} = [200.8987 153.3239; 406.0814 154.3870; 405.0183 284.0880; 197.7093 287.2774];
cornerSets{2} = [150.9319 166.0814; 457.1113 176.7126; 484.7525 386.1478; 99.9020 398.9053];
cornerSets{3} = [217.9086 160.7658; 449.6694 179.9020; 440.1013 335.1179; 197.7093 315.9817];
cornerSets{4} = [198.7724 194.7857; 418.8389 191.5963; 430.5332 339.3704; 186.0150 344.6860];
cornerSets{5} = [233.8555 234.1213; 430.5332 211.7957; 447.5432 342.5598; 242.3605 370.2010];
cornerSets{6} = [134.9850 128.8721; 468.8056 83.1578; 485.8156 287.2774; 184.9518 351.0648];

% Load persistent parameters
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'darkSignal.mat'), 'darkSignal');
load(fullfile(tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'arducamB0392cameraIntrinsics.mat'), 'arducamB0392cameraIntrinsics');

% Preallocate storage for plotting
allDgain = [];
allSignal = [];
allSNR = [];

% Loop over validation datasets
for vv = 1:length(plotOrder)
    valIdx = plotOrder(vv);
    
    % Load the raw camera frame and AGC settings
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'macbethColorCheck', ...
        valSetOptions{valIdx}, valNumOptions{valIdx}, 'lightLogger', 'close_AGCandMS_01.mat');
    load(dataFileName, 'worldFrame', 'AGCSettings');
    
    % Get pixel bounding boxes for the 24 patches
    rawPixelIndices = extractCheckerPixels(worldFrame * AGCSettings.Dgain, ...
        arducamB0392cameraIntrinsics.results.Intrinsics, cornerSets{valIdx});
    
    % Extract signal and noise per patch
    for r = 1:nRows
        for c = 1:nColumns
            idx = rawPixelIndices{r, c};
            
            % Extract the raw 8-bit counts for this specific patch
            rawPixels = double(worldFrame(idx));
            
            % Subtract baseline to isolate the true photo-signal
            photoSignal = rawPixels - darkSignal;
            
            % Calculate spatial mean and standard deviation (noise proxy)
            mu = mean(photoSignal, 'omitnan');
            sigma = std(photoSignal, 0, 'omitnan');
            
            % Prevent divide-by-zero for perfectly flat clamped patches
            if sigma == 0
                sigma = eps;
            end
            
            % Store metrics
            allDgain(end+1, 1) = AGCSettings.Dgain;
            allSignal(end+1, 1) = mu;
            allSNR(end+1, 1) = mu / sigma;
        end
    end
end

% =========================================================================
% PLOT RESULTS
% =========================================================================
figure('Name', 'Pre-Gain Signal and SNR Analysis', 'WindowStyle', 'Docked', ...
       'Position', [100, 100, 1000, 450]);
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "tight");

% Add small random horizontal jitter so the 24 patches don't perfectly overlap
jitter = (rand(size(allDgain)) - 0.5) * 0.2;

% --- TILE 1: Pre-Gain Photo-Signal ---
nexttile; hold on;
scatter(allDgain + jitter, allSignal, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
yline(0, 'k--', 'LineWidth', 1.5);
yline(5, 'r:', 'Critical Quantization Boundary (5 counts)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

xlabel('Digital Gain (Dgain)');
ylabel('Pre-Gain Photo-Signal (Raw Counts - Dark Signal)');
title('Absolute Signal Level per Patch');
grid on; box on;
xlim([0, max(allDgain) + 1]);

% --- TILE 2: Signal-to-Noise Ratio (SNR) ---
nexttile; hold on;
scatter(allDgain + jitter, allSNR, 40, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
yline(1, 'r--', 'Noise = Signal (SNR = 1)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

set(gca, 'YScale', 'log');
xlabel('Digital Gain (Dgain)');
ylabel('SNR (\mu / \sigma)');
title('Patch Signal-to-Noise Ratio');
grid on; box on;
xlim([0, max(allDgain) + 1]);

fprintf('SNR calculations complete. Displaying figure.\n');