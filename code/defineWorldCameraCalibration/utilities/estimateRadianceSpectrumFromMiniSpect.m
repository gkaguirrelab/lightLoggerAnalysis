function [spectralRadiance,S] = estimateRadianceSpectrumFromMinispect(miniSpectValues)
% Spectral Reconstruction from 9-Channel AS7341 Sensor using Tikhonov Regularization
% This function estimates the mean environmental radiance from the
% minispect values
%{
load("/Users/aguirre/Documents/MATLAB/projects/lightLoggerAnalysis/data/exampleWorldCameraImages/outdoor_AGCandMS_01.mat")
miniSpectValues = minispectValue.AS(1:9);
[spectralRadiance,S] = estimateRadianceSpectrumFromMiniSpect(miniSpectValues);
plot(SToWls(S),spectralRadiance)
%}

showPlots = true;

persistent miniSpectT miniSpectWls
if isempty(miniSpectT)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'ASM7341_spectralSensitivity.mat');
    load(paramFileName,'T');
    miniSpectWls = T.wl;
    miniSpectT = table2array(T(:,["F1" "F2" "F3","F4","F5","F6","F7","F8","Clear"]))';
    % Sensitivity is max-normalized, matching the calibration dot product
    miniSpectT = miniSpectT ./ max(miniSpectT')';
end

nChannels = size(miniSpectT,1);

persistent miniSpectKVals
if isempty(miniSpectKVals)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'minispectRadianceWeights.mat');
    load(paramFileName,"fitObj");
    miniSpectKVals = fitObj.coeff(1:nChannels,2);
end

N = length(miniSpectWls);
V = double(miniSpectValues');
k = double(miniSpectKVals);
A = miniSpectT; % Reverted: A is the max-normalized sensitivity matrix

%% 2. Construct the Forward Model
% Calculate calibrated integrated radiance (dot product) for each channel
y = V ./ 10.^k;

%% 3. Setup the Tikhonov Regularization Matrices
% Create the second-derivative matrix D to enforce smoothness
D = diff(eye(N), 2);

% Dynamic alpha scaling: Adjust regularization inversely to signal magnitude (SNR proxy)
baseAlpha = 1e5;
referenceIntensity = 0.05; % Replace with the empirical mean(y) from your typical calibration exposure

% Inverse scaling enforces heavier smoothing for dim (low SNR) scenes
alpha = baseAlpha * (referenceIntensity / (mean(y) + eps));

% Combine the data fidelity and smoothing terms for lsqlin
C = [A; sqrt(alpha) * D];
d = [y; zeros(N-2, 1)];

%% 4. Solve using Constrained Linear Least Squares
lb = zeros(N, 1);
ub = [];

options = optimoptions('lsqlin', 'Display', 'off');
x_est = lsqlin(C, d, [], [], [], [], lb, ub, [], options);

%% 5. Visualization
if showPlots
    figure;
    hold on; grid on;

    % Plot the estimated continuous spectrum (per-bin spectral radiance)
    plot(miniSpectWls, x_est, 'b-', 'LineWidth', 2, 'DisplayName', 'Reconstructed Spectrum');

    % Convert integrated radiance (y) to effective mean spectral radiance for visualization
    y_mean = y ./ sum(miniSpectT, 2);

    % Plot the effective mean radiances at the peak wavelength of each channel
    [~, max_idx] = max(miniSpectT, [], 2);
    plot(miniSpectWls(max_idx), y_mean, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Effective Mean Measurement');

    % Overlay normalized sensitivities lightly in the background for context, scaled to y_mean
    for i = 1:nChannels
        plot(miniSpectWls, miniSpectT(i,:) * max(y_mean), 'k--', 'Color', [0.7 0.7 0.7 0.5], 'HandleVisibility', 'off');
    end

    xlim([400 800]);
    xlabel('Wavelength (nm)');
    ylabel('Spectral Radiance');
    title('Environmental Spectrum Reconstruction');
    legend('Location', 'best');

    % Verification plot (comparing integrated values directly)
    figure
    y_predicted = A * x_est;
    bar(1:nChannels, [y, y_predicted], 'FaceColor', 'flat');
    grid on;
    xlabel('Sensor Channel (1-8 Narrow, 9 Clear)');
    ylabel('Integrated Radiance (Dot Product)');
    title('Verification: Measured vs. Predicted Integrated Radiance');
    legend('Measured (y)', 'Predicted (A*x_{est})', 'Location', 'northwest');
end

% Prepare variables for return
spectralRadiance = x_est;
S = WlsToS(miniSpectWls);

end