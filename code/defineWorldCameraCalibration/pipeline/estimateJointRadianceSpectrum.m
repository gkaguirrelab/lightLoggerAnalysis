function [spectralRadiance, S, fVal, fitErrors] = estimateJointRadianceSpectrum(minispectValues, cameraRGBRadiance)
% Spectral Reconstruction combining 7 of the AS7341 Channels and 3 IMX219
% Channels. Uses joint Tikhonov Regularization to estimate mean
% environmental radiance.

showPlots = false;

%% 1. Load and Format ASM7341 Sensitivities
persistent miniSpectT miniSpectWls miniSpectKVals
if isempty(miniSpectT)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'ASM7341_spectralSensitivity.mat');
    load(paramFileName,'T');
    miniSpectWls = T.wl;
    miniSpectT = table2array(T(:,["F1" "F2" "F3","F4","F5","F6","F7"]))';
    % Sensitivity is max-normalized
    miniSpectT = miniSpectT ./ max(miniSpectT')';
end

if isempty(miniSpectKVals)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'), 'derived', 'minispectRadianceWeights.mat');
    load(paramFileName,"fitObj");
    nChannelsAS = size(miniSpectT,1);
    miniSpectKVals = fitObj.coeff(1:nChannelsAS,2);
end

%% 2. Load and Format IMX219 Sensitivities
persistent cameraT
if isempty(cameraT)
    paramFileNameCam = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'), 'data', 'IMX219_spectralSensitivity.mat');
    camData = load(paramFileNameCam,'T');
    camWls = camData.T{:, 1};
    camSens = [camData.T.red, camData.T.green, camData.T.blue];

    % Spline to ASM7341 wavelength support
    camSensSplined = SplineRaw(camWls, camSens, WlsToS(miniSpectWls));

    % Max-normalize to match the camera calibration integration logic
    cameraT = (camSensSplined ./ max(camSensSplined, [], 1))'; 
end


%% 3. Construct the Unified Forward Model
N = length(miniSpectWls);

% AS7341 calibrated integrated radiance
V = double(minispectValues(1:size(miniSpectT,1))');
k = double(miniSpectKVals);
y_as = V ./ 10.^k;

% IMX219 calibrated integrated radiance (ensure column vector)
y_cam = double(cameraRGBRadiance(:));

% Concatenate modalities
A_joint = [miniSpectT; cameraT];
y_joint = [y_as; y_cam];


%% 4. Setup Joint Tikhonov Regularization (with Shrinkage)
% Create the second-derivative matrix D to enforce smoothness
D = diff(eye(N), 2);

% Create the identity matrix I_mat to enforce shrinkage (minimize total radiance)
I_mat = eye(N);

% Dynamic alpha scaling based on joint signal magnitude
baseAlpha = 1e7;
referenceIntensity = 0.05; 
intensityScale = (referenceIntensity / (mean(y_joint) + eps));

% Alpha (smoothness) scales INVERSELY to enforce more smoothing in dim/noisy light
alpha = baseAlpha * intensityScale;

% Beta (shrinkage) remains CONSTANT to provide a penalty strictly proportional 
% to the signal magnitude across all light levels. 
% (Note: You will likely need to tune this down; 1e2 to 1e4 is a good starting range).
beta = 1e1; 

% Combine the data fidelity, smoothing, and shrinkage terms
C = [A_joint; sqrt(alpha) * D; sqrt(beta) * I_mat];
d = [y_joint; zeros(N-2, 1); zeros(N, 1)];


%% 5. Solve using Constrained Linear Least Squares
lb = zeros(N, 1);

% Initialize upper bounds to infinity
ub = inf(N, 1);

% Constrain the spectrum to zero power at or below 300 nm and at or above 800 nm.
% Because lb is 0, setting ub to 0 strictly forces the solution to 0 at these indices.
zeroPowerIdx = (miniSpectWls <= 200) | (miniSpectWls >= 1000);
ub(zeroPowerIdx) = 0;

options = optimoptions('lsqlin', 'Display', 'off');
x_est = lsqlin(C, d, [], [], [], [], lb, ub, [], options);

y_predicted = A_joint * x_est;
fitErrors = y_joint - y_predicted;
fVal = norm(fitErrors);


%% 6. Visualization
if showPlots
    % --- DIAGNOSTIC PLOT ---
    figure('Name', 'Estimated Spectrum and Scaled Sensitivities');
    hold on; grid on;
    
    % Plot the estimated continuous mean radiance spectrum
    plot(miniSpectWls, x_est, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Estimated Spectrum');
    
    numChannels = length(y_predicted);
    
    for i = 1:numChannels
        % Identify lambda max for color mapping
        sensitivitySpectrum = A_joint(i, :);
        [~, maxIdx] = max(sensitivitySpectrum);
        lambdaMax = miniSpectWls(maxIdx);
        
        % Obtain RGB color corresponding to lambda max
        lineColor = wavelength2color(lambdaMax);
        
        % Convert integrated radiance back to an effective mean spectral radiance
        % so it plots on the same scale as x_est
        effectiveMeanRadiance = y_predicted(i) / sum(sensitivitySpectrum);
        scaledSensitivity = sensitivitySpectrum * effectiveMeanRadiance;
        
        % Format line thickness and labels based on the sensor chip
        if i <= size(miniSpectT, 1)
            chName = sprintf('ASM7341 Ch %d (\\lambda_{max}=%d)', i, round(lambdaMax));
            lineWidth = 1.0; % Thin solid line for AS chip
        else
            rgbNames = {'IMX219 Red', 'IMX219 Green', 'IMX219 Blue'};
            chName = sprintf('%s (\\lambda_{max}=%d)', rgbNames{i - size(miniSpectT, 1)}, round(lambdaMax));
            lineWidth = 3.0; % Thick solid line for IMX chip
        end
        
        plot(miniSpectWls, scaledSensitivity, '-', 'Color', lineColor, ...
            'LineWidth', lineWidth, 'DisplayName', chName);
    end
    
    % Restrict the X axis to the relevant visible/NIR domain
    xlim([400 800]);
    
    xlabel('Wavelength (nm)');
    ylabel('Radiance (W/m^2/sr/nm)');
    title('Estimated Radiance Spectrum & Scaled Sensor Sensitivities');
    legend('Location', 'eastoutside');
    hold off;

    % --- EXISTING BAR PLOT ---
    figure('Name', 'Measured vs Predicted Integrated Radiance');
    bar(1:length(y_joint), [y_joint, y_predicted], 'FaceColor', 'flat');
    grid on;
    xlabel('Sensor Channel (1-8 ASM7341, 9-11 IMX219 RGB)');
    ylabel('Integrated Radiance (Dot Product)');
    title('Verification: Measured vs. Predicted Integrated Radiance');
    legend('Measured (y)', 'Predicted (A*x_{est})', 'Location', 'northwest');
end

spectralRadiance = x_est;
S = WlsToS(miniSpectWls);
end

%% LOCAL FUNCTIONS

function color = wavelength2color(lambda)
    % Converts a given wavelength (in nm) to an approximate RGB color triplet.
    if lambda >= 380 && lambda < 440
        R = -(lambda - 440) / (440 - 380);
        G = 0.0;
        B = 1.0;
    elseif lambda >= 440 && lambda < 490
        R = 0.0;
        G = (lambda - 440) / (490 - 440);
        B = 1.0;
    elseif lambda >= 490 && lambda < 510
        R = 0.0;
        G = 1.0;
        B = -(lambda - 510) / (510 - 490);
    elseif lambda >= 510 && lambda < 580
        R = (lambda - 510) / (580 - 510);
        G = 1.0;
        B = 0.0;
    elseif lambda >= 580 && lambda < 645
        R = 1.0;
        G = -(lambda - 645) / (645 - 580);
        B = 0.0;
    elseif lambda >= 645 && lambda <= 780
        R = 1.0;
        G = 0.0;
        B = 0.0;
    else
        R = 0.0; G = 0.0; B = 0.0; 
    end

    % Adjust intensity near the limits of human vision
    if lambda >= 380 && lambda < 420
        factor = 0.3 + 0.7*(lambda - 380) / (420 - 380);
    elseif lambda >= 420 && lambda < 701
        factor = 1.0;
    elseif lambda >= 701 && lambda <= 780
        factor = 0.3 + 0.7*(780 - lambda) / (780 - 700);
    else
        factor = 0.3;
    end

    color = [R, G, B] * factor;
    
    % Cap edge cases extending into UV or IR for plotting visibility
    if lambda > 780
       color = [0.5 0 0]; % Deep red
    elseif lambda < 380
       color = [0.3 0 0.5]; % Deep violet
    end
end