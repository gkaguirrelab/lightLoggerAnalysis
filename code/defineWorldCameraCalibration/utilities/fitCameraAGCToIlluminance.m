function [predictedIlluminance, best_p, sensitivity, illuminance] = fitCameraAGCToIlluminance(cameraScore)
% Fit, save, and evaluate the empirical AGC-to-illuminance model.
%
% The fit is saved as cameraAGCToIlluminanceFit in
% derived/cameraAGCToIlluminanceFit.mat. Python loads that artifact directly
% rather than starting MATLAB and repeating the fit during video processing.
% Nonfinite and nonpositive cameraScore inputs return NaN.

if nargin < 1
    cameraScore = [];
end

% Locate the project from this file's utilities/define/code nesting.
analysisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(analysisDir)));
sourceDataPath = fullfile( ...
    projectRoot,"data","empircalAGCAndIlluminance.mat");

% Load and clean the selected environmental calibration points.
calibrationData = load(sourceDataPath,"empiralAGC");
x_data = calibrationData.empiralAGC.cameraScoreLinear;
y_data = calibrationData.empiralAGC.msIlluminance;
valid_idx = isfinite(x_data) & isfinite(y_data) & (x_data > 0) & (y_data > 0);
sensitivity = x_data(valid_idx);
illuminance = y_data(valid_idx);

% Transform the calibration data into log10 space.
log_x = log10(sensitivity);
log_y = log10(illuminance);

% Define the continuous two-slope model in log10 space.
piecewise_model = @(p, x) (x < p(4)) .* (p(1) .* x + p(2)) + (x >= p(4)) .* (p(1) * p(4) + p(2) + p(3) .* (x - p(4)));

% Minimize the L1 residual so extreme observations do not dominate the fit.
robust_func = @(p) sum(abs(log_y - piecewise_model(p, log_x)));
mid_x = mean(log_x);
initial_guess = [-0.5, 5, -1.5, mid_x];
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
best_p = fminsearch(robust_func, initial_guess, options);

% Save both human-readable named coefficients and the original parameter
% vector so other languages can consume a stable calibration contract.
cameraAGCToIlluminanceFit = struct();
cameraAGCToIlluminanceFit.parameterVector = best_p(:)';
cameraAGCToIlluminanceFit.slopeBelowBreakpoint = best_p(1);
cameraAGCToIlluminanceFit.interceptBelowBreakpoint = best_p(2);
cameraAGCToIlluminanceFit.slopeAboveBreakpoint = best_p(3);
cameraAGCToIlluminanceFit.log10CameraScoreBreakpoint = best_p(4);
cameraAGCToIlluminanceFit.sampleCount = numel(sensitivity);
cameraAGCToIlluminanceFit.cameraScoreRange = ...
    [min(sensitivity),max(sensitivity)];
cameraAGCToIlluminanceFit.illuminanceRange = ...
    [min(illuminance),max(illuminance)];
cameraAGCToIlluminanceFit.sourceDataFile = ...
    'data/empircalAGCAndIlluminance.mat';
cameraAGCToIlluminanceFit.model = ...
    'continuous two-slope piecewise linear model in log10 space';
cameraAGCToIlluminanceFit.objective = ...
    'minimum sum of absolute log10 illuminance residuals';
if isfield(calibrationData.empiralAGC,"sharedLagSeconds")
    cameraAGCToIlluminanceFit.sharedLagSeconds = ...
        calibrationData.empiralAGC.sharedLagSeconds;
end

fitOutputPath = fullfile( ...
    projectRoot,"derived","cameraAGCToIlluminanceFit.mat");
save(fitOutputPath,"cameraAGCToIlluminanceFit");
fprintf('Saved AGC-to-illuminance fit to %s\n',fitOutputPath);

% Evaluate only valid camera scores and preserve NaN for invalid inputs.
predictedIlluminance = nan(size(cameraScore));
validCameraScore = isfinite(cameraScore) & (cameraScore > 0);
predictedIlluminance(validCameraScore) = 10.^(piecewise_model(best_p, log10(cameraScore(validCameraScore))));

end
