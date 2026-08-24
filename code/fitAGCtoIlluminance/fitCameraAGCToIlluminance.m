function [predictedIlluminance, best_p, sensitivity, illuminance] = fitCameraAGCToIlluminance(cameraScore)
% Fit the empirical AGC-to-illuminance model and evaluate camera scores.
%
% This function is the shared implementation used by both
% fitEmpircalAGCtoIlluminance.m and Python's video_to_illuminance. It loads
% the selected camera-score / illuminance point cloud, performs the robust
% continuous two-slope fit in log10 space, and evaluates the fitted model for
% the supplied camera scores. Nonfinite and nonpositive inputs return NaN.

if nargin < 1
    cameraScore = [];
end

% Load and clean the selected environmental calibration points.
analysisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(analysisDir));
calibrationData = load(fullfile(projectRoot, "data", "empircalAGC.mat"));
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

% Evaluate only valid camera scores and preserve NaN for invalid inputs.
predictedIlluminance = nan(size(cameraScore));
validCameraScore = isfinite(cameraScore) & (cameraScore > 0);
predictedIlluminance(validCameraScore) = 10.^(piecewise_model(best_p, log10(cameraScore(validCameraScore))));

end
