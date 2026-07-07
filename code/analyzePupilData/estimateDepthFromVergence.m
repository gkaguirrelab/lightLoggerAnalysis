function [depthTable, summaryStats] = estimateDepthFromVergence(recordingDir, calibVergenceDeg, taskName)
% estimateDepthFromVergence
%
% Estimate relative depth over time from 3D optical-axis separation.
%
% Inputs:
%   recordingDir        folder containing 3d_eye_states.csv
%   calibVergenceDeg    median 3D optical-axis separation when target = 1 m
%   taskName            name used in plot title
%
% Outputs:
%   depthTable          table with time, vergence, and estimated depth
%   summaryStats        table with summary depth statistics
%
% Calibration assumption:
%   depth_m = calibVergenceDeg ./ currentVergenceDeg
%
% Because calibVergenceDeg corresponds to 1 meter, values larger than
% calibVergenceDeg estimate closer than 1 m, and values smaller than
% calibVergenceDeg estimate farther than 1 m.

%% Defaults
if nargin < 3
    taskName = 'Recording';
end

smoothWindow = 15;

%% Load eye-state data
eyeFile = fullfile(recordingDir, '3d_eye_states.csv');

eye = readtable(eyeFile, 'VariableNamingRule', 'preserve');

t = (eye.("timestamp [ns]") - eye.("timestamp [ns]")(1)) * 1e-9;

%% Optical-axis vectors
vLeft = [ ...
    eye.("optical axis left x"), ...
    eye.("optical axis left y"), ...
    eye.("optical axis left z")];

vRight = [ ...
    eye.("optical axis right x"), ...
    eye.("optical axis right y"), ...
    eye.("optical axis right z")];

%% 3D angular separation between left/right optical axes
dotLR = sum(vLeft .* vRight, 2);
normL = vecnorm(vLeft, 2, 2);
normR = vecnorm(vRight, 2, 2);

cosTheta = dotLR ./ (normL .* normR);
cosTheta = max(min(cosTheta, 1), -1);

vergenceDeg = acosd(cosTheta);

%% Convert vergence to estimated depth
estimatedDepth_m = calibVergenceDeg ./ vergenceDeg;

% Remove impossible/extreme values for visualization
estimatedDepth_m(vergenceDeg <= 0) = NaN;
estimatedDepth_m(estimatedDepth_m > 10) = NaN;

%% Smooth for plotting
vergenceDeg_s = smoothdata(vergenceDeg, ...
    'movmedian', smoothWindow, 'omitnan');

estimatedDepth_s = smoothdata(estimatedDepth_m, ...
    'movmedian', smoothWindow, 'omitnan');

%% Plot estimated depth over time
figure;
plot(t, estimatedDepth_s, 'LineWidth', 1.2);
yline(1, '--', '1 m calibration depth');

xlabel('Time from recording start (s)');
ylabel('Estimated depth (m)');
title(['Estimated Depth Over Time - ' taskName]);
grid on

%% Optional plot: vergence over time
figure;
plot(t, vergenceDeg_s, 'LineWidth', 1.2);
yline(calibVergenceDeg, '--', '1 m calibration vergence');

xlabel('Time from recording start (s)');
ylabel('3D optical-axis separation (deg)');
title(['Vergence Signal Over Time - ' taskName]);
grid on

%% Output tables
depthTable = table( ...
    t, ...
    vergenceDeg, ...
    estimatedDepth_m, ...
    'VariableNames', {'Time_s', ...
                      'Vergence3D_deg', ...
                      'EstimatedDepth_m'} ...
);

summaryStats = table( ...
    median(estimatedDepth_m, 'omitnan'), ...
    mean(estimatedDepth_m, 'omitnan'), ...
    std(estimatedDepth_m, 'omitnan'), ...
    prctile(estimatedDepth_m, 25), ...
    prctile(estimatedDepth_m, 75), ...
    sum(~isnan(estimatedDepth_m)), ...
    'VariableNames', {'MedianDepth_m', ...
                      'MeanDepth_m', ...
                      'SDDepth_m', ...
                      'P25Depth_m', ...
                      'P75Depth_m', ...
                      'NSamples'} ...
);

disp(summaryStats)

end