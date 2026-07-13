%% debug_darkApertureDifference.m
% Visualize left-right aperture differences during the dark recording.

clear; close all; clc

darkDir = '/Users/sophiamirabal/Downloads/migraine_timeseriesData/m_dark_data/flic_1029_dark_2-69b12aae';

eye = readtable(fullfile(darkDir,'3d_eye_states.csv'), ...
    'VariableNamingRule','preserve');

t = (eye.("timestamp [ns]") - eye.("timestamp [ns]")(1))*1e-9;

apLeft  = eye.("eyelid aperture left [mm]");
apRight = eye.("eyelid aperture right [mm]");

%% Smooth a little so single-frame noise doesn't dominate
smoothWindow = 15;

apLeft_s  = smoothdata(apLeft,'movmedian',smoothWindow,'omitnan');
apRight_s = smoothdata(apRight,'movmedian',smoothWindow,'omitnan');

diffSigned = apLeft_s - apRight_s;
diffAbs    = abs(diffSigned);

%% Plot 1: left and right aperture

figure;
plot(t,apLeft_s,'LineWidth',1.5); hold on
plot(t,apRight_s,'LineWidth',1.5);

xlabel('Time (s)');
ylabel('Eye openness (mm)');
title('Dark Recording: Left and Right Eye Openness');
legend('Left','Right');
grid on

%% Plot 2: signed difference

figure;
plot(t,diffSigned,'k','LineWidth',1.5);
yline(0,'--');

xlabel('Time (s)');
ylabel('Left - Right (mm)');
title('Signed Eye Openness Difference');
grid on

%% Plot 3: absolute difference

figure;
plot(t,diffAbs,'r','LineWidth',1.5);

xlabel('Time (s)');
ylabel('|Left - Right| (mm)');
title('Absolute Eye Openness Difference');
grid on

%% Report largest disagreements

[sortedDiff,idx] = sort(diffAbs,'descend');

disp(table( ...
    t(idx(1:20)), ...
    apLeft(idx(1:20)), ...
    apRight(idx(1:20)), ...
    diffSigned(idx(1:20)), ...
    sortedDiff(1:20), ...
    'VariableNames', ...
    {'Time_s','Left_mm','Right_mm','SignedDiff_mm','AbsDiff_mm'}));


threshold = 2;   % mm

idxBad = diffAbs > threshold;

figure;
plot(t, apLeft_s, 'b'); hold on
plot(t, apRight_s, 'r');

scatter(t(idxBad), apLeft_s(idxBad), 25, 'b', 'filled');
scatter(t(idxBad), apRight_s(idxBad), 25, 'r', 'filled');

xlabel('Time (s)');
ylabel('Eye openness (mm)');
title('Frames where left/right aperture differs by >2 mm');
legend('Left','Right');
grid on

%% Average aperture difference over the dark recording

meanSignedDiff = mean(diffSigned,'omitnan');
stdSignedDiff  = std(diffSigned,'omitnan');

meanAbsDiff = mean(diffAbs,'omitnan');
stdAbsDiff  = std(diffAbs,'omitnan');

fprintf('\nAverage eye-openness difference (dark recording)\n');
fprintf('-----------------------------------------------\n');
fprintf('Mean (Left - Right): %.3f mm\n', meanSignedDiff);
fprintf('SD   (Left - Right): %.3f mm\n', stdSignedDiff);
fprintf('Mean |Left - Right|: %.3f mm\n', meanAbsDiff);
fprintf('SD   |Left - Right|: %.3f mm\n', stdAbsDiff);

%% Distribution of aperture differences

figure;
histogram(diffSigned,40);

xlabel('Left - Right eye openness (mm)');
ylabel('Count');
title(sprintf(['Dark Recording\n' ...
    'Mean = %.2f mm,  Mean |diff| = %.2f mm'], ...
    meanSignedDiff, meanAbsDiff));

grid on;
xline(meanSignedDiff,'r','LineWidth',2);

%% Mean eye openness over the dark recording

meanLeft = mean(apLeft,'omitnan');
meanRight = mean(apRight,'omitnan');

fprintf('\nMean eye openness (dark recording)\n');
fprintf('---------------------------------\n');
fprintf('Left eye : %.3f mm\n', meanLeft);
fprintf('Right eye: %.3f mm\n', meanRight);

figure;

bar([meanLeft, meanRight]);

xticklabels({'Left eye','Right eye'});
ylabel('Mean eye openness (mm)');
title('Average Eye Openness During Dark Recording');
grid on

%% Show exact dark-recording timepoints used to normalize the Read task
% These reproduce the current baseline-selection rules in
% plotParticipantState:
%   1. Valid left and right aperture measurements
%   2. Horizontal gaze within +/- 10 degrees
%   3. Left-right aperture difference <= 2 mm
%   4. Top 5% of joint left/right openness
%
% The resulting dark baselines are the values used to normalize the
% separate left and right eye-aperture signals in the Read recording.

%% Load dark gaze data
darkGaze = readtable(fullfile(darkDir, 'gaze.csv'), ...
    'VariableNamingRule', 'preserve');

eyeTime_ns = double(eye.("timestamp [ns]"));
gazeTime_ns = double(darkGaze.("timestamp [ns]"));

darkAzimuth = darkGaze.("azimuth [deg]");

% Align gaze samples to the eye-state timestamps
azimuthAtEye = interp1( ...
    gazeTime_ns, ...
    darkAzimuth, ...
    eyeTime_ns, ...
    'nearest', ...
    NaN);

%% Baseline-selection settings
baselinePercentile = 95;
maxAzimuth_deg = 10;
maxEyeDifference_mm = 2;

%% First-stage filtering
validDark = ...
    isfinite(apLeft) & ...
    isfinite(apRight) & ...
    isfinite(azimuthAtEye) & ...
    abs(azimuthAtEye) <= maxAzimuth_deg;

validIdx = find(validDark);

apLeft_valid = apLeft(validDark);
apRight_valid = apRight(validDark);
azimuth_valid = azimuthAtEye(validDark);

if isempty(apLeft_valid)
    error(['No dark-recording samples survived the valid-data and ', ...
           'horizontal-gaze filters.']);
end

%% Joint openness ranking
% Rank each eye separately so the two signals remain separate but are
% placed on comparable percentile scales.
rankLeft = tiedrank(apLeft_valid) ./ numel(apLeft_valid);
rankRight = tiedrank(apRight_valid) ./ numel(apRight_valid);

% A sample ranks highly only when both eyes are highly open.
jointOpennessRank = min(rankLeft, rankRight);

%% Reject large left-right differences
eyeDifference_valid = abs(apLeft_valid - apRight_valid);

%% Select top jointly open samples
jointThreshold = prctile( ...
    jointOpennessRank, ...
    baselinePercentile);

selectedValid = ...
    jointOpennessRank >= jointThreshold & ...
    eyeDifference_valid <= maxEyeDifference_mm;

% Return the selected samples to full-recording indexing
selectedMask = false(size(apLeft));
selectedMask(validIdx(selectedValid)) = true;

selectedTimes_s = t(selectedMask);

%% Calculate the actual separate baselines
darkBaselineLeft_mm = mean(apLeft(selectedMask), 'omitnan');
darkBaselineRight_mm = mean(apRight(selectedMask), 'omitnan');

fprintf('\nDark normalization samples used for the Read task\n');
fprintf('-------------------------------------------------\n');
fprintf('Number of selected samples: %d\n', sum(selectedMask));
fprintf('Left dark baseline:  %.3f mm\n', darkBaselineLeft_mm);
fprintf('Right dark baseline: %.3f mm\n', darkBaselineRight_mm);

%% Plot aperture and selected normalization timepoints
figure;

plot(t, apLeft_s, ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Left eye');
xlim([0 150]);
hold on

plot(t, apRight_s, ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Right eye');
xlim([0 150]);

% Draw one vertical line at every selected dark-recording timestamp
for k = 1:numel(selectedTimes_s)
    xline(selectedTimes_s(k), ...
        ':', ...
        'Color', [0.4 0.4 0.4], ...
        'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% Dummy object so the vertical lines appear once in the legend
pSelected = plot(nan, nan, ...
    ':', ...
    'Color', [0.4 0.4 0.4], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Selected normalization timepoint');

xlabel('Time from dark recording start (s)');
ylabel('Eye openness (mm)');
title('Dark Timepoints Used to Normalize for FLIC_1029');

legend('Location', 'best');
grid on

% Print the selected timestamps for normalization
disp('Selected normalization timepoints (s):');
disp(selectedTimes_s);
