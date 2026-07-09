%% visualizePupilCloudData.m
% Visualize Pupil Cloud eye-state, gaze, blink, and saccade outputs.

clear; close all; clc

%% File names
eyeFile      = '3d_eye_states.csv';
gazeFile     = 'gaze.csv';
blinksFile   = 'blinks.csv';
saccadesFile = 'saccades.csv';

%% Load data
eye      = readtable(eyeFile, 'VariableNamingRule', 'preserve');
gaze     = readtable(gazeFile, 'VariableNamingRule', 'preserve');
blinks   = readtable(blinksFile, 'VariableNamingRule', 'preserve');
saccades = readtable(saccadesFile, 'VariableNamingRule', 'preserve');

%% Time vectors, in seconds from recording start
t0 = eye.("timestamp [ns]")(1);

tEye = (eye.("timestamp [ns]") - t0) * 1e-9;
tGaze = (gaze.("timestamp [ns]") - t0) * 1e-9;

blinkStart = (blinks.("start timestamp [ns]") - t0) * 1e-9;
blinkEnd   = (blinks.("end timestamp [ns]") - t0) * 1e-9;

saccadeStart = (saccades.("start timestamp [ns]") - t0) * 1e-9;
saccadeEnd   = (saccades.("end timestamp [ns]") - t0) * 1e-9;

%% Extract eye-state variables
apertureL = eye.("eyelid aperture left [mm]");
apertureR = eye.("eyelid aperture right [mm]");

pupilL = eye.("pupil diameter left [mm]");
pupilR = eye.("pupil diameter right [mm]");

apertureDiff = apertureL - apertureR;
pupilDiff = pupilL - pupilR;

%% Optional smoothing for visualization
smoothWindow = 15;

apertureL_s = smoothdata(apertureL, 'movmedian', smoothWindow, 'omitnan');
apertureR_s = smoothdata(apertureR, 'movmedian', smoothWindow, 'omitnan');

pupilL_s = smoothdata(pupilL, 'movmedian', smoothWindow, 'omitnan');
pupilR_s = smoothdata(pupilR, 'movmedian', smoothWindow, 'omitnan');

apertureDiff_s = apertureL_s - apertureR_s;
pupilDiff_s = pupilL_s - pupilR_s;

%% 1. Eye openness
figure;
plot(tEye, apertureL_s, 'LineWidth', 1.2); hold on
plot(tEye, apertureR_s, 'LineWidth', 1.2);
xlabel('Time from recording start (s)');
ylabel('Eyelid aperture (mm)');
title('Eye Openness: Left vs Right Eyelid Aperture');
legend({'Left eye', 'Right eye'}, 'Location', 'best');
grid on

%% 2. Pupil size
figure;
plot(tEye, pupilL_s, 'LineWidth', 1.2); hold on
plot(tEye, pupilR_s, 'LineWidth', 1.2);
xlabel('Time from recording start (s)');
ylabel('Pupil diameter (mm)');
title('Pupil Size: Left vs Right Pupil Diameter');
legend({'Left eye', 'Right eye'}, 'Location', 'best');
grid on

%% 3a. Eye openness asymmetry
figure;
plot(tEye, apertureDiff_s, 'LineWidth', 1.2); hold on
yline(0, '--');
xlabel('Time from recording start (s)');
ylabel('Left - Right aperture (mm)');
title('Eye Openness Asymmetry');
grid on

%% 3b. Pupil size asymmetry
figure;
plot(tEye, pupilDiff_s, 'LineWidth', 1.2); hold on
yline(0, '--');
xlabel('Time from recording start (s)');
ylabel('Left - Right pupil diameter (mm)');
title('Pupil Diameter Asymmetry');
grid on

%% 4. Blink periods overlaid on aperture
figure;
plot(tEye, apertureL_s, 'LineWidth', 1.2); hold on
plot(tEye, apertureR_s, 'LineWidth', 1.2);

yl = ylim;

for i = 1:height(blinks)
    patch([blinkStart(i) blinkEnd(i) blinkEnd(i) blinkStart(i)], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.8 0.8 0.8], ...
          'FaceAlpha', 0.25, ...
          'EdgeColor', 'none');
end

% Put aperture lines back on top of blink shading
plot(tEye, apertureL_s, 'LineWidth', 1.2);
plot(tEye, apertureR_s, 'LineWidth', 1.2);

xlabel('Time from recording start (s)');
ylabel('Eyelid aperture (mm)');
title('Blink Periods Overlaid on Eyelid Aperture');
legend({'Left aperture', 'Right aperture', 'Blink period'}, 'Location', 'best');
grid on

%% 5. Gaze: combined vs monocular left/right
gazeX = gaze.("gaze x [px]");
gazeY = gaze.("gaze y [px]");

gazeLeftX = gaze.("gaze mono left x [px]");
gazeLeftY = gaze.("gaze mono left y [px]");

gazeRightX = gaze.("gaze mono right x [px]");
gazeRightY = gaze.("gaze mono right y [px]");

figure;
plot(tGaze, gazeX, 'LineWidth', 1.1); hold on
plot(tGaze, gazeLeftX, 'LineWidth', 1.1);
plot(tGaze, gazeRightX, 'LineWidth', 1.1);
xlabel('Time from recording start (s)');
ylabel('Gaze x position (px)');
title('Horizontal Gaze: Combined vs Left/Right Monocular');
legend({'Combined gaze', 'Left eye gaze', 'Right eye gaze'}, 'Location', 'best');
grid on

figure;
plot(tGaze, gazeY, 'LineWidth', 1.1); hold on
plot(tGaze, gazeLeftY, 'LineWidth', 1.1);
plot(tGaze, gazeRightY, 'LineWidth', 1.1);
xlabel('Time from recording start (s)');
ylabel('Gaze y position (px)');
title('Vertical Gaze: Combined vs Left/Right Monocular');
legend({'Combined gaze', 'Left eye gaze', 'Right eye gaze'}, 'Location', 'best');
grid on

%% 6. Gaze disagreement between eyes
gazeDisagreement = sqrt((gazeLeftX - gazeRightX).^2 + ...
                        (gazeLeftY - gazeRightY).^2);

figure;
plot(tGaze, gazeDisagreement, 'LineWidth', 1.2);
xlabel('Time from recording start (s)');
ylabel('Left-right gaze disagreement (px)');
title('Difference Between Left- and Right-Eye Gaze Estimates');
grid on

%% 7. Saccades overlaid on combined gaze x
figure;
plot(tGaze, gazeX, 'LineWidth', 1.1); hold on

yl = ylim;

for i = 1:height(saccades)
    patch([saccadeStart(i) saccadeEnd(i) saccadeEnd(i) saccadeStart(i)], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.8 0.8 0.8], ...
          'FaceAlpha', 0.25, ...
          'EdgeColor', 'none');
end

plot(tGaze, gazeX, 'LineWidth', 1.1);

xlabel('Time from recording start (s)');
ylabel('Gaze x position (px)');
title('Saccade Periods Overlaid on Combined Horizontal Gaze');
legend({'Combined gaze x', 'Saccade period'}, 'Location', 'best');
grid on

%% Optical-axis angles: horizontal, vertical, and vergence
% Pupil Labs coordinate system:
% X = horizontal
% Y = vertical
% Z = forward/depth
%
% The optical-axis columns are 3D direction vectors for each eye.
% We convert those vectors into horizontal and vertical angles.

axisLX = eye.("optical axis left x");
axisLY = eye.("optical axis left y");
axisLZ = eye.("optical axis left z");

axisRX = eye.("optical axis right x");
axisRY = eye.("optical axis right y");
axisRZ = eye.("optical axis right z");

% Horizontal angle: left-right rotation relative to forward Z axis
hAngleLeft_deg  = atan2d(axisLX, axisLZ);
hAngleRight_deg = atan2d(axisRX, axisRZ);

% Vertical angle: up-down rotation relative to forward Z axis
vAngleLeft_deg  = atan2d(axisLY, axisLZ);
vAngleRight_deg = atan2d(axisRY, axisRZ);

% HORIZONTAL VERTICAL DEG DIFFERENCE
horizontalAngleDiff_deg = hAngleLeft_deg - hAngleRight_deg;
verticalAngleDiff_deg = vAngleLeft_deg - vAngleRight_deg;

% horizontalAngleDiff_deg = abs(hAngleLeft_deg - hAngleRight_deg);
% verticalAngleDiff_deg   = abs(vAngleLeft_deg - vAngleRight_deg);

% Optional smoothing for visualization only
angleSmoothWindow = 15;

hAngleLeft_s  = smoothdata(hAngleLeft_deg,  'movmedian', angleSmoothWindow, 'omitnan');
hAngleRight_s = smoothdata(hAngleRight_deg, 'movmedian', angleSmoothWindow, 'omitnan');

vAngleLeft_s  = smoothdata(vAngleLeft_deg,  'movmedian', angleSmoothWindow, 'omitnan');
vAngleRight_s = smoothdata(vAngleRight_deg, 'movmedian', angleSmoothWindow, 'omitnan');

horizontalAngleDiff_s = smoothdata(horizontalAngleDiff_deg, ...
    'movmedian', angleSmoothWindow, 'omitnan');

verticalAngleDiff_s = smoothdata(verticalAngleDiff_deg, ...
    'movmedian', angleSmoothWindow, 'omitnan');

%% Plot: horizontal eye angles
figure;
plot(tEye, hAngleLeft_s, 'LineWidth', 1.2); hold on
plot(tEye, hAngleRight_s, 'LineWidth', 1.2);
xlabel('Time from recording start (s)');
ylabel('Horizontal angle (deg)');
title('Horizontal Eye Angle from Optical Axis - Read Task');
legend({'Left eye', 'Right eye'}, 'Location', 'best');
grid on

%% Plot: vertical eye angles
figure;
plot(tEye, vAngleLeft_s, 'LineWidth', 1.2); hold on
plot(tEye, vAngleRight_s, 'LineWidth', 1.2);
xlabel('Time from recording start (s)');
ylabel('Vertical angle (deg)');
title('Vertical Eye Angle from Optical Axis - Read Task');
legend({'Left eye', 'Right eye'}, 'Location', 'best');
grid on

%% Plot: horizontal vergence candidate and vertical difference
figure;
plot(tEye, horizontalAngleDiff_s, 'LineWidth', 1.2); hold on
plot(tEye, verticalAngleDiff_s, 'LineWidth', 1.2);
yline(0, '--');
xlabel('Time from recording start (s)');
ylabel('Left - Right angle difference (deg)');
title('Horizontal Angle Difference and Vertical Angle Difference - Read Task');
legend({'Horizontal difference / candidate vergence', ...
    'Vertical difference'}, ...
    'Location', 'best');
grid on

%% Plot: absolute horizontal and vertical angle differences
% figure;
% plot(tEye, horizontalAngleDiff_s, 'LineWidth', 1.2); hold on
% plot(tEye, verticalAngleDiff_s, 'LineWidth', 1.2);
% xlabel('Time from recording start (s)');
% ylabel('Absolute angle difference (deg)');
% title('Absolute Horizontal and Vertical Eye-Angle Differences');
% legend({'|Horizontal left-right difference| / candidate vergence magnitude', ...
%     '|Vertical left-right difference|'}, ...
%     'Location', 'best');
% grid on

%% Summary
fprintf('\nSummary:\n');
fprintf('Mean left eyelid aperture: %.3f mm\n', mean(apertureL, 'omitnan'));
fprintf('Mean right eyelid aperture: %.3f mm\n', mean(apertureR, 'omitnan'));
fprintf('Mean left pupil diameter: %.3f mm\n', mean(pupilL, 'omitnan'));
fprintf('Mean right pupil diameter: %.3f mm\n', mean(pupilR, 'omitnan'));
fprintf('Mean left-right aperture difference: %.3f mm\n', mean(apertureDiff, 'omitnan'));
fprintf('Mean left-right pupil difference: %.3f mm\n', mean(pupilDiff, 'omitnan'));
fprintf('Mean left-right gaze disagreement: %.3f px\n', mean(gazeDisagreement, 'omitnan'));

fprintf('\nOptical Axis Angle Difference Summary:\n');
fprintf('Mean absolute horizontal difference: %.3f deg\n', ...
    mean(horizontalAngleDiff_deg, 'omitnan'));
fprintf('SD absolute horizontal difference: %.3f deg\n', ...
    std(horizontalAngleDiff_deg, 'omitnan'));
fprintf('Mean absolute vertical difference: %.3f deg\n', ...
    mean(verticalAngleDiff_deg, 'omitnan'));
fprintf('SD absolute vertical difference: %.3f deg\n', ...
    std(verticalAngleDiff_deg, 'omitnan'));


%% Compare signed horizontal angle difference across two recordings

walkInDir = '/Users/sophiamirabal/Downloads/PupilData_Analysis/walkIndoor_timeseriesData/flic_1034_walkindoor_1-86ecaaf3';
readDir = '/Users/sophiamirabal/Downloads/PupilData_Analysis/read_timeseriesData/flic_1034_read_1-a5ad020a';

walkInEye = readtable(fullfile(walkInDir, '3d_eye_states.csv'), ...
    'VariableNamingRule', 'preserve');

readEye = readtable(fullfile(readDir, '3d_eye_states.csv'), ...
    'VariableNamingRule', 'preserve');

tWalkIn = (walkInEye.("timestamp [ns]") - walkInEye.("timestamp [ns]")(1)) * 1e-9;
tRead = (readEye.("timestamp [ns]") - readEye.("timestamp [ns]")(1)) * 1e-9;

hLeftWalkIn = atan2d(walkInEye.("optical axis left x"), ...
    walkInEye.("optical axis left z"));
hRightWalkIn = atan2d(walkInEye.("optical axis right x"), ...
    walkInEye.("optical axis right z"));

horizDiffWalkIn = hLeftWalkIn - hRightWalkIn;

hLeftRead = atan2d(readEye.("optical axis left x"), ...
    readEye.("optical axis left z"));
hRightRead = atan2d(readEye.("optical axis right x"), ...
    readEye.("optical axis right z"));

horizDiffRead = hLeftRead - hRightRead;

horizDiffWalkIn_s = smoothdata(horizDiffWalkIn, ...
    'movmedian', angleSmoothWindow, 'omitnan');

horizDiffRead_s = smoothdata(horizDiffRead, ...
    'movmedian', angleSmoothWindow, 'omitnan');

figure;
plot(tWalkIn, horizDiffWalkIn_s, 'LineWidth', 1.2); hold on
plot(tRead, horizDiffRead_s, 'LineWidth', 1.2);
yline(0, '--');

xlabel('Time from recording start (s)');
ylabel('Horizontal angle difference, Left - Right (deg)');
title('Signed Horizontal Angle Difference Across Recordings');
legend({'Walk indoor task', 'Read task'}, 'Location', 'best');
grid on

fprintf('\nRecording Comparison: Signed Horizontal Angle Difference\n');
fprintf('Walk indoor task mean ± SD: %.3f ± %.3f deg\n', ...
    mean(horizDiffWalkIn, 'omitnan'), std(horizDiffWalkIn, 'omitnan'));

fprintf('Read task mean ± SD: %.3f ± %.3f deg\n', ...
    mean(horizDiffRead, 'omitnan'), std(horizDiffRead, 'omitnan'));