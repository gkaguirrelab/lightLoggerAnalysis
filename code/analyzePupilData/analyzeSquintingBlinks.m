%% analyzeSquintingBlinks.m
% Analyze left/right squinting from Neon 3d_eye_states.csv for one recording.

clear; close all; clc

%% Choose one recording
recName = 'Walk Indoor';
recDir = '/Users/sophiamirabal/Downloads/PupilData_Analysis/walkIndoor_timeseriesData/flic_1034_walkindoor_1-86ecaaf3';

eyeFile = fullfile(recDir, '3d_eye_states.csv');
blinkFile = fullfile(recDir, 'blinks.csv');

%% Load data
eye = readtable(eyeFile, 'VariableNamingRule', 'preserve');

t0 = eye.("timestamp [ns]")(1);
t = (eye.("timestamp [ns]") - t0) * 1e-9;

apertureL = eye.("eyelid aperture left [mm]");
apertureR = eye.("eyelid aperture right [mm]");

pupilL = eye.("pupil diameter left [mm]");
pupilR = eye.("pupil diameter right [mm]");

%% Optional: load blink periods
hasBlinkFile = isfile(blinkFile);

isBlink = false(size(t));

if hasBlinkFile
    blinks = readtable(blinkFile, 'VariableNamingRule', 'preserve');

    blinkStart = (blinks.("start timestamp [ns]") - t0) * 1e-9;
    blinkEnd   = (blinks.("end timestamp [ns]") - t0) * 1e-9;

    for b = 1:height(blinks)
        isBlink = isBlink | (t >= blinkStart(b) & t <= blinkEnd(b));
    end
end

%% Remove blinks from squint baseline/statistics
apertureL_noBlink = apertureL;
apertureR_noBlink = apertureR;

apertureL_noBlink(isBlink) = NaN;
apertureR_noBlink(isBlink) = NaN;

%% Smooth aperture for visualization
smoothWindow = 15;

apertureL_s = smoothdata(apertureL, 'movmedian', smoothWindow, 'omitnan');
apertureR_s = smoothdata(apertureR, 'movmedian', smoothWindow, 'omitnan');

%% Define open-eye baseline
% 90th percentile = typical "open" aperture, excluding blinks.
baselineL = prctile(apertureL_noBlink, 90);
baselineR = prctile(apertureR_noBlink, 90);

%% Squint index
% 0 = normal/open
% larger value = more closed/squinted
squintL = 1 - apertureL ./ baselineL;
squintR = 1 - apertureR ./ baselineR;

squintAsymmetry = squintL - squintR;

squintL(isBlink) = NaN;
squintR(isBlink) = NaN;

squintMean = mean([squintL squintR], 2, 'omitnan');

squintL_s = smoothdata(squintL, 'movmedian', smoothWindow, 'omitnan');
squintR_s = smoothdata(squintR, 'movmedian', smoothWindow, 'omitnan');
squintMean_s = smoothdata(squintMean, 'movmedian', smoothWindow, 'omitnan');

%% Squint asymmetry: left vs right
squintAsymmetry = squintL - squintR;
squintAsymmetry_s = smoothdata(squintAsymmetry, ...
    'movmedian', smoothWindow, 'omitnan');

asymThreshold = 0.15; % difference in squint index

leftMoreSquint = squintAsymmetry > asymThreshold;
rightMoreSquint = squintAsymmetry < -asymThreshold;

leftMoreSquint(isBlink) = false;
rightMoreSquint(isBlink) = false;

%% Threshold squinting
squintThreshold = 0.25; % aperture reduced by at least 25%
isSquintL = squintL > squintThreshold;
isSquintR = squintR > squintThreshold;
isSquintEither = isSquintL | isSquintR;
isSquintBoth = isSquintL & isSquintR;

%% Plot 1: raw/smoothed eyelid aperture
figure;
plot(t, apertureL_s, 'LineWidth', 1.2); hold on
plot(t, apertureR_s, 'LineWidth', 1.2);

if hasBlinkFile
    yl = ylim;
    for b = 1:height(blinks)
        patch([blinkStart(b) blinkEnd(b) blinkEnd(b) blinkStart(b)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.8 0.8 0.8], ...
              'FaceAlpha', 0.25, ...
              'EdgeColor', 'none');
    end
    plot(t, apertureL_s, 'LineWidth', 1.2);
    plot(t, apertureR_s, 'LineWidth', 1.2);
end

xlabel('Time from recording start (s)');
ylabel('Eyelid aperture (mm)');
title(['Left/Right Eyelid Aperture - ' recName]);
legend({'Left aperture','Right aperture','Blink period'}, 'Location', 'best');
grid on

%% Plot 2: squint index left vs right
figure;
plot(t, squintL_s, 'LineWidth', 1.2); hold on
plot(t, squintR_s, 'LineWidth', 1.2);
plot(t, squintMean_s, 'w', 'LineWidth', 1.4);
yline(squintThreshold, '--', 'Squint threshold');

xlabel('Time from recording start (s)');
ylabel('Squint index');
title(['Squint Index - ' recName]);
legend({'Left eye','Right eye','Mean','Threshold'}, 'Location', 'best');
grid on

%% Plot 3: pupil diameter with squint index
figure;
plot(t, pupilL, 'LineWidth', 1.1); hold on
plot(t, pupilR, 'LineWidth', 1.1);
xlabel('Time from recording start (s)');
ylabel('Pupil diameter (mm)');
title(['Pupil Diameter During Recording - ' recName]);
legend({'Left pupil','Right pupil'}, 'Location', 'best');
grid on

%% Plot 4: left-right squint asymmetry
figure;
plot(t, squintAsymmetry_s, 'LineWidth', 1.3); hold on
yline(0, '--');
yline(asymThreshold, '--', 'Left more squint');
yline(-asymThreshold, '--', 'Right more squint');

xlabel('Time from recording start (s)');
ylabel('Squint asymmetry index (Left - Right)');
title(['Left-Right Squint Asymmetry - ' recName]);
grid on

%% Summary
summaryStats = table( ...
    string(recName), ...
    baselineL, ...
    baselineR, ...
    median(apertureL_noBlink, 'omitnan'), ...
    median(apertureR_noBlink, 'omitnan'), ...
    mean(squintL, 'omitnan'), ...
    mean(squintR, 'omitnan'), ...
    mean(squintMean, 'omitnan'), ...
    100 * mean(isSquintL, 'omitnan'), ...
    100 * mean(isSquintR, 'omitnan'), ...
    100 * mean(isSquintEither, 'omitnan'), ...
    100 * mean(isSquintBoth, 'omitnan'), ...
    100 * mean(isBlink, 'omitnan'), ...
    mean(squintAsymmetry, 'omitnan'), ...
    std(squintAsymmetry, 'omitnan'), ...
    100 * mean(leftMoreSquint, 'omitnan'), ...
    100 * mean(rightMoreSquint, 'omitnan'), ...
    'VariableNames', {'Recording', ...
                      'BaselineApertureLeft_mm', ...
                      'BaselineApertureRight_mm', ...
                      'MedianApertureLeft_mm', ...
                      'MedianApertureRight_mm', ...
                      'MeanSquintIndexLeft', ...
                      'MeanSquintIndexRight', ...
                      'MeanSquintIndexBothEyes', ...
                      'PercentTimeSquintLeft', ...
                      'PercentTimeSquintRight', ...
                      'PercentTimeSquintEitherEye', ...
                      'PercentTimeSquintBothEyes', ...
                      'PercentTimeBlinking', ...
                      'MeanSquintAsymmetryLeftMinusRight', ...
                      'SDSquintAsymmetryLeftMinusRight', ...
                      'PercentTimeLeftMoreSquinted', ...
                      'PercentTimeRightMoreSquinted'} ...
);

disp(summaryStats)

outFile = ['squintingSummary_' recName '.csv'];
outFile = strrep(outFile, ' ', '_');
writetable(summaryStats, outFile);

fprintf('\nSaved summary to: %s\n', outFile);