%% analyzeSquinting_darkBaseline_allTasks.m
% Use dark recording as open-eye baseline, then compare aperture across tasks.

clear; close all; clc

%% Recordings
recordings = {
    'Dark',         '/Users/sophiamirabal/Downloads/timeseriesData/dark_data/flic_0021_dark_1-878d6eb2';
    'Read',         '/Users/sophiamirabal/Downloads/timeseriesData/read_data/flic_0021_read_1-286bdbca';
    'Chat',         '/Users/sophiamirabal/Downloads/timeseriesData/chat_data/flic_0021_chat_1-1e66beb5';
    'Walk Indoor',  '/Users/sophiamirabal/Downloads/timeseriesData/walkIndoor_data/flic_0021_walkindoor_1-a061c060';
    'Walk Outdoor', '/Users/sophiamirabal/Downloads/timeseriesData/walkOutdoor_data/flic_0021_walkoutdoor_1-8cead05d';
    'Walk Biopond', '/Users/sophiamirabal/Downloads/timeseriesData/walkBiopond_data/flic_0021_walkbiopond_1-f90adb81';
    'Sit Biopond',  '/Users/sophiamirabal/Downloads/timeseriesData/sitBiopond_data/flic_0021_sitbiopond_1-dd70c117';
};

%% Trim settings
% Ignore starts/ends where camera is being adjusted, task not underway, etc.
trimStart_s = 90;
trimEnd_s   = 90;

smoothWindow = 15;
baselinePercentile = 95; % use high percentile from dark as "open eye"

%% First: get dark baseline
darkDir = recordings{1,2};
darkEye = readtable(fullfile(darkDir, '3d_eye_states.csv'), ...
    'VariableNamingRule', 'preserve');

t0_dark = darkEye.("timestamp [ns]")(1);
tDark = (darkEye.("timestamp [ns]") - t0_dark) * 1e-9;

apDarkL = darkEye.("eyelid aperture left [mm]");
apDarkR = darkEye.("eyelid aperture right [mm]");

blinkFileDark = fullfile(darkDir, 'blinks.csv');
isBlinkDark = false(size(tDark));

if isfile(blinkFileDark)
    blinksDark = readtable(blinkFileDark, 'VariableNamingRule', 'preserve');
    blinkStart = (blinksDark.("start timestamp [ns]") - t0_dark) * 1e-9;
    blinkEnd   = (blinksDark.("end timestamp [ns]") - t0_dark) * 1e-9;

    for b = 1:height(blinksDark)
        isBlinkDark = isBlinkDark | ...
            (tDark >= blinkStart(b) & tDark <= blinkEnd(b));
    end
end

idxDarkTrim = tDark >= trimStart_s & ...
              tDark <= max(tDark) - trimEnd_s & ...
              ~isBlinkDark;

darkBaselineL = prctile(apDarkL(idxDarkTrim), baselinePercentile);
darkBaselineR = prctile(apDarkR(idxDarkTrim), baselinePercentile);

fprintf('\nDark baseline aperture:\n');
fprintf('Left eye %.3f mm\n', darkBaselineL);
fprintf('Right eye %.3f mm\n', darkBaselineR);

summary = table();
allLeftRightTables = cell(size(recordings,1),1);

%% Analyze each recording relative to dark baseline
figure; hold on
taskColors = lines(size(recordings,1));

for r = 1:size(recordings,1)

    recName = recordings{r,1};
    recDir  = recordings{r,2};

    eye = readtable(fullfile(recDir, '3d_eye_states.csv'), ...
        'VariableNamingRule', 'preserve');

    t0 = eye.("timestamp [ns]")(1);
    t = (eye.("timestamp [ns]") - t0) * 1e-9;

    apertureL = eye.("eyelid aperture left [mm]");
    apertureR = eye.("eyelid aperture right [mm]");

    %% Blink mask
    isBlink = false(size(t));
    blinkFile = fullfile(recDir, 'blinks.csv');

    if isfile(blinkFile)
        blinks = readtable(blinkFile, 'VariableNamingRule', 'preserve');
        blinkStart = (blinks.("start timestamp [ns]") - t0) * 1e-9;
        blinkEnd   = (blinks.("end timestamp [ns]") - t0) * 1e-9;

        for b = 1:height(blinks)
            isBlink = isBlink | ...
                (t >= blinkStart(b) & t <= blinkEnd(b));
        end
    end

    %% Trim usable middle of task
    idxUse = t >= trimStart_s & ...
             t <= max(t) - trimEnd_s & ...
             ~isBlink;

    %% Relative aperture to dark baseline
    relativeApertureL = apertureL ./ darkBaselineL;
    relativeApertureR = apertureR ./ darkBaselineR;

    squintL = 1 - relativeApertureL;
    squintR = 1 - relativeApertureR;
    squintMean = mean([squintL squintR], 2, 'omitnan');

    relativeApertureMean = mean([relativeApertureL relativeApertureR], 2, 'omitnan');

    %% Sliding-window relative aperture

    windowSize_s = 5;      % seconds
    stepSize_s = 1;        % seconds
    
    windowStarts = trimStart_s:stepSize_s:(max(t)-trimEnd_s-windowSize_s);
    
    windowCenter = windowStarts + windowSize_s/2;
    
    relativeApertureWindow = NaN(size(windowCenter));
    relativeApertureLeftWindow = NaN(size(windowCenter));
    relativeApertureRightWindow = NaN(size(windowCenter));
    
    for w = 1:numel(windowStarts)
    
        idxWindow = ...
            t >= windowStarts(w) & ...
            t < windowStarts(w)+windowSize_s & ...
            idxUse;
    
        relativeApertureWindow(w) = ...
            median(relativeApertureMean(idxWindow),'omitnan');
    
        relativeApertureLeftWindow(w) = ...
            median(relativeApertureL(idxWindow),'omitnan');
    
        relativeApertureRightWindow(w) = ...
            median(relativeApertureR(idxWindow),'omitnan');
    
    end
    
    plot(windowCenter, ...
         relativeApertureWindow, ...
         'Color',taskColors(r,:), ...
         'LineWidth',2);
    
    leftRightTable = table( ...
        repmat(string(recName), numel(windowCenter), 1), ...
        windowCenter(:), ...
        relativeApertureLeftWindow(:), ...
        relativeApertureRightWindow(:), ...
        relativeApertureWindow(:), ...
        'VariableNames', {'Recording', ...
                          'WindowCenter_s', ...
                          'RelativeApertureLeft', ...
                          'RelativeApertureRight', ...
                          'RelativeApertureMean'} ...
    );
    
    allLeftRightTables{r} = leftRightTable;

    %% Summary row
    newRow = table( ...
        string(recName), ...
        median(relativeApertureL(idxUse), 'omitnan'), ...
        median(relativeApertureR(idxUse), 'omitnan'), ...
        median(relativeApertureMean(idxUse), 'omitnan'), ...
        median(squintL(idxUse), 'omitnan'), ...
        median(squintR(idxUse), 'omitnan'), ...
        median(squintMean(idxUse), 'omitnan'), ...
        100 * mean(squintL(idxUse) > 0.25, 'omitnan'), ...
        100 * mean(squintR(idxUse) > 0.25, 'omitnan'), ...
        100 * mean(isBlink, 'omitnan'), ...
        sum(idxUse), ...
        'VariableNames', {'Recording', ...
                          'MedianRelativeApertureLeft', ...
                          'MedianRelativeApertureRight', ...
                          'MedianRelativeApertureMean', ...
                          'MedianSquintLeft', ...
                          'MedianSquintRight', ...
                          'MedianSquintMean', ...
                          'PercentTimeSquintLeft', ...
                          'PercentTimeSquintRight', ...
                          'PercentTimeBlinking', ...
                          'NSamplesUsed'} ...
    );

    summary = [summary; newRow];

end

yline(1, '--', 'Dark baseline');
yline(0.75, '--', '25% reduction');

xlabel('Time from recording start (s)');
ylabel('Aperture relative to dark baseline');
title('Eyelid Aperture Relative to Dark Open-Eye Baseline');
legend(recordings(:,1), 'Location', 'best');
grid on

%% Plot separate left and right eye openness

combinedLeftRightTable = vertcat(allLeftRightTables{:});

figure; hold on

for r = 1:size(recordings,1)

    recName = string(recordings{r,1});
    idx = combinedLeftRightTable.Recording == recName;

    thisColor = taskColors(r,:);

    plot(combinedLeftRightTable.WindowCenter_s(idx), ...
        combinedLeftRightTable.RelativeApertureLeft(idx), ...
        '--', ...
        'Color', thisColor, ...
        'LineWidth', 1.5);

    plot(combinedLeftRightTable.WindowCenter_s(idx), ...
        combinedLeftRightTable.RelativeApertureRight(idx), ...
        '-', ...
        'Color', thisColor, ...
        'LineWidth', 1.5);

end

yline(1, '--', 'Dark baseline');
yline(0.75, '--', '25% reduction');

xlabel('Time from recording start (s)');
ylabel('Aperture relative to dark baseline');
title('Separate Left and Right Eye Openness Relative to Dark Baseline');

legend({'Dark left','Dark right', ...
    'Read left','Read right', ...
    'Chat left','Chat right', ...
    'Walk Indoor left','Walk Indoor right', ...
    'Walk Outdoor left','Walk Outdoor right', ...
    'Walk Biopond left','Walk Biopond right', ...
    'Sit Biopond left','Sit Biopond right'}, ...
    'Location', 'best');

grid on

%% Plot summary bar chart: left, right, and mean aperture

figure;

barData = [ ...
    summary.MedianRelativeApertureLeft, ...
    summary.MedianRelativeApertureRight, ...
    summary.MedianRelativeApertureMean];

bar(barData, 'grouped');

xticks(1:height(summary));
xticklabels(summary.Recording);
xtickangle(30);

ylabel('Median relative aperture');
title('Median Eyelid Aperture Relative to Dark Baseline');

yline(1, '--', 'Dark baseline');
yline(0.75, '--', '25% reduction');

legend({'Left eye','Right eye','Mean of both'}, ...
       'Location','best');

grid on

%% Save
disp(summary)

outFile = 'squinting_darkBaseline_allTasks_summary.csv';
writetable(summary, outFile);

fprintf('\nSaved summary to: %s\n', outFile);