%% plotSlidingWindowDepthVergence.m
% Run estimateDepthFromVergence for multiple activities and plot
% 5-second sliding-window depth and vergence summaries.

clear; close all; clc

%% Recording folders and 1 m calibration values
% calibVergenceDeg = median 3D optical-axis separation during known 1 m fixation

recordings = {
    'Dark',         '/Users/sophiamirabal/Downloads/timeseriesData/dark_data/flic_0021_dark_1-878d6eb2',                    XX;
    'Read',         '/Users/sophiamirabal/Downloads/timeseriesData/read_data/flic_0021_read_1-286bdbca',                    1.7038;
    'Chat',         '/Users/sophiamirabal/Downloads/timeseriesData/chat_data/flic_0021_chat_1-1e66beb5',                    7.1289;
    'Walk Indoor',  '/Users/sophiamirabal/Downloads/timeseriesData/walkIndoor_data/flic_0021_walkindoor_1-a061c060'         XX;
    'Walk Outdoor', '/Users/sophiamirabal/Downloads/timeseriesData/walkOutdoor_data/flic_0021_walkoutdoor_1-8cead05d',      XX;
    'Walk Biopond', '/Users/sophiamirabal/Downloads/timeseriesData/walkBiopond_data/flic_0021_walkbiopond_1-f90adb81',      XX; 
    'Sit Biopond',  '/Users/sophiamirabal/Downloads/timeseriesData/sitBiopond_data/flic_0021_walkbiopond_1-f90adb81',       XX;
};

%% Sliding window settings
windowSize_s = 5;      % 5-second basket/window
stepSize_s = 1;        % move window forward every 1 second
useMedian = true;      % median is more robust to spikes than mean

allWindowTables = cell(size(recordings,1),1);

%% Run depth estimation and sliding-window summaries
for r = 1:size(recordings,1)

    recName = recordings{r,1};
    recDir  = recordings{r,2};
    calibVergenceDeg = recordings{r,3};

    fprintf('\nProcessing %s...\n', recName);

    [depthTable, summaryStats] = estimateDepthFromVergence( ...
        recDir, ...
        calibVergenceDeg, ...
        recName);

    t = depthTable.Time_s;
    depth = depthTable.EstimatedDepth_m;
    vergence = depthTable.Vergence3D_deg;

    tStart = min(t);
    tEnd = max(t);

    windowStarts = tStart:stepSize_s:(tEnd - windowSize_s);
    nWindows = numel(windowStarts);

    windowCenter = NaN(nWindows,1);
    depthWindow = NaN(nWindows,1);
    vergenceWindow = NaN(nWindows,1);
    nSamplesWindow = NaN(nWindows,1);

    for w = 1:nWindows

        winStart = windowStarts(w);
        winEnd = winStart + windowSize_s;

        idx = t >= winStart & t < winEnd;

        windowCenter(w) = winStart + windowSize_s/2;
        nSamplesWindow(w) = sum(idx & ~isnan(depth) & ~isnan(vergence));

        if useMedian
            depthWindow(w) = median(depth(idx), 'omitnan');
            vergenceWindow(w) = median(vergence(idx), 'omitnan');
        else
            depthWindow(w) = mean(depth(idx), 'omitnan');
            vergenceWindow(w) = mean(vergence(idx), 'omitnan');
        end
    end

    windowTable = table( ...
        repmat(string(recName), nWindows, 1), ...
        windowStarts(:), ...
        windowStarts(:) + windowSize_s, ...
        windowCenter, ...
        depthWindow, ...
        vergenceWindow, ...
        nSamplesWindow, ...
        'VariableNames', {'Recording', ...
                          'WindowStart_s', ...
                          'WindowEnd_s', ...
                          'WindowCenter_s', ...
                          'DepthWindow_m', ...
                          'VergenceWindow_deg', ...
                          'NSamples'} ...
    );

    allWindowTables{r} = windowTable;

end

%% Combine all activities
combinedWindowTable = vertcat(allWindowTables{:});

%% Plot 1: sliding-window estimated depth
figure; hold on

for r = 1:size(recordings,1)

    recName = string(recordings{r,1});
    idx = combinedWindowTable.Recording == recName;

    plot(combinedWindowTable.WindowCenter_s(idx), ...
         combinedWindowTable.DepthWindow_m(idx), ...
         'LineWidth', 1.5);

end

yline(1, '--', '1 m calibration depth');

xlabel('Time from recording start (s)');
ylabel('Estimated depth (m)');
title('5-second Sliding-Window Estimated Depth');
legend(recordings(:,1), 'Location', 'best');
grid on

%% Plot 2: sliding-window vergence
figure; hold on

for r = 1:size(recordings,1)

    recName = string(recordings{r,1});
    idx = combinedWindowTable.Recording == recName;

    plot(combinedWindowTable.WindowCenter_s(idx), ...
         combinedWindowTable.VergenceWindow_deg(idx), ...
         'LineWidth', 1.5);

end

xlabel('Time from recording start (s)');
ylabel('3D optical-axis separation (deg)');
title('5-second Sliding-Window Vergence');
legend(recordings(:,1), 'Location', 'best');
grid on

%% Save combined sliding-window table
outFile = 'slidingWindow_depthVergence_5s.csv';
writetable(combinedWindowTable, outFile);

fprintf('\nSaved sliding-window table to: %s\n', outFile);