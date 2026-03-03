function plotGazeAccuracyNeon(startIds, endIds, filePaths)
%{
    startIds = [3, 45, 5, 36, 43];
    endIds = [82, 163, 52, 108, 120];
    filePaths = { ...
        '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-13/2026-03-03_15-39-46-286215af/gaze.csv', ...
        '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-14/2026-03-03_16-14-03-169f6719/gaze.csv', ...
        '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-11/2026-03-03_15-22-52-7675bb73/gaze.csv', ...
        '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-12/2026-03-03_15-31-40-6c501c11/gaze.csv', ...
        '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-10/2026-03-02_15-20-42-4a242290/gaze.csv' ...
    };
    plotGazeAccuracyNeon(startIds, endIds, filePaths)
%}

numParticipants = length(startIds);
maxTargets = 34;
allActualGaze = NaN(maxTargets, 2, numParticipants);
baseTargets = [0,0; -15,15; -15,-15; 15,15; 15,-15; 0,15; 0,-15; -15,0; 15,0; ...
    -7.5,7.5; -7.5,-7.5; 7.5,7.5; 7.5,-7.5; 0,7.5; 0,-7.5; -7.5,0; 7.5,0];
masterIntended = [baseTargets; baseTargets];

% 1. Setup the main UI figure and Tab Group
hFig = figure('Position', [100 100 900 800], 'Name', 'Gaze Accuracy Analysis');
tGroup = uitabgroup('Parent', hFig);

fprintf('Starting analysis for %d participants...\n', numParticipants);
for i = 1:numParticipants
    fprintf('Processing participant %d...\n', i);
    [actGaze, intGaze] = processParticipant(startIds(i), endIds(i), filePaths{i});

    biasX = mean(actGaze(:,1) - intGaze(:,1), 'omitnan');
    biasY = mean(actGaze(:,2) - intGaze(:,2), 'omitnan');
    actGaze_corrected = actGaze - [biasX, biasY];

    numFound = size(actGaze_corrected, 1);
    allActualGaze(1:numFound, :, i) = actGaze_corrected;

    % 2. Create a Tab for this participant
    hTab = uitab('Parent', tGroup, 'Title', sprintf('P%d', i));
    drawGazePlot(hTab, intGaze, actGaze_corrected, sprintf('Participant %d (Avg. Offset: X=%.2f, Y=%.2f)', i, biasX, biasY));
end

% 3. Create a final Tab for the Grand Average
grandAvgActual = mean(allActualGaze, 3, 'omitnan');
avgTab = uitab('Parent', tGroup, 'Title', 'Grand Avg');
drawGazePlot(avgTab, masterIntended, grandAvgActual, 'Grand Average Across All Participants');
end

function drawGazePlot(parent, intendedGaze, actualGaze, figTitle)
% Create axes specifically inside the passed tab
ax = axes('Parent', parent);
hold(ax, 'on');

numFound = size(actualGaze, 1);
xOffsets = actualGaze(:,1) - intendedGaze(:,1);
yOffsets = actualGaze(:,2) - intendedGaze(:,2);
totalErr = mean(sqrt(xOffsets.^2 + yOffsets.^2), 'omitnan');

% Plotting
hInt = scatter(ax, intendedGaze(:,1), intendedGaze(:,2), 150, 'r', 'x', 'LineWidth', 2);
hAct = scatter(ax, actualGaze(:,1), actualGaze(:,2), 100, 'cyan', 'filled', 'MarkerEdgeColor', 'k');

for i = 1:numFound
    if ~isnan(actualGaze(i,1))
        line(ax, [intendedGaze(i,1), actualGaze(i,1)], [intendedGaze(i,2), actualGaze(i,2)], ...
            'Color', [0.5 0.5 0.5], 'LineStyle', '-');
        text(ax, actualGaze(i,1) + 0.6, actualGaze(i,2) + 0.6, num2str(i), ...
            'Color', 'k', 'FontSize', 8);
    end
end

title(ax, {figTitle, sprintf('Mean Error: %.2f°', totalErr)});
xlabel(ax, 'Azimuth (degrees)'); ylabel(ax, 'Elevation (degrees)');
legend(ax, [hInt, hAct], {'Intended Target', 'Actual Gaze'}, 'Location', 'northeastoutside');
grid(ax, 'on'); axis(ax, 'equal');
xlim(ax, [-22 22]); ylim(ax, [-22 22]);
end

% (Keep your existing processParticipant function exactly as it was)

function [actualGaze, intendedGaze] = processParticipant(startId, endId, filePath)
% 1. FORCE IMPORT: Read as text to bypass "NaN" heuristic trap
opts = detectImportOptions(filePath, 'FileType', 'text');
opts = setvartype(opts, 'string');
data = readtable(filePath, opts);

% 2. DYNAMIC COLUMN FINDING
headers = lower(data.Properties.VariableNames);
idx_fId = find(contains(headers, 'fixation'), 1);
idx_ts = find(contains(headers, 'timestamp'), 1);
idx_az = find(contains(headers, 'azimuth'), 1);
idx_el = find(contains(headers, 'elevation'), 1);

if isempty(idx_fId) || isempty(idx_ts) || isempty(idx_az) || isempty(idx_el)
    error('Could not map columns. Found these headers: %s', strjoin(headers, ', '));
end

% 3. EXTRACT AND CONVERT
fId = str2double(data{:, idx_fId});
ts_ns = str2double(data{:, idx_ts});
az = str2double(data{:, idx_az});
el = str2double(data{:, idx_el});

% 4. EXTRACT RANGE
startIndex = find(fId == startId, 1, 'first');
endIndex = find(fId == endId, 1, 'last');

if isempty(startIndex) || isempty(endIndex)
    error('Fixation IDs %d-%d not found in this file.', startId, endId);
end

ts_ns = ts_ns(startIndex:endIndex);
fId = fId(startIndex:endIndex);
az = az(startIndex:endIndex);
el = el(startIndex:endIndex);

% 5. CORE PROCESSING LOGIC
baseTargets = [0,0; -15,15; -15,-15; 15,15; 15,-15; 0,15; 0,-15; -15,0; 15,0; ...
    -7.5,7.5; -7.5,-7.5; 7.5,7.5; 7.5,-7.5; 0,7.5; 0,-7.5; -7.5,0; 7.5,0];
intendedGaze = [baseTargets; baseTargets];

isValid = ~isnan(fId);
diffValid = diff([0; isValid; 0]);
blockStarts = find(diffValid == 1);
blockEnds = find(diffValid == -1) - 1;

actualGaze = [];

% Initialize with first block
currTargetAz = az(blockStarts(1):blockEnds(1));
currTargetEl = el(blockStarts(1):blockEnds(1));
targetStartTime = ts_ns(blockStarts(1));

for i = 2:length(blockStarts)
    thisBlockAz = median(az(blockStarts(i):blockEnds(i)), 'omitnan');
    thisBlockEl = median(el(blockStarts(i):blockEnds(i)), 'omitnan');

    prevMedAz = median(currTargetAz, 'omitnan');
    prevMedEl = median(currTargetEl, 'omitnan');

    dist = sqrt((thisBlockAz - prevMedAz)^2 + (thisBlockEl - prevMedEl)^2);
    timeElapsed = (ts_ns(blockStarts(i)) - targetStartTime) / 1e9;

    if dist > 3.5 && timeElapsed > 2.5
        actualGaze = [actualGaze; median(currTargetAz, 'omitnan'), median(currTargetEl, 'omitnan')];
        currTargetAz = az(blockStarts(i):blockEnds(i));
        currTargetEl = el(blockStarts(i):blockEnds(i));
        targetStartTime = ts_ns(blockStarts(i));
    else
        currTargetAz = [currTargetAz; az(blockStarts(i):blockEnds(i))];
        currTargetEl = [currTargetEl; el(blockStarts(i):blockEnds(i))];
    end
end

actualGaze = [actualGaze; median(currTargetAz, 'omitnan'), median(currTargetEl, 'omitnan')];

% Truncate to match intended targets
if size(actualGaze, 1) > 34
    actualGaze = actualGaze(1:34, :);
end

intendedGaze = intendedGaze(1:size(actualGaze, 1), :);
end