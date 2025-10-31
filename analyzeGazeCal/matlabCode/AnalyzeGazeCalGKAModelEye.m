%% Initial Data Loading and Cleaning
load('/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/sam_gazecal_106_pupilData.mat')

% --- CONTROL PARAMETER ---
x_axis_mode = 'seconds'; % Set to 'seconds' or 'frames'

% --- Fixed Parameters ---
fps = 120;
dotTime_s = 3.488;         % Duration of each target
startFrame = 1.0313e+04;       % Frame where the targets start (ORIGINAL REFERENCE)
time_offset_s = 0.341667; % TIME OFFSET to set T=0 at the start of the first dot
preTaskCutoff_s = 1;      % Time to display before the start frame (1 second)
RMSECutoff = 2;           
azimuthColor = [0 0.447 0.741]; 
elevationColor = [0.85 0.325 0.098]; 

% --- Gaze Analysis Parameters ---
window_s = 1.0;                  % Analysis window duration
windowFrames = round(window_s * fps);
halfWindow = floor(windowFrames / 2);

% --- Plotting Colors for Averaged Gaze ---
repeat1Color = [0.8 0.2 0.2];     % Light red for Repeat 1 (Dots 1-5)
repeat2Color = [0.5 0.0 0.0];     % Dark red for Repeat 2 (Dots 6-10)

% --- Data Filtering and Time Alignment ---

% 1. Create indices for bad data (logical OR)
badRMSEIdx = pupilData.sceneConstrained.ellipses.RMSE > RMSECutoff;
badFitAtBoundIdx = pupilData.sceneConstrained.eyePoses.fitAtBound > 0;
badIdx = badRMSEIdx | badFitAtBoundIdx; 
filteredGazeX = pupilData.sceneConstrained.eyePoses.values(:, 1);
filteredGazeY = pupilData.sceneConstrained.eyePoses.values(:, 2);

% 2. Define the plotting range based on the startFrame and cutoff
framesBeforeStart = round(preTaskCutoff_s * fps);
plotStartFrame = max(1, startFrame - framesBeforeStart);
nTotalFrames = size(pupilData.sceneConstrained.eyePoses.values, 1);

% 3. Create time/frame vectors for plotting
timeVector_frames_raw = (plotStartFrame : nTotalFrames)';
timeVector_s_raw = (timeVector_frames_raw - startFrame) / fps;
timeVector_s = timeVector_s_raw - time_offset_s; 

% 4. Select the correct X-vector for plotting and define the label
if strcmp(x_axis_mode, 'frames')
    X_PLOT_VECTOR = timeVector_frames_raw;
    X_LABEL = 'Frame Number (Absolute)';
else % Default to 'seconds'
    X_PLOT_VECTOR = timeVector_s;
    X_LABEL = 'Time (seconds, T=0 at First Dot Start)';
end

% 5. Slice the actual gaze data and the bad data mask to the new time range
gazeX_full_sliced = filteredGazeX(plotStartFrame:nTotalFrames);
gazeY_full_sliced = filteredGazeY(plotStartFrame:nTotalFrames);
badIdx_sliced = badIdx(plotStartFrame:nTotalFrames);

% --- Prepare Gaze Target Data ---
% 5-dot sequence repeated twice = 10 total targets
deg_positions_raw = [-1,1].*[0, 0; -20, 20; -20, -20; 20, 20; 20, -20; ...
            0, 20; 0, -20; -20, 0; 20, 0;...
            -10, 10; -10, -10; 10, 10; 10, -10; ...
            0, 10; 0, -10; -10, 0; 10, 0];
deg_positions = repmat(deg_positions_raw, 2, 1); 
nTotalDots = size(deg_positions, 1);

% --- Create Known Target Time Series (Step Function for Figure 1) ---
nPlotFrames = length(X_PLOT_VECTOR);
knownTargets = zeros(nPlotFrames, 2);

framesPerDot = round(dotTime_s * fps);
fullDataEndFrame = startFrame + (nTotalDots * framesPerDot) - 1;

currentFrame_abs = startFrame; 

for i = 1:nTotalDots
    % ... (Target series creation loop - unchanged)
    idxStart_abs = currentFrame_abs;
    idxEnd_abs = min(currentFrame_abs + framesPerDot - 1, fullDataEndFrame);
    
    if idxEnd_abs >= plotStartFrame
        idxStart_sliced = max(1, idxStart_abs - plotStartFrame + 1);
        idxEnd_sliced = idxEnd_abs - plotStartFrame + 1;
        idxEnd_sliced = min(idxEnd_sliced, nPlotFrames);
        
        if idxStart_sliced <= idxEnd_sliced
            numFramesInSegment = idxEnd_sliced - idxStart_sliced + 1;
            targetValue = deg_positions(i, :);
            knownTargets(idxStart_sliced:idxEnd_sliced, :) = repmat(targetValue, numFramesInSegment, 1); 
        end
    end
    currentFrame_abs = idxEnd_abs + 1;
    if currentFrame_abs > nTotalFrames; break; end
end


%% -------------------------------------------------------------------
%  SECTION 2: GAZE AVERAGING AND COMPARISON PLOT
% -------------------------------------------------------------------

% Preallocate arrays to store results
avgMeasuredGaze = zeros(nTotalDots, 2); % [Avg Gaze X, Avg Gaze Y]

% Recalculate the absolute center frame for each dot presentation
currentFrame_abs = startFrame; 

for i = 1:nTotalDots
    
    % 1. Determine Target Window
    targetStart_abs = currentFrame_abs;
    targetEnd_abs = targetStart_abs + framesPerDot - 1;
    
    % 2. Find Center of Target Presentation
    centerFrame = round((targetStart_abs + targetEnd_abs) / 2);
    
    % 3. Define 1-Second Analysis Window centered on the presentation
    windowStart_abs = max(1, centerFrame - halfWindow);
    windowEnd_abs = min(nTotalFrames, centerFrame + halfWindow - 1);
    
    % 4. Slice Gaze Data for the Window
    windowIndices = windowStart_abs : windowEnd_abs;
    
    % Get raw data for the window
    gazeX_window = filteredGazeX(windowIndices);
    gazeY_window = filteredGazeY(windowIndices);
    
    % Get bad data flags for the window
    badIdx_window = badIdx(windowIndices);
    
    % 5. Apply Filter (Set Bad Data to NaN)
    gazeX_window(badIdx_window) = NaN;
    gazeY_window(badIdx_window) = NaN;
    
    % 6. Calculate Average (excluding NaNs)
    avgMeasuredGaze(i, 1) = nanmean(gazeX_window);
    avgMeasuredGaze(i, 2) = nanmean(gazeY_window);
    
    % Update for next target
    currentFrame_abs = targetEnd_abs + 1;
end

% --- Create Figure 2: Known vs. Measured ---
figure('Units','normalized','Position',[0.55 0.1 0.4 0.8]);
ax3 = axes;
hold on;
title(ax3, 'Gaze Calibration Accuracy: Known vs. Measured');
xlabel(ax3, 'Horizontal Gaze (Azimuth, deg)');
ylabel(ax3, 'Vertical Gaze (Elevation, deg)');
grid on;
axis equal; % Ensure X and Y scales are visually matched

% Plot Known Targets (All 10 dots, using 'o' markers)
plot(deg_positions(:, 1), deg_positions(:, 2), 'o', ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ... % Grey for known targets
    'MarkerFaceColor', [0.8 0.8 0.8], ... % Light fill
    'MarkerSize', 8, ...
    'DisplayName', 'Known Target Position');

% Plot Averaged Measured Gaze - Repeat 1 (Dots 1-5)
plot(avgMeasuredGaze(1:5, 1), avgMeasuredGaze(1:5, 2), 'x', ...
    'MarkerEdgeColor', repeat1Color, ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'DisplayName', 'Averaged Gaze (Repeat 1)');

% Plot Averaged Measured Gaze - Repeat 2 (Dots 6-10)
plot(avgMeasuredGaze(6:10, 1), avgMeasuredGaze(6:10, 2), 'x', ...
    'MarkerEdgeColor', repeat2Color, ...
    'LineWidth', 2, ...
    'MarkerSize', 10, ...
    'DisplayName', 'Averaged Gaze (Repeat 2)');

% Add legend and adjust limits slightly
legend(ax3, 'Location', 'best');
padding = 5; % Add a 5-degree buffer to the plot limits
xlim([min(deg_positions(:, 1)) - padding, max(deg_positions(:, 1)) + padding]);
ylim([min(deg_positions(:, 2)) - padding, max(deg_positions(:, 2)) + padding]);

% --- Figure 1 Plotting (Original time series, unchanged) ---

figure('Units','normalized','Position',[0.1 0.1 0.4 0.8]);
t = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- Tile 1 (Top): Actual Gaze Angles ---
ax1 = nexttile;
hold on;
gazeX_plot = gazeX_full_sliced;
gazeX_plot(badIdx_sliced) = NaN;
plot(X_PLOT_VECTOR, gazeX_plot, 'Color', azimuthColor, 'LineWidth', 1, 'DisplayName', 'Gaze X (Filtered Azimuth)');

gazeY_plot = gazeY_full_sliced;
gazeY_plot(badIdx_sliced) = NaN;
plot(X_PLOT_VECTOR, gazeY_plot, 'Color', elevationColor, 'LineWidth', 1, 'DisplayName', 'Gaze Y (Filtered Elevation)');

hold off;
title(ax1, 'Measured Gaze Angles');
ylabel(ax1, 'Gaze Angle (deg)');
legend(ax1, 'Location', 'best');
grid on;

% --- Tile 2 (Bottom): Known Gaze Targets ---
ax2 = nexttile;
hold on;
plot(X_PLOT_VECTOR, knownTargets(:, 1), 'Color', azimuthColor, 'LineWidth', 2, 'DisplayName', 'Target X (Azimuth)');
plot(X_PLOT_VECTOR, knownTargets(:, 2), 'Color', elevationColor, 'LineWidth', 2, 'DisplayName', 'Target Y (Elevation)');
hold off;
title(ax2, 'Known Gaze Target Sequence');
xlabel(ax2, X_LABEL);
ylabel(ax2, 'Target Angle (deg)');
legend(ax2, 'Location', 'best');
grid on;

% --- Match Axes Properties ---
linkaxes([ax1, ax2], 'x');
linkaxes([ax1, ax2], 'y');