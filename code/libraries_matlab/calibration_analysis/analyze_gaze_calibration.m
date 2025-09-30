function analyze_gaze_calibration(confidence_cutoff, dot1_start_time)
    arguments
        % Set a less conservative default cutoff (e.g., 0.80)
        confidence_cutoff (1,1) double = 0.80
        % UPDATED: Setting the default Dot 1 start time to 0 seconds.
        dot1_start_time (1,1) double = 0.0
    end 
    % Load in the extracted pupil features
    dropBoxBaseDir=getpref('combiExperiments','dropboxBaseDir');
    % CORRECTED: Using the correct file path
    pupil_features = load(fullfile(dropBoxBaseDir, '/FLIC_data/lightLogger/SM_gaze_cal_pupil_features0929.mat')).pupil_features;
    
    % UPDATED: Restored the center position [0, 0] at the start of the list. 
    deg_positions = [ ...
        0, 0;   -20, 20;   -20, -20;   20, 20;   20, -20; ...
        0, 20;   0, -20;   -20,   0;   20, 0; ...
        -15, 15;  15, 15;   -15, -15;   15, -15; ...
        -10,  10;   -10, -10;    10,  10;    10, -10; ...
        0,  10;    0, -10;   -10,   0;     10,   0; ...
        -5,   5;   5,   5;   -5,  -5;     5,  -5]; % Now 24 unique positions
        
    % Retrieve the features we want to analyze 
    [gaze_angles, confidence_measures, pupil_t] = flatten_features(pupil_features);
    
    % Ensure we generated the correct number of points 
    assert(size(gaze_angles, 1) == size(pupil_t, 1), "Incorrect number of t values for datapoints"); 
    
    % Display the data after making desired transformations 
    [avg_gaze_angles, positionsRepeated, experimentEndTime, dotDuration, window_starts, window_ends] = display_adjusted_results(gaze_angles, confidence_measures, pupil_t, deg_positions, confidence_cutoff, dot1_start_time); 
    
    % Create the new Gaze Calibration Figure
    create_calibration_figure(avg_gaze_angles, positionsRepeated, dotDuration);
    
    return ; 
end 

% -------------------------------------------------------------------------
% LOCAL FUNCTIONS START HERE
% -------------------------------------------------------------------------

% Local function to extract the desired features from the list 
function [gaze_angles, confidence_measures, pupil_t] = flatten_features(pupil_features)
    % Initialize return arrays for the features we will 
    % flatten
    gaze_angles = nan(size(pupil_features, 1), 2); 
    confidence_measures = nan(size(pupil_features, 1), 2); 
    pupil_t = nan(size(pupil_features, 1), 1); 
    % Iterate over the pupil features structs 
    % and extract the ones we care about 
    for ii = 1:numel(pupil_features)
        % Retrieve a given frame's features 
        frame_features = pupil_features{ii};
        % Extract timestamp 
        pupil_t(ii, :) = frame_features.timestamp; 
        
        % Save phi and theta
        gaze_angles(ii, :) = rad2deg([frame_features.phi, frame_features.theta]); 
        % Save the confidence metrics
        confidence_measures(ii, :) = [frame_features.confidence, frame_features.model_confidence];
    end 
    return ; 
end 

% -------------------------------------------------------------------------

% Local function to display results after clipping to 
% only gaze calibration period
% mean centering
% removing low confidence points
function [avg_gaze_angles, positionsRepeated, experimentEndTime, dotDuration, window_starts, window_ends] = display_adjusted_results(gaze_angles_raw, confidence_measures_raw, pupil_t_raw, deg_positions, confidence_cutoff, dot1_start_time)
    
    % --- EXPERIMENT PARAMETERS INITIALIZATION & TIMING RECALCULATION ---
    repeats = 3;
    % Dot duration is 2.45 seconds
    dotDuration = 2.45; 
    interRepDelay = 0; 
    nDotsPerRep = size(deg_positions, 1); 
    
    % Calculate target appearance times (relative to the first dot's start at t=0)
    appearanceTimes = getDotAppearanceTimes_updated(dotDuration, nDotsPerRep, repeats, interRepDelay);
    
    % Initialize positionsRepeated (72 total targets)
    positionsRepeated = repmat(deg_positions, repeats, 1);
    
    % --- TIMING ADJUSTMENT (To align Dot 1 Start to dot1_start_time) ---
    
    % RAW START TIME: Based on user input (67.833s)
    raw_time_of_dot1_start = 67.833; 
    
    % Calculate the overall offset to apply to all raw timestamps:
    % New Offset = Raw Time of Dot 1 - Desired Plot Time of Dot 1 (now 0.0s by default)
    time_offset_to_apply = raw_time_of_dot1_start - dot1_start_time;
    
    % We need to cut the data starting from some point before Dot 1.
    % We'll use the raw time corresponding to the plot time of -8.0s for the cut-off.
    cut_off_plot_time = -8.0;
    cut_off_raw_time = time_offset_to_apply + cut_off_plot_time;
    
    % 1. Cut out the padding seconds of recording
    [~, start_idx] = min(abs(pupil_t_raw - cut_off_raw_time));
    pupil_t_all = pupil_t_raw(start_idx:end, :); 
    gaze_angles_all = gaze_angles_raw(start_idx:end, :); 
    confidence_measures_all = confidence_measures_raw(start_idx:end, :); 
    
    % Correct time for plotting so that the *entire plot* aligns with the desired dot1_start_time
    pupil_t_corrected_all = pupil_t_all - time_offset_to_apply; 
    
    % Adjust appearance times (relative to the new corrected timebase)
    appearanceTimes_adjusted = appearanceTimes + dot1_start_time; % appearanceTimes starts at 0, adding dot1_start_time (now 0) shifts all dot times.
    
    % Determine the expected end time of the experiment (for plot limits)
    expectedTaskEnd = appearanceTimes_adjusted(end) + dotDuration;
    experimentEndTime = expectedTaskEnd + 2; % Add a small padding (2s) for visualization
    
    % --- DATA CLEANING AND CENTERING ---
    
    % Store the count BEFORE filtering for confidence
    initial_data_count = size(gaze_angles_raw, 1); % Use raw data for this count
    
    % Apply confidence filter for *averaging and for gaze angle plots*
    good_idx_for_filtering = (confidence_measures_all(:, 1) >= confidence_cutoff) & (confidence_measures_all(:, 2) >= confidence_cutoff);
    
    % Use only high-confidence points for mean-centering and averaging
    gaze_angles_high_confidence = gaze_angles_all(good_idx_for_filtering, :);
    pupil_t_high_confidence = pupil_t_corrected_all(good_idx_for_filtering, :);
    
    % Store the count AFTER filtering
    filtered_data_count = size(gaze_angles_high_confidence, 1);
    
    % Calculate the percentage of data points tossed from the *averaging/transformation*
    tossed_count = initial_data_count - filtered_data_count;
    tossed_percentage = (tossed_count / initial_data_count) * 100;
    
    % 3. Subtract the mean theta from the theta vector, and the mean phi from the phi vector
    % Calculate mean from *all* available data (or filtered, based on preference - here using initial segment)
    mean_phi_theta = mean(gaze_angles_raw(start_idx:end, :), 1); 

    % Apply mean-centering to the high-confidence gaze angles
    gaze_angles_high_confidence_centered = gaze_angles_high_confidence;
    gaze_angles_high_confidence_centered(:, 1) = gaze_angles_high_confidence_centered(:, 1) - mean_phi_theta(1);
    gaze_angles_high_confidence_centered(:, 2) = gaze_angles_high_confidence_centered(:, 2) - mean_phi_theta(2); 
    
    % Calculate the average gaze angles for each dot
    [midpoints_shifted, avg_gaze_angles_full] = average_gaze_middle_segment_updated(gaze_angles_high_confidence_centered, pupil_t_high_confidence, appearanceTimes_adjusted, dotDuration);
    
    % --- DATA RETENTION: KEEP ALL 72 DOTS ---
    avg_gaze_angles = avg_gaze_angles_full; 
    midpoints = midpoints_shifted; 
    
    % --- PLOT MODIFICATION FOR STEP FUNCTION (Tile 2) ---
    
    full_positions_repeated = repmat(deg_positions, repeats, 1);
    
    % Step times use the adjusted appearance times (which start at dot1_start_time)
    stepTimes = reshape([appearanceTimes_adjusted, appearanceTimes_adjusted + dotDuration - 0.0001]', [], 1);
    
    temp_phi = [full_positions_repeated(:, 1)'; full_positions_repeated(:, 1)'];
    stepPositionsPhi = temp_phi(:);
    
    temp_theta = [full_positions_repeated(:, 2)'; full_positions_repeated(:, 2)'];
    stepPositionsTheta = temp_theta(:);
    
    % --- CALCULATE WINDOW EDGES FOR PLOTTING ---
    middle_segment_duration = 1;
    start_offset = (dotDuration - middle_segment_duration) / 2;
    % Appearance times are already shifted (start at dot1_start_time)
    plotted_appearance_times = appearanceTimes_adjusted; 
    
    window_starts = plotted_appearance_times + start_offset;
    window_ends = window_starts + middle_segment_duration;
    
    % --- PLOTTING (Main Figure: 4 Tiles) ---
    figure; 
    tiled_fig_handle = tiledlayout(4, 1); 
    % Update Title with new dot duration and t=0 alignment
    tiled_fig_handle.Title.String = sprintf("Transformed Data (Dot 1 Start: %.3fs, %.2fs Dot Duration)", dot1_start_time, dotDuration); 
    tiled_fig_handle.Title.FontWeight = "bold"; 
    
    % Tile 1: Gaze Angle by Time
    nexttile; 
    title("Gaze Angle by Time (Relative to Corrected T=0) [High Confidence Points Only]"); 
    hold on; 
    
    % Get y-limits before plotting data
    y_min = min(min(gaze_angles_high_confidence_centered)) - 2;
    y_max = max(max(gaze_angles_high_confidence_centered)) + 2;

    % --- ADD TRANSPARENT GREY PATCHES (BANDS) ---
    grayColor = [0.8 0.8 0.8]; % Light Gray
    alpha = 0.3; % Transparency (30%)
    
    for i = 1:length(window_starts)
        % Define patch vertices (bottom-left, top-left, top-right, bottom-right)
        x = [window_starts(i), window_starts(i), window_ends(i), window_ends(i)];
        y = [y_min, y_max, y_max, y_min];
        
        % Plot patch and set its properties
        patch(x, y, grayColor, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    % Re-set y limits since patch creation can change them
    ylim([y_min, y_max]); 
    
    % Plot data points (must be plotted AFTER patches for correct layering)
    plot(pupil_t_high_confidence(:, 1), gaze_angles_high_confidence_centered(:, 1), 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Phi'); 
    plot(pupil_t_high_confidence(:, 1), gaze_angles_high_confidence_centered(:, 2), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show; 
    xlim([-5, experimentEndTime]); % Synchronized X-axis
    hold off; 
    
    % Tile 2: Target Position Timecourse (STEP FUNCTION)
    nexttile;
    title(sprintf("Target Position (Step Function: %.2fs Duration)", dotDuration));
    hold on;
    stairs(stepTimes, stepPositionsPhi, 'r', 'LineWidth', 1.5, 'DisplayName', 'Phi'); 
    stairs(stepTimes, stepPositionsTheta, 'b', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    legend('Phi', 'Theta');
    % XLim starts at 0 as requested for Tile 2
    xlim([0, experimentEndTime]); 
    hold off;
    
    % Tile 3: Avg. Gaze Angle in window
    nexttile
    title('Avg. Gaze Angle in middle 1s window');
    hold on;
    plot(midpoints, avg_gaze_angles(:, 1), 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Phi');
    plot(midpoints(:, 1), avg_gaze_angles(:, 2), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show;
    xlim([-5, experimentEndTime]); % Synchronized X-axis
    hold off;
    
    % Tile 4: Pupil Detection Confidence
    nexttile;
    
    title_str = sprintf("Pupil Detection Confidence (Filter $\\geq$ %.2f, %.2f%% tossed from averaging)", confidence_cutoff, tossed_percentage);
    title(title_str);
    
    hold on; 
    % Plot all confidence points (from the selected time segment)
    plot(pupil_t_corrected_all(:, 1), confidence_measures_all(:, 1), 'kx', 'MarkerSize', 4, 'DisplayName', 'Confidence (all points)');
    plot(pupil_t_corrected_all(:, 1), confidence_measures_all(:, 2), 'gx', 'MarkerSize', 4, 'DisplayName', 'Model Confidence (all points)');

    % Add horizontal line at confidence cutoff
    yline(confidence_cutoff, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Cutoff (%.2f)', confidence_cutoff));
    
    xlabel("Time [s]"); 
    ylabel("Confidence"); 
    ylim([0, 1]);
    legend show; 
    xlim([-5, experimentEndTime]); % Synchronized X-axis
    hold off; 
    
    return ;
end 

% -------------------------------------------------------------------------

% NEW FUNCTION: Creates the separate Gaze Calibration figure
function create_calibration_figure(avg_gaze_angles, positionsRepeated, dotDuration)
    figure; % Create new figure
    
    title('Gaze Calibration Plot (\Phi vs. \Theta)');
    hold on;

    % --- Plot Target Locations (o) ---
    % Target list is positionsRepeated (ALL 72 elements)
    plot(positionsRepeated(:, 1), positionsRepeated(:, 2), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Target Position (o)');

    % --- Plot Recovered Gaze Angles (x) with different colors per repeat ---
    
    % All 72 points are kept. Each repeat has 24 points.
    nDotsPerRep_full = 24; % Base number of dots per repeat
    
    % Indices for the 72 data points:
    idx_rep1_end = nDotsPerRep_full; % 24 points (1-24)
    idx_rep2_end = idx_rep1_end + nDotsPerRep_full; % 48 points (1-48)
    
    dots_rep1 = avg_gaze_angles(1 : idx_rep1_end, :); 
    dots_rep2 = avg_gaze_angles(idx_rep1_end + 1 : idx_rep2_end, :); 
    dots_rep3 = avg_gaze_angles(idx_rep2_end + 1 : end, :); 
    
    % Shades of Red for the three repeats (R, G, B)
    color_rep1 = [1.0, 0.5, 0.5]; % Light Red
    color_rep2 = [1.0, 0.0, 0.0]; % Medium Red
    color_rep3 = [0.7, 0.0, 0.0]; % Dark Red

    plot(dots_rep1(:, 1), dots_rep1(:, 2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', color_rep1, 'DisplayName', 'Recovered Gaze (Rep 1)');
    plot(dots_rep2(:, 1), dots_rep2(:, 2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', color_rep2, 'DisplayName', 'Recovered Gaze (Rep 2)');
    plot(dots_rep3(:, 1), dots_rep3(:, 2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', color_rep3, 'DisplayName', 'Recovered Gaze (Rep 3)');

    xlabel('Phi Angle [deg]');
    ylabel('Theta Angle [deg]');
    
    % Ensure plot is centered at 0,0 and axes limits encompass all data
    all_x = [positionsRepeated(:, 1); avg_gaze_angles(:, 1)];
    all_y = [positionsRepeated(:, 2); avg_gaze_angles(:, 2)];
    
    % Find the maximum absolute value to ensure symmetry around 0
    max_range = max(abs([all_x; all_y])) + 5; 
    
    xlim([-max_range, max_range]);
    ylim([-max_range, max_range]);
    axis equal; % Ensures correct angular representation
    grid on;
    legend show;
    hold off;
end

% -------------------------------------------------------------------------

% Local function for calculating dot appearance times 
function appearanceTimes = getDotAppearanceTimes_updated(dotDuration, nDotsPerRep, repeats, interRepDelay)
    totalDots = nDotsPerRep * repeats;
    appearanceTimes = zeros(totalDots, 1);
    time = 0;  
    for rep = 1:repeats
        for ii = 1:nDotsPerRep
            appearanceTimes((rep-1)*nDotsPerRep + ii) = time;  
            time = time + dotDuration;                         
        end
        time = time + interRepDelay;                          
    end
end

% -------------------------------------------------------------------------

% Local function for averaging gaze angles 
function [midpoints, avg_gaze_angles] = average_gaze_middle_segment_updated(gaze_angles_transformed, pupil_t_corrected, appearanceTimes_adjusted, dotDuration)
    n_targets = length(appearanceTimes_adjusted);
    avg_gaze_angles = nan(n_targets, 2);
    midpoints = nan(n_targets, 1);
    
    % Use a 1s window, centered within the dotDuration 
    middle_segment_duration = 1;
    % The start offset is calculated based on the dotDuration
    start_offset = (dotDuration - middle_segment_duration) / 2;
    
    for ii = 1:n_targets
        % Segment start/end are calculated relative to appearanceTimes_adjusted(ii), 
        % which is the corrected plot time.
        segment_start = appearanceTimes_adjusted(ii) + start_offset;
        segment_end = segment_start + middle_segment_duration;
        
        % Find indices of gaze data within this segment
        % pupil_t_corrected is also on the corrected time scale.
        idx_in_segment = pupil_t_corrected >= segment_start & pupil_t_corrected < segment_end;
        
        % The midpoint is on the corrected time scale
        midpoints(ii) = (segment_start + segment_end) / 2;
        
        if any(idx_in_segment)
            avg_gaze_angles(ii, :) = mean(gaze_angles_transformed(idx_in_segment, :), 1);
        else
            avg_gaze_angles(ii, :) = [NaN, NaN];
        end 
    end
end