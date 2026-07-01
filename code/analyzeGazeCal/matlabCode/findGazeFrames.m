function [frame_list] = findGazeFrames(start_time, gaze_targets_deg, perimeterFile, target_dur_s, onset_dur_s, confidence_cutoff, min_perimeter_points)
% Identify representative pupil perimeter frames for each gaze calibration target
%
% Syntax:
%   frame_list = findGazeFrames(start_time, gaze_targets_deg, perimeterFile)
%   frame_list = findGazeFrames(start_time, gaze_targets_deg, perimeterFile, target_dur_s, onset_dur_s, confidence_cutoff, min_perimeter_points)
%
% Description:
%   For each gaze calibration target, this function identifies the single
%   video frame whose pupil center is closest to the median pupil position
%   within a 1.5-second analysis window centered on the expected fixation
%   interval. Uses a 120 fps rate assumption. Handles skipped targets by
%   detecting large inter-frame gaps and adjusting subsequent timing
%   windows. Prints detailed diagnostics when a target fails to yield a
%   valid frame.
%
% Inputs:
%   start_time            - 1x3 double. Time of the first fixation dot in
%                           [minutes, seconds, milliseconds] format.
%   gaze_targets_deg      - Nx2 double. Target positions in degrees of
%                           visual angle [azimuth, elevation]. N determines
%                           the number of targets to search for.
%   perimeterFile         - Char. Full path to the .mat file containing
%                           the pupil perimeter data struct.
%   target_dur_s          - Scalar double (default: 3.43). Duration in
%                           seconds each dot was presented.
%   onset_dur_s           - Scalar double (default: 2.6). Duration of the
%                           first fixation, which may be shorter due to
%                           eye opening.
%   confidence_cutoff     - Scalar double (default: 0.7). Minimum
%                           confidence value for perimeter points to be
%                           considered valid.
%   min_perimeter_points  - Scalar double (default: 8). Minimum number of
%                           high-confidence perimeter points required for
%                           a frame to be a valid candidate.
%
% Outputs:
%   frame_list            - Nx1 double. Frame indices for each target. NaN
%                           entries indicate targets where no valid frame
%                           was found.
%
% Examples:
%{
    start_time = [1, 25, 400];
    gaze_targets_deg = [0 0; -15 15; -15 -15; 15 15];
    perimeterFile = '/path/to/perimeter.mat';
    frame_list = findGazeFrames(start_time, gaze_targets_deg, perimeterFile, 3.34, 2.6, 0.8, 8);
%}
    arguments
        start_time (1,3) double      % time of first dot in [minute, second, millisecond] format. Human observer should determine this with IINA for now.
        gaze_targets_deg (:, 2) double         % N x 2 array of target positions [phi, theta] (Used for N, but not position values)
        perimeterFile char
        target_dur_s (1,1) double = 3.43       % Duration (s) each dot was presented (optional, default 3.43s)
        onset_dur_s = 2.6 % how much shorter is the first fixation than intended?
        confidence_cutoff = 0.7; %confidencce of perimeter points 
        min_perimeter_points = 8;% number of points at above confidence required to be a "good frame"
    end
    % Hardcoded Parameters
    FPS = 120;
    analysis_window_duration = 1.5; % seconds

    % ADJUSTMENT PARAMETER: How many dot durations (frames) of "slack" is allowed.
    % If the difference between selected frames is greater than this, we assume a skipped dot.
    % We use 1.5 * target_dur_frames as a generous threshold.
    FRAME_DIFFERENCE_THRESHOLD_FACTOR = 1.5;
    
    % convert time to frame number
    start_frame = time2frame(start_time, FPS);
    onset_dur_frames = round(onset_dur_s * FPS); 

    % --- 1. DATA LOADING AND PREPARATION ---
    pupil_features_struct = load(perimeterFile, 'perimeter');
    pupil_features = pupil_features_struct.perimeter.data;
    
    % Call flatten_features to filter frames and extract valid data
    [median_pupil_centers, ~, ~, frame_idx] = flatten_features(pupil_features, confidence_cutoff, min_perimeter_points);
    
    % --- 2. CONVERT TIME TO FRAMES ---
    target_dur_frames = round(target_dur_s * FPS);
    analysis_window_frames = round(analysis_window_duration * FPS);
    
    % --- 3. TIMING AND WINDOW SETUP (in Frames) ---
    N = size(gaze_targets_deg, 1);
    frame_list = nan(N, 1);
    % Create per-target durations
    target_dur_each = repmat(target_dur_frames, N, 1); % default durations
    target_dur_each(1) = onset_dur_frames; % shorten first target
    % Determine the start frame for each target
    dot_start_frames = start_frame + [0; cumsum(target_dur_each(1:end-1))];
    % Calculate the start and end of the 1.5s analysis window for each target (in frames)
    center_offset_frames = round(target_dur_each / 2);
    half_window_frames = round(analysis_window_frames / 2);
    
    % **The dot_start_frames and related window variables will be dynamically updated**
    window_starts_frame = dot_start_frames + center_offset_frames - half_window_frames;
    window_ends_frame = window_starts_frame + analysis_window_frames - 1; 
    
    % Set the frame difference threshold
    frame_diff_threshold = round(target_dur_frames * FRAME_DIFFERENCE_THRESHOLD_FACTOR);
    
    % Initialize the frame index for the previously selected target
    last_selected_frame = start_frame; % Use the start of the first dot as the previous reference
    
    % --- 4. FRAME SELECTION LOOP (with skipped target correction) ---
    for i = 1:N
        f_start = window_starts_frame(i);
        f_end = window_ends_frame(i);
        
        % a) Segment: Find all valid (pre-filtered) frames within the analysis window
        idx_window = (frame_idx >= f_start) & (frame_idx <= f_end); 
        
        centers_window = median_pupil_centers(idx_window, :);
        frame_idx_window = frame_idx(idx_window);
        
        % Check if a valid frame was found (the rest of the code is unchanged from here)
        if isempty(centers_window)
            warning('No valid detected frames found for target %d in the window [Frame %d to %d].', i, f_start, f_end);
            
            % *** DIAGNOSTIC OUTPUT HERE ***
            diagnoseFailure(pupil_features, f_start, f_end, confidence_cutoff, min_perimeter_points);
            
            % If a target is missed, we cannot reliably update 'last_selected_frame'
            % and the skip correction should be handled by the next *successful* selection.
            continue; 
        end
        
        % b) Find the Frame Closest to the Median
        median_x_across_frames = median(centers_window(:, 1));
        median_y_across_frames = median(centers_window(:, 2));
        median_center_across_frames = [median_x_across_frames, median_y_across_frames];
        distances = hypot(centers_window(:, 1) - median_center_across_frames(1), centers_window(:, 2) - median_center_across_frames(2));
        [~, min_dist_sub_idx] = min(distances);
        selected_frame_index = frame_idx_window(min_dist_sub_idx);
        
        % --- NEW: SKIPPED TARGET CORRECTION LOGIC ---
        frame_diff = selected_frame_index - last_selected_frame;
        
        % The expected duration is target_dur_frames, unless it's the first target, 
        % which is shortened by onset_delay_frames.
        expected_frame_dur = target_dur_each(i);
        
        % Check if the frame difference is significantly larger than expected
        if frame_diff > (expected_frame_dur + frame_diff_threshold)
            
            % Calculate how many full targets were skipped
            num_targets_skipped = floor( (frame_diff - expected_frame_dur) / target_dur_frames);
            
            if num_targets_skipped > 0
                warning('Large time gap detected between selected frames (Target %d). Frame jump of %d (Expected: ~%d). Assuming %d target(s) were skipped.', ...
                        i, frame_diff, expected_frame_dur, num_targets_skipped);
                
                % Calculate the total adjustment needed for all subsequent targets
                adjustment_frames = num_targets_skipped * target_dur_frames;
                
                % Apply the adjustment to the start frames for all *remaining* targets (including the current one)
                dot_start_frames(i:end) = dot_start_frames(i:end) - adjustment_frames;
                
                % Re-calculate the windows based on the adjusted dot start times
                window_starts_frame = dot_start_frames + center_offset_frames - half_window_frames;
                window_ends_frame = window_starts_frame + analysis_window_frames - 1; 
                
                % Since we applied the correction *on* this successful frame selection, 
                % we need to re-run the loop for the current index 'i' to find the 
                % *correct* frame within the *new* (earlier) analysis window.
                % The 'continue' will jump to the next iteration, which will be the 
                % current index 'i' again since 'i' increments at the loop start.
                
                % To re-run the current index, we must decrement 'i' here.
                i = i - 1; 
                continue; % Restart current target search with the adjusted window
            end
        end
        
        % If loop successfully completes or no significant skip was found:
        frame_list(i) = selected_frame_index;
        last_selected_frame = selected_frame_index; % Update for the next comparison
    end
end
% -------------------------------------------------------------------------
% --- LOCAL FUNCTIONS ---
function [median_pupil_centers, frame_confidence, pupil_t, frame_idx] = flatten_features(pupil_features, confidence_cutoff, min_perimeter_points)
% Internal helper to flatten features.
%
% Syntax:
%   median_pupil_centers, frame_confidence, pupil_t, frame_idx = flatten_features(pupil_features, confidence_cutoff, min_perimeter_points)
%
% Description:
%   This local helper function internal helper to flatten features within its parent workflow.
% Inputs:
%   pupil_features           - Input used by the function.
%   confidence_cutoff        - Input used by the function.
%   min_perimeter_points     - Input used by the function.
%
% Outputs:
%   median_pupil_centers     - Output produced by the function.
%   frame_confidence         - Output produced by the function.
%   pupil_t                  - Output produced by the function.
%   frame_idx                - Output produced by the function.
%
% Examples:
%{
    % See findGazeFrames.m for usage context.
%}

    if ~iscell(pupil_features) 
        return;
    end

    N_frames = numel(pupil_features);

    median_pupil_centers = nan(N_frames, 2);
    frame_confidence = nan(N_frames, 1); 

    pupil_t = nan(N_frames, 1); 
    frame_idx = nan(N_frames, 1);

    % FIELD NAMES
    Xp_field = 'Xp';
    Yp_field = 'Yp';
    Confidence_field = 'confidence';

    for ii = 1:N_frames

        frame_idx(ii) = ii; 

        try
            % SIMPLIFIED ACCESS: pupil_features{ii} is now the data struct.
            frame_data_struct = pupil_features{ii}; 
        catch
            continue;
        end

        if isempty(frame_data_struct) || ~isstruct(frame_data_struct)
             continue; 
        end

        % Check for required fields
        if ~isfield(frame_data_struct, Xp_field) || ~isfield(frame_data_struct, Yp_field) || ~isfield(frame_data_struct, Confidence_field)
             continue;
        end

        % --- CONFIDENCE CHECK LOGIC ---
        perimeter_confidences = frame_data_struct.(Confidence_field);

        num_high_confidence_points = sum(perimeter_confidences > confidence_cutoff);

        if num_high_confidence_points < min_perimeter_points
            continue;
        end
        % --- END CONFIDENCE CHECK LOGIC ---

        % a) Calculate the Median Center (ONLY if confidence check passed)
        median_x = median(frame_data_struct.(Xp_field));
        median_y = median(frame_data_struct.(Yp_field));

        median_pupil_centers(ii, :) = [median_x, median_y];

        % b) Extract frame-level metadata (using perimeter confidence since there's no top-level field)
        frame_confidence(ii) = median(perimeter_confidences); 
    end

    % Clean up NaN rows (frames that failed to load OR failed the confidence check)
    valid_idx = ~all(isnan(median_pupil_centers), 2);
    median_pupil_centers = median_pupil_centers(valid_idx, :);

    frame_confidence = frame_confidence(valid_idx); 
    pupil_t = pupil_t(valid_idx); 
    frame_idx = frame_idx(valid_idx);
end

% -------------------------------------------------------------------------
function diagnoseFailure(pupil_features, f_start, f_end, confidence_cutoff, min_perimeter_points)
% Internal helper to diagnose failure.
%
% Syntax:
%   diagnoseFailure(pupil_features, f_start, f_end, confidence_cutoff, min_perimeter_points)
%
% Description:
%   This local helper function internal helper to diagnose failure within its parent workflow.
% Inputs:
%   pupil_features           - Input used by the function.
%   f_start                  - Input used by the function.
%   f_end                    - Input used by the function.
%   confidence_cutoff        - Input used by the function.
%   min_perimeter_points     - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See findGazeFrames.m for usage context.
%}

    window_indices = f_start:f_end;

    % Get the size of the total data set for bounds checking
    max_frame_idx = numel(pupil_features);

    fprintf('--- Diagnosis for Target Window [Frame %d to %d]: ---\n', f_start, f_end);
    fprintf('Check: Requires > %d points with confidence > %.2f\n', min_perimeter_points, confidence_cutoff);

    % FIELD NAME
    Confidence_field = 'confidence';

    for ii = window_indices
        if ii > max_frame_idx
             fprintf('Frame %d: Outside of available data range (Max frame in log is %d).\n', ii, max_frame_idx);
             continue;
        end

        try
            % SIMPLIFIED ACCESS: pupil_features{ii} is now the data struct.
            frame_data_struct = pupil_features{ii}; 
        catch
             fprintf('Frame %d: DATA LOAD FAILED (e.g., frame struct is missing or corrupted).\n', ii);
             continue;
        end

        if isempty(frame_data_struct) || ~isstruct(frame_data_struct)
             fprintf('Frame %d: DATA LOAD FAILED (empty/not a struct).\n', ii);
             continue;
        end

        if ~isfield(frame_data_struct, Confidence_field)
             fprintf('Frame %d: Missing ''confidence'' field. Cannot perform confidence check.\n', ii);
             continue;
        end

        perimeter_confidences = frame_data_struct.(Confidence_field);
        num_points_total = numel(perimeter_confidences);
        num_high_confidence_points = sum(perimeter_confidences > confidence_cutoff);

        % The failure reason
        if num_high_confidence_points < min_perimeter_points
            if num_high_confidence_points == 0
                 fprintf('Frame %d: FAILED. All %d perimeter points were below the %.2f cutoff.\n', ...
                         ii, num_points_total, confidence_cutoff);
            else
                 fprintf('Frame %d: FAILED. Only %d/%d points had confidence > %.2f (Required: %d).\n', ...
                         ii, num_high_confidence_points, num_points_total, confidence_cutoff, min_perimeter_points);
            end
        else
            fprintf('Frame %d: PASSED confidence check.\n', ii);
        end
    end
    fprintf('---------------------------------------------------\n');
end

