function [frame_list] = findGazeFrames(start_time, gaze_targets_deg, target_dur_s)
%FINDGAZEFRAMES Automated process to identify representative frames from gaze calibration.
%   Uses the MEDIAN of the Xp/Yp points to find a robust center 
%   frame for each target, based on a 120 fps rate.
%
%   NOTE: When a target fails to yield a valid frame (returns NaN), 
%   a detailed diagnostic is printed to the command window.
% Example
%{ 
    gaze_targets_deg = [ ...
            0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
            0, 15; 0, -15; -15, 0; 15, 0;...
            -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
            0, 10; 0, -7.5; -7.5, 0; 7.5, 0];
    findGazeFrames([1, 23, 683], gaze_targets_deg, 3.267)
%}
    arguments
        start_time (1,1) double      % time of first dot in [minute, second, millisecond] format. Human observer should determine this with IINA for now.
        gaze_targets_deg (:, 2) double         % N x 2 array of target positions [phi, theta] (Used for N, but not position values)
        target_dur_s (1,1) double = 3.43       % Duration (s) each dot was presented (optional, default 3.43s)
    end
    % Hardcoded Parameters
    FPS = 120;
    analysis_window_duration = 1.5; % seconds

    % CONFIDENCE PARAMETERS
    confidence_cutoff = 0.7; 
    min_perimeter_points = 8;

    % convert time to frame number
    start_frame = estimateFrameFromTime(FPS, start_time);
    
    % --- 1. DATA LOADING AND PREPARATION ---
    data_file_path = '/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/sam_gazecal_106.mat';
    
    % UPDATED DATA LOADING based on user feedback
    pupil_features_struct = load(data_file_path, 'perimeter');
    pupil_features = pupil_features_struct.perimeter.data;
    
    % Call flatten_features to filter frames and extract valid data
    [median_pupil_centers, ~, ~, frame_idx] = flatten_features(pupil_features, confidence_cutoff, min_perimeter_points);
    
    % --- 2. CONVERT TIME TO FRAMES ---
    target_dur_frames = round(target_dur_s * FPS);
    analysis_window_frames = round(analysis_window_duration * FPS);
    
    % --- 3. TIMING AND WINDOW SETUP (in Frames) ---
    N = size(gaze_targets_deg, 1); 
    frame_list = nan(N, 1); 
    
    % Determine the start frame for each target
    dot_start_frames = start_frame + (0:N-1)' * target_dur_frames;
    
    % Calculate the start and end of the 1.5s analysis window for each target (in frames)
    center_offset_frames = round(target_dur_frames / 2);
    half_window_frames = round(analysis_window_frames / 2);
    window_starts_frame = dot_start_frames + center_offset_frames - half_window_frames;
    window_ends_frame = window_starts_frame + analysis_window_frames - 1; 
    
    % --- 4. FRAME SELECTION LOOP ---
    for i = 1:N
        f_start = window_starts_frame(i);
        f_end = window_ends_frame(i);
        
        % a) Segment: Find all valid (pre-filtered) frames within the analysis window
        idx_window = (frame_idx >= f_start) & (frame_idx <= f_end); 
        
        centers_window = median_pupil_centers(idx_window, :);
        frame_idx_window = frame_idx(idx_window);
        
        % Skip if no valid frames found in the window 
        if isempty(centers_window)
            warning('No valid detected frames found for target %d in the window [Frame %d to %d].', i, f_start, f_end);
            
            % *** DIAGNOSTIC OUTPUT HERE ***
            diagnoseFailure(pupil_features, f_start, f_end, confidence_cutoff, min_perimeter_points);
            
            continue;
        end
        % b) Find the Median (X_median and Y_median) across all frames in the window
        median_x_across_frames = median(centers_window(:, 1));
        median_y_across_frames = median(centers_window(:, 2));
        
        % c) Find the Frame Closest to the Median
        median_center_across_frames = [median_x_across_frames, median_y_across_frames];
        distances = hypot(centers_window(:, 1) - median_center_across_frames(1), centers_window(:, 2) - median_center_across_frames(2));
        [~, min_dist_sub_idx] = min(distances);
        selected_frame_index = frame_idx_window(min_dist_sub_idx);
        
        frame_list(i) = selected_frame_index;
    end
end
% -------------------------------------------------------------------------
% --- LOCAL FUNCTIONS ---
function [median_pupil_centers, frame_confidence, pupil_t, frame_idx] = flatten_features(pupil_features, confidence_cutoff, min_perimeter_points)
    % UPDATED: Now assumes pupil_features is the direct cell array of perimeter structs.
    
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
% Helper function to print detailed reasons for why frames failed the confidence check.
    
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