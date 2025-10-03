function [frame_list] = findGazeFrames(t_start_s, gaze_targets_deg, target_dur_s)
%FINDGAZEFRAMES Automated process to identify representative frames from gaze calibration.
%   Uses the MEDIAN of the 16 Xp/Yp perimeter points to find a robust center 
%   frame for each target, based on a 120 fps rate.
% Example
%{ 
    gaze_targets_deg = [0, 0; 0, 20; 0, -20; -20, 0; 20, 0;...
        0, 0; 0, 20; 0, -20; -20, 0; 20, 0;];
    findGazeFrames(88.5370, gaze_targets_deg)
%}

    arguments
        t_start_s (1,1) double                 % Time (s) the experiment task (first dot ONSET) began.
        gaze_targets_deg (:, 2) double         % N x 2 array of target positions [phi, theta] (Used for N, but not position values)
        target_dur_s (1,1) double = 3.02       % Duration (s) each dot was presented (optional, default 3.02s)
    end

    % Hardcoded Parameters
    FPS = 120;
    analysis_window_duration = 1.5; % seconds
    confidence_cutoff = 0.8; % Applies to the FRAME-LEVEL confidence

    % --- 1. DATA LOADING AND PREPARATION ---

    % UPDATED FILE PATH
    data_file_path = '/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/SM_gaze_cal_pupil_features101_dark.mat';
    
    % Load the data structure. Access the field that contains the cell array.
    pupil_features_struct = load(data_file_path, 'pupil_features');
    pupil_features = pupil_features_struct.pupil_features;

    % Local function now returns MEDIAN PIXEL CENTER and the timestamp
    [median_pupil_centers, confidence_measures, pupil_t, frame_idx] = flatten_features(pupil_features);
    
    % --- 2. CONVERT TIME TO FRAMES ---
    
    % Convert the task start time (t_start_s) to the starting frame number (1-based)
    % This assumes the video starts at time t=0
    start_frame = round(t_start_s * FPS) + 1;
    
    % Convert duration parameters to frames
    target_dur_frames = round(target_dur_s * FPS);
    analysis_window_frames = round(analysis_window_duration * FPS);
    
    % --- 3. TIMING AND WINDOW SETUP (in Frames) ---

    N = size(gaze_targets_deg, 1); % N is now 10
    frame_list = nan(N, 1); % The list of output frame indices (original 1-based index)
    
    % Determine the start frame for each target
    % Target 1 starts at start_frame
    dot_start_frames = start_frame + (0:N-1)' * target_dur_frames;
    
    % Calculate the start and end of the 1.5s analysis window for each target (in frames)
    center_offset_frames = round(target_dur_frames / 2);
    half_window_frames = round(analysis_window_frames / 2);

    window_starts_frame = dot_start_frames + center_offset_frames - half_window_frames;
    window_ends_frame = window_starts_frame + analysis_window_frames;

    % --- 4. FRAME SELECTION LOOP ---

    for i = 1:N
        f_start = window_starts_frame(i);
        f_end = window_ends_frame(i);

        % a) Segment: Find all frames within the analysis window and confidence cutoff
        % NOTE: We filter on the 'frame_idx' which is the original 1-based index (frame number)
        idx_window = (frame_idx >= f_start) & (frame_idx < f_end) & ...
                     (confidence_measures > confidence_cutoff); % Filter by the single frame-level confidence
        
        centers_window = median_pupil_centers(idx_window, :);
        frame_idx_window = frame_idx(idx_window);

        % Skip if no valid frames found in the window
        if isempty(centers_window)
            warning('No valid frames found for target %d in the window [Frame %d to %d].', i, f_start, f_end);
            continue;
        end

        % b) Find the Median (X_median and Y_median) across all frames in the window
        % This is the core logic you requested: median position across the set of frames
        median_x_across_frames = median(centers_window(:, 1));
        median_y_across_frames = median(centers_window(:, 2));
        
        % c) Find the Frame Closest to the Median
        
        % Calculate Euclidean distance from each frame's median center to the overall median center
        median_center_across_frames = [median_x_across_frames, median_y_across_frames];
        
        % Distance calculation using hypot
        distances = hypot(centers_window(:, 1) - median_center_across_frames(1), centers_window(:, 2) - median_center_across_frames(2));
        
        % Find the index within the *centers_window* array that has the minimum distance
        [~, min_dist_sub_idx] = min(distances);
        
        % The actual frame index (frame number) to return
        selected_frame_index = frame_idx_window(min_dist_sub_idx);
        
        frame_list(i) = selected_frame_index;
    end
end

% -------------------------------------------------------------------------

% --- LOCAL FUNCTIONS ---

% Local function to extract the desired features from the list
function [median_pupil_centers, frame_confidence, pupil_t, frame_idx] = flatten_features(pupil_features)
    % Calculates the median (X, Y) center from the 16 Xp/Yp perimeter points for each frame.
    
    N_frames = numel(pupil_features);
    median_pupil_centers = nan(N_frames, 2);
    frame_confidence = nan(N_frames, 1);
    
    % Since the timestamp field is missing, we will use the frame index (ii)
    % for all time-based calculations, which is handled by frame_idx.
    pupil_t = nan(N_frames, 1); 
    frame_idx = nan(N_frames, 1);

    % Iterate over the pupil features structs
    for ii = 1:N_frames
        
        % The index ii is the frame number (1-based)
        frame_idx(ii) = ii; 

        frame_data_struct = pupil_features{ii};

        % Robustly check if the element is a struct AND non-empty before dot indexing
        if isempty(frame_data_struct) || ~isstruct(frame_data_struct) || isempty(frame_data_struct)
             continue;
        end
        
        % Now that we know it's a non-empty struct, we can safely check its fields.
        if ~isfield(frame_data_struct, 'data') || isempty(frame_data_struct.data)
             continue;
        end
        
        % Access the actual data (which is a struct array inside a cell)
        data = frame_data_struct.data{1};
        
        % Ensure data fields exist
        if ~isfield(data, 'Xp') || ~isfield(data, 'Yp')
             continue;
        end

        % a) Calculate the Median Center (our equivalent of gaze angles)
        % Using the MEDIAN of the 16 perimeter points for X and Y separately
        median_x = median(data.Xp);
        median_y = median(data.Yp);
        
        median_pupil_centers(ii, :) = [median_x, median_y];
        
        % b) Extract frame-level metadata
        
        % REMOVED: pupil_t(ii) = frame_data_struct.timestamp;
        % The frame index 'ii' will be used for time instead.
        
        % The confidence field is now outside the 'data' struct, directly in the frame struct
        if isfield(frame_data_struct, 'confidence')
             frame_confidence(ii) = frame_data_struct.confidence;
        else
             % If the frame-level confidence is not available, use the median of the perimeter confidences
             frame_confidence(ii) = median(data.confidence); 
        end
    end
    
    % Clean up NaN rows (frames where no detection was made)
    valid_idx = ~isnan(frame_confidence);
    median_pupil_centers = median_pupil_centers(valid_idx, :);
    frame_confidence = frame_confidence(valid_idx);
    
    % The pupil_t array is now just NaNs, but we filter out the corresponding frames
    % to keep the array lengths consistent with the centers and confidences.
    pupil_t = pupil_t(valid_idx); 
    frame_idx = frame_idx(valid_idx);
end