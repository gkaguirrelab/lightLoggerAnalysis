function analyze_gaze_calibration()

    arguments
    
    end 

    % Load in the extracted pupil features
    dropBoxBaseDir=getpref('combiExperiments','dropboxBaseDir');
    pupil_features = load(fullfile(dropBoxBaseDir, '/FLIC_data/lightLogger/gaze_cal_pupil_features')).pupil_features;

    % Load in the target positions in degrees (will make this dynamic in the future)
    deg_positions = [ ...
        0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
        0, 20;   0, -20;   -20,   0;   20, 0; ...
        -15, 15;  15, 15;   -15, -15;   15, -15; ...

        -10,  10;   -10, -10;    10,  10;    10, -10; ...
        0,  10;    0, -10;   -10,   0;     10,   0; ...
        -5,   5;   5,   5;   -5,  -5;     5,  -5;     0,   0];


    % Retrieve the features we want to analyze 
    [gaze_angles, confidence_measures, pupil_t] = flatten_features(pupil_features);
    
    % Ensure we generated the correct number of points 
    assert(size(gaze_angles, 1) == size(pupil_t, 1), "Incorrect number of t values for datapoints"); 

    % Display the features in their raw state 
    %display_raw_results(gaze_angles, confidence_measures, pupil_t, deg_positions);

    % Display the data after making desired transformations 
    % (splicing out extra parts of video, mean center, remove low confidence points)
    display_adjusted_results(gaze_angles, confidence_measures, pupil_t, deg_positions); 


    return ; 

end 


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

% Local function to display raw results before adjustments
function display_raw_results(gaze_angles, confidence_measures, pupil_t, deg_positions)
    % First, make a plot of the raw features. No clipping, no confidence threshold, no mean subtraction 
    figure; 
    tiled_fig_handle = tiledlayout(1, 3); 
    tiled_fig_handle.Title.String = "Transformed Data"; 
    tiled_fig_handle.Title.FontWeight = "bold"; 

    % Move to the first tile
    % and plot the gaze angles over time and the degree positions
    % of the targets 
    nexttile; 
    title("Gaze Angle by Time"); 
    hold on; 
    plot(pupil_t(:, 1), gaze_angles(:, 1), '-o', 'DisplayName', 'Phi'); 
    plot(pupil_t(:, 1), gaze_angles(:, 2), '-o', 'DisplayName', 'Theta');
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show; 
    hold off; 

    
    % Move to the second tile
    % and plot the gaze angles over time and the degree positions
    % of the targets 
    nexttile; 
    title("Target Position"); 
    hold on; 
    plot(deg_positions(:, 1), 'x', 'DisplayName', 'Phi'); 
    plot(deg_positions(:, 2), 'x', 'DisplayName', 'Theta'); 
    xlabel("Target Number"); 
    ylabel("Angle [deg]"); 
    legend show; 
    hold off;

    % Move to the last tile and plot the confidence of fits 
    nexttile; 
    title("Pupil Detection Confidence"); 
    hold on; 
    plot(pupil_t(:, 1), confidence_measures(:, 1), 'x', 'DisplayName', 'Confidence');
    plot(pupil_t(:, 1), confidence_measures(:, 2), 'x', 'DisplayName', 'Model Confidence');
    xlabel("Time [s]"); 
    ylabel("Confidence"); 
    ylim([0, 1]);
    legend show; 
    hold off; 

    return ;
end 

% Local function to display results after clipping to 
% only gaze calibration period
% mean centering
% removing low confidence points
function display_adjusted_results(gaze_angles, confidence_measures, pupil_t, deg_positions)
    % Now, plot the adjusted features after transformation has been applied
    figure; 
    tiled_fig_handle = tiledlayout(4, 1); 
    tiled_fig_handle.Title.String = "Transformed Data"; 
    tiled_fig_handle.Title.FontWeight = "bold"; 

    % 1. Trim out just the gaze calibration portion of the video 
    start_idx = 10800; 
    end_idx = size(gaze_angles, 1) - 3600; 
    pupil_t = pupil_t(start_idx:end_idx, :); 
    %in first recording, it looks like the first point starts around 95
    %seconds. May need to fiddle with this for other recordings, which is
    %not ideal. Find the point closest to 95 sec 
    [diff, idx] = min(abs(pupil_t - 95));
    pupil_t = pupil_t(idx:end, :); 
    pupil_t_corrected = pupil_t - pupil_t(1);
    gaze_angles = gaze_angles((start_idx +idx):end_idx, :); 
    confidence_measures = confidence_measures((start_idx + idx):end_idx, :); 

    % Initialize variable for transforming gaze angles
    gaze_angles_transformed = gaze_angles;

    % 2. Remove low confidence frames
    good_idx = (confidence_measures(:, 1) > 0.90) & (confidence_measures(:, 2) > 0.9);
    gaze_angles_transformed = gaze_angles_transformed(good_idx, :);
    confidence_measures = confidence_measures(good_idx, :);
    pupil_t_corrected = pupil_t_corrected(good_idx, :); 

    % 3. Subtract the mean theta from the theta vector, and the mean phi from the phi vector
    mean_phi_theta = mean(gaze_angles, 1); 
    gaze_angles_transformed(:, 1) = gaze_angles_transformed(:, 1) - mean_phi_theta(1);
    gaze_angles_transformed(:, 2) = gaze_angles_transformed(:, 2) - mean_phi_theta(2); 


    % Move to the first tile
    % and plot the gaze angles over time and the degree positions
    % of the targets 
    nexttile; 
    title("Gaze Angle by Frame"); 
    hold on; 
    plot(pupil_t_corrected(:, 1), gaze_angles_transformed(:, 1), 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Phi'); 
    plot(pupil_t_corrected(:, 1), gaze_angles_transformed(:, 2), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show; 
    hold off; 

    % Move to the second tile
    % and average gaze angle for 2-second time-bins
    nexttile
    title('Avg. Gaze Angle in window');
    hold on;
    window_duration = 2; % seconds, for now
    [midpoints, avg_gaze_angles] = average_gaze_middle_segment(gaze_angles_transformed, pupil_t_corrected, window_duration);
    
    plot(midpoints, avg_gaze_angles(:, 1), 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Phi');
    plot(midpoints(:, 1), avg_gaze_angles(:, 2), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show;
    hold off;

    % Move to the third tile
    % and plot the gaze angles over time and the degree positions
    % of the targets
    nexttile;
    % Repeat each position 3 times
    repeats = 3;
    positionsRepeated = repmat(deg_positions, repeats, 1);

    % Create time vector: each point lasts 1 second
    dotTime = 1; %seconds (as defined in runGazeCalibrationStimulus, but the timing was actually 2x this and then another 1 between repeats) 
    nDotsPerRep = size(deg_positions,1);
    appearanceTimes = getDotAppearanceTimes(dotTime, nDotsPerRep, repeats);

    % Plot
    title("Target Position");
    hold on;
    plot(appearanceTimes, positionsRepeated(:,1), 'ro-', 'LineWidth', 1.5); hold on;
    plot(appearanceTimes, positionsRepeated(:,2), 'bo-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    legend('Phi', 'Theta');
    hold off;

    % Move to the last tile and plot the confidence of fits
    nexttile;
    title("Pupil Detection Confidence"); 
    hold on; 
    plot(pupil_t_corrected(:, 1), confidence_measures(:, 1), 'x', 'DisplayName', 'Confidence');
    plot(pupil_t_corrected(:, 1), confidence_measures(:, 2), 'x', 'DisplayName', 'Model Confidence');
    xlabel("Time [s]"); 
    ylabel("Confidence"); 
    ylim([0, 1]);
    legend show; 
    hold off; 

    return ;
end 

function appearanceTimes = getDotAppearanceTimes(dotTime, nDots, repetitions)
% getDotAppearanceTimes returns a vector of timestamps when each dot appears on screen
%
% Inputs:
%   dotTime     - time in seconds for each wait period (pre-flip and post-flip)
%   nDots       - number of dots per repetition
%   repetitions - number of repetitions
%
% Output:
%   appearanceTimes - (nDots * repetitions) x 1 vector of appearance times in seconds

    totalDots = nDots * repetitions;
    appearanceTimes = zeros(totalDots, 1);

    time = 0;  % initialize time counter

    for rep = 1:repetitions
        for ii = 1:nDots
            time = time + dotTime;  % wait before flip (dot drawn but not visible yet)
            appearanceTimes((rep-1)*nDots + ii) = time;  % dot appears (flip time)
            time = time + dotTime;  % wait after flip (dot visible)
        end
        time = time + dotTime;  % wait between repetitions
    end
end

function [midpoints, avg_gaze_angles] = average_gaze_middle_segment(gaze_angles_transformed, pupil_t_corrected, window_duration)
    % average_gaze_middle_segment averages gaze angles over the middle segment of time windows
    %
    % Inputs:
    %   gaze_angles_transformed: Nx2 matrix of gaze angles [phi, theta]
    %   pupil_t_corrected: Nx1 vector of timestamps (seconds)
    %   window_duration: duration of the full time window in seconds (default 2s)
    %
    % Outputs:
    %   midpoints: Mx1 vector of midpoint times for each averaged window
    %   avg_gaze_angles: Mx2 matrix of average gaze angles [phi, theta] per middle segment
    
    if nargin < 3
        window_duration = 2; % default 2-second windows
    end

    start_time = pupil_t_corrected(1);
    end_time = pupil_t_corrected(end);

    % Define full window edges
    edges = start_time:window_duration:end_time;

    n_windows = length(edges) - 1;
    avg_gaze_angles = nan(n_windows, 2);
    midpoints = nan(n_windows, 1);

    for ii = 1:n_windows
        % Full window boundaries
        t_start = edges(ii);
        t_end = edges(ii+1);

        % Define middle 1-second segment within the full window
        segment_start = t_start + (window_duration - 1)/2; % 0.5 sec after start if window_duration=2
        segment_end = segment_start + 1;

        % Midpoint time of the middle segment
        midpoints(ii) = (segment_start + segment_end) / 2;

        % Find indices within the middle 1-second segment
        idx_in_segment = pupil_t_corrected >= segment_start & pupil_t_corrected < segment_end;

        if any(idx_in_segment)
            avg_gaze_angles(ii, :) = mean(gaze_angles_transformed(idx_in_segment, :), 1);
        else
            avg_gaze_angles(ii, :) = [NaN, NaN]; % no data for that segment
        end
    end
end

