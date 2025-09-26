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
    good_idx = (confidence_measures(:, 1) > 0.90) & (confidence_measures(:, 2) > 0.90);
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
    % and plot the gaze angles over time and the degree positions
    % of the targets
    nexttile;
    % Repeat each position 3 times
    repeats = 3;
    positionsRepeated = repmat(deg_positions, repeats, 1);

    % Create time vector: each point lasts 1 second
    dotDuration = 2.2; %seconds approximately. measured from world camera recording 
    nDotsPerRep = size(deg_positions,1);
    interRepDelay = 2; % seconds between repetitions, blank screen.
    appearanceTimes = getDotAppearanceTimes(dotDuration, nDotsPerRep, repeats, interRepDelay);

    % Plot
    title("Target Position");
    hold on;
    plot(appearanceTimes, positionsRepeated(:,1), 'ro-', 'LineWidth', 1.5); hold on;
    plot(appearanceTimes, positionsRepeated(:,2), 'bo-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    legend('Phi', 'Theta');
    hold off;

     % Move to the third tile
    % and average gaze angle for middle of when dots should have been
    % presented
    nexttile
    title('Avg. Gaze Angle in window');
    hold on;
    [midpoints, avg_gaze_angles] = average_gaze_middle_segment(gaze_angles_transformed, pupil_t_corrected, appearanceTimes, dotDuration);
    
    plot(midpoints, avg_gaze_angles(:, 1), 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Phi');
    plot(midpoints(:, 1), avg_gaze_angles(:, 2), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Theta');
    xlabel("Time [s]")
    ylabel("Angle [deg]")
    legend show;
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

function appearanceTimes = getDotAppearanceTimes(dotDuration, nDotsPerRep, repeats, interRepDelay)
% getDotAppearanceTimes returns a vector of timestamps when each dot appears on screen
%
% Inputs:
%   dotDuration   - duration (in seconds) each dot is visible (~2.15s)
%   nDots         - number of dots per repetition
%   repetitions   - number of repetitions
%   interRepDelay - delay (in seconds) between repetitions (e.g., 2 seconds)
%
% Output:
%   appearanceTimes - (nDots * repetitions) x 1 vector of appearance times in seconds

    totalDots = nDotsPerRep * repeats;
    appearanceTimes = zeros(totalDots, 1);

    time = 0;  % initialize time counter

    for rep = 1:repeats
        for ii = 1:nDotsPerRep
            appearanceTimes((rep-1)*nDotsPerRep + ii) = time;  % dot appears at current time
            time = time + dotDuration;                   % advance time by dot duration
        end
        time = time + interRepDelay;                      % add delay between repetitions
    end
end

function [midpoints, avg_gaze_angles] = average_gaze_middle_segment(gaze_angles_transformed, pupil_t_corrected, appearanceTimes, dotDuration)
    % average_gaze_middle_segment_updated averages gaze angles over the middle 1-second segment of each dot's appearance
    %
    % Inputs:
    %   gaze_angles_transformed: Nx2 matrix of gaze angles [phi, theta]
    %   pupil_t_corrected: Nx1 vector of timestamps (seconds)
    %   appearanceTimes: Mx1 vector of timestamps when each dot appears
    %   dotDuration: duration (in seconds) each dot is visible
    %
    % Outputs:
    %   midpoints: Mx1 vector of midpoint times for each averaged segment
    %   avg_gaze_angles: Mx2 matrix of average gaze angles [phi, theta] per middle segment
    
    n_targets = length(appearanceTimes);
    avg_gaze_angles = nan(n_targets, 2);
    midpoints = nan(n_targets, 1);
    
    % Define the middle 1-second segment
    middle_segment_duration = 1;
    % Calculate the offset to find the middle of the dot's appearance
    % For a dotDuration of 2.15s, the middle 1s segment starts at (2.15-1)/2 = 0.575s
    start_offset = (dotDuration - middle_segment_duration) / 2;
    
    for ii = 1:n_targets
        % Calculate the start and end of the middle 1-second segment for the current dot
        segment_start = appearanceTimes(ii) + start_offset;
        segment_end = segment_start + middle_segment_duration;
        
        % Calculate the midpoint of the segment
        midpoints(ii) = (segment_start + segment_end) / 2;
        
        % Find indices of gaze data within this segment
        idx_in_segment = pupil_t_corrected >= segment_start & pupil_t_corrected < segment_end;
        
        if any(idx_in_segment)
            % Average the gaze angles for the data in the segment
            avg_gaze_angles(ii, :) = mean(gaze_angles_transformed(idx_in_segment, :), 1);
        else
            % If no data is found, set the average to NaN
            avg_gaze_angles(ii, :) = [NaN, NaN];
        end
    end
end

