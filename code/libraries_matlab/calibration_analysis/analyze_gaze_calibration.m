function analyze_gaze_calibration()

    arguments
    
    end 

    % Load in the extracted pupil features 
    pupil_features = load('/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis/gaze_cal_pupil_features').pupil_features; 

    % Load in the target positions in degrees (will make this dynamic in the future)
    deg_positions = [ ...
        0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
        0, 20;   0, -20;   -20,   0;   20, 0; ...
        -15, 15;  15, 15;   -15, -15;   15, -15; ...

        -10,  10;   -10, -10;    10,  10;    10, -10; ...
        0,  10;    0, -10;   -10,   0;     10,   0; ...
        -5,   5;   5,   5;   -5,  -5;     5,  -5;     0,   0];


    % Retrieve the features we want to analyze 
    [gaze_angles, confidence_measures] = flatten_features(pupil_features); 

    % Display the features in their raw state 
    display_raw_results(gaze_angles, deg_positions, confidence_measures);

    % Display the data after making desired transformations 
    % (splicing out extra parts of video, mean center, remove low confidence points)
    display_adjusted_results(gaze_angles, deg_positions, confidence_measures); 


    return ; 

end 


% Local function to extract the desired features from the list 
function [gaze_angles, confidence_measures] = flatten_features(pupil_features)
    % Initialize return arrays for the features we will 
    % flatten
    gaze_angles = nan(size(pupil_features, 1), 2); 
    confidence_measures = nan(size(pupil_features, 1), 2); 

    % Iterate over the pupil features structs 
    % and extract the ones we care about 
    for ii = 1:numel(pupil_features)
        % Retrieve a given frame's features 
        frame_features = pupil_features{ii};
        
        % Save phi and theta
        gaze_angles(ii, :) = rad2deg([frame_features.phi, frame_features.theta]); 

        % Save the confidence metrics
        confidence_measures(ii, :) = [frame_features.confidence, frame_features.model_confidence];
    end 

    return ; 
end 

% Local function to display raw results before adjustments
function display_raw_results(gaze_angles, deg_positions, confidence_measures)
    % First, make a plot of the raw features. No clipping, no confidence threshold, no mean subtraction 
    figure; 
    tiled_fig_handle = tiledlayout(1, 3); 
    tiled_fig_handle.Title.String = "Transformed Data"; 
    tiled_fig_handle.Title.FontWeight = "bold"; 

    % Move to the first tile
    % and plot the gaze angles over time and the degree positions
    % of the targets 
    nexttile; 
    title("Gaze Angle by Frame"); 
    hold on; 
    plot(gaze_angles(:, 1), 'x', 'DisplayName', 'Phi'); 
    plot(gaze_angles(:, 2), 'x', 'DisplayName', 'Theta');
    xlabel("Time [frame]")
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
    plot(confidence_measures(:, 1), 'x', 'DisplayName', 'Confidence');
    plot(confidence_measures(:, 2), 'x', 'DisplayName', 'Model Confidence');
    xlabel("Time [frame]"); 
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
function display_adjusted_results(gaze_angles, deg_positions, confidence_measures)
    % Now, plot the adjusted features after transformation has been applied
    figure; 
    tiled_fig_handle = tiledlayout(1, 3); 
    tiled_fig_handle.Title.String = "Raw Data (No splicing, thresholding, or value adjustment)"; 
    tiled_fig_handle.Title.FontWeight = "bold"; 

    % 1. Trim out just the gaze calibration portion of the video 
    start_idx = 10800; 
    end_idx = size(gaze_angles, 1) - 3600; 
    gaze_angles = gaze_angles(start_idx:end_idx, :); 
    confidence_measures = confidence_measures(start_idx:end_idx, :); 

    % 2. Subtract the mean theta from the theta vector, and the mean phi from the phi vector
    mean_phi_theta = mean(gaze_angles, 1); 
    gaze_angles_transformed = gaze_angles;
    gaze_angles_transformed(:, 1) = gaze_angles_transformed(:, 1) - mean_phi_theta(1);
    gaze_angles_transformed(:, 2) = gaze_angles_transformed(:, 2) - mean_phi_theta(2); 

    % 3. Remove low confidence frames
    good_idx = (confidence_measures(:, 1) > 0.90) & (confidence_measures(:, 2) > 0.9);
    gaze_angles_transformed = gaze_angles_transformed(good_idx, :);
    confidence_measures = confidence_measures(good_idx, :);


    % Move to the first tile
    % and plot the gaze angles over time and the degree positions
    % of the targets 
    nexttile; 
    title("Gaze Angle by Frame"); 
    hold on; 
    plot(gaze_angles_transformed(:, 1), 'x', 'DisplayName', 'Phi'); 
    plot(gaze_angles_transformed(:, 2), 'x', 'DisplayName', 'Theta');
    xlabel("Time [frame]")
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
    plot(confidence_measures(:, 1), 'x', 'DisplayName', 'Confidence');
    plot(confidence_measures(:, 2), 'x', 'DisplayName', 'Model Confidence');
    xlabel("Time [frame]"); 
    ylabel("Confidence"); 
    ylim([0, 1]);
    legend show; 
    hold off; 

    return ;
end 