function analyze_gaze_calibration()

    arguments
    
    end 

    % Load in the extracted pupil features 
    pupil_features = load('/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis').pupil_features; 

    % Extract just the theta and pi angles over time and convert to degrees 
    visual_angles = nan(size(pupil_features, 1), 2); 
    for ii = 1:numel{pupil_features}
        % Retrieve a given frame's features 
        frame_features = pupil_features{ii};
        
        % Save phi and theta
        visual_angles(ii, :) = rad2deg([frame_features.phi, frame_features.theta]); 
    end 



    % Load in the target positions in degrees (will make this dynamic in the future)
    degPositions = [ ...
        0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
        0, 20;   0, -20;   -20,   0;   20, 0; ...
        -15, 15;  15, 15;   -15, -15;   15, -15; ...

        -10,  10;   -10, -10;    10,  10;    10, -10; ...
        0,  10;    0, -10;   -10,   0;     10,   0; ...
        -5,   5;   5,   5;   -5,  -5;     5,  -5;     0,   0];



end 