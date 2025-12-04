close all; 
clear all;

subjectIDs = {"2001", "2002", "2003", "2004", "2005", "2006"};

for ii = 6:6 % numel(subjectIDs)
    subjectID = "FLIC_" + subjectIDs{ii}
    
    % Load in the world camera intrinscis 
    world_camera_intrinsics = load("/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat");

    % Load in the target positions in their intended angle form 
    target_pos_ang_intended = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_runData.mat", subjectID, subjectID)).taskData.gaze_target_positions_deg;

    % Load in the targets in their screen positions 
    % NOTE: If you do not have this, consult extract_gaze_stimulus.py
    target_pos_screen = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_gazeTargetsScreen", subjectID, subjectID)).gaze_targets;
    
    output_path = sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_perspectiveProjection.mat", subjectID, subjectID)
    
    % Calculate the perspective projection 
    perspective_projection = calculate_perspective_transform_w2e(world_camera_intrinsics, target_pos_ang_intended, target_pos_screen);

    % save the projection onto dropbox 
    save(output_path, "perspective_projection");
    
end 
