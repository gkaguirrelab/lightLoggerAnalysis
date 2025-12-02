close all; 
clear all;

subjectIDs = {"2001", "2003", "2004"};
april_tag_start_ends = {[9688, 14182], [9019, 11622], [10668, 13416]};  % For whole video, do [1, inf] 
experiment_portion_start_ends = {[20495, 49295], [17942, 46742], [18906, 47706]};
manual_offsets = {[0, 0], [0, 0], [0, 0]}

activity = "lunch";
video_type = "virtuallyFoveatedVideoAprilTag";

for ii = 1:numel(subjectIDs)
    subjectID = "FLIC_" + subjectIDs{ii}

    if(subjectID ~= "FLIC_2001")
        continue; 
    end 
    
    % First, we will define a path to the playable video of the world camera we want to virtually foveate
    % and its original chunks, to get the respective timestamps of all the sensors 
    world_video = sprintf("/Volumes/GKA spare/scriptedIndoorOutdoorVideos/%s/%s/temporalFrequency/W.avi", subjectID, activity)
    path_to_recording_chunks = sprintf("/Volumes/EXTERNAL_1/scriptedIndoorOutdoorVideos/%s/%s/temporalFrequency", subjectID, activity)
    
    % Load in the gaze angles and the constant offset
    gaze_angles = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/%s/temporalFrequency/%s_%s_pupilData_contrast1x25gamma1.mat", subjectID, activity, subjectID, activity)).pupilData.smoothPupilTime_02.eyePoses.values; 
    offsets = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_SceneGeometryMetadata.mat", subjectID, subjectID)).gazeOffset;

    % Define the output path where this video will write to 
    output_path = sprintf("/Users/zacharykelly/%s_%s_%s.avi", subjectID, activity, video_type) 

    % Load in the camera intrinscis of the world camera 
    path_to_intrinsics = "/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 

    % Load in the perspective projection object used to transform sensor positions to eye coordinates 
    % NOTE: If you do not have this, please consult calculate_perspective_transform_w2e.m 
    path_to_perspective_projection = sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_perspectiveProjection.mat", subjectID, subjectID)
    
    if(contains(lower(video_type), "apriltag"))
        start_end = april_tag_start_ends{ii};
    
    else
        start_end = experiment_portion_start_ends{ii};

    end 
       
    manual_offset = manual_offsets{ii}

    virtuallyFoveateVideo(world_video, gaze_angles, offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, "num_frames_to_process", start_end, "verbose", true, "manual_offset", manual_offset);

    
end 
