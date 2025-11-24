close all; 
clear all;

subjectIDs = {"2005"};
start_ends = {[10640, 10640]};
manual_offsets = { [-7.5, -7.5] }

video_type = "testing";

for ii = 1:numel(subjectIDs)
    subjectID = "FLIC_" + subjectIDs{ii}
    
    % First, we will define a path to the playable video of the world camera we want to virtually foveate
    % and its original chunks, to get the respective timestamps of all the sensors 
    world_video = sprintf("/Volumes/GKA Spare/scriptedIndoorOutdoorVideos/%s/walkIndoor/temporalFrequency/W.avi", subjectID)
    path_to_recording_chunks = sprintf("/Volumes/EXTERNAL_1/%s/walkIndoor/temporalFrequency", subjectID)

    % Load in the gaze angles and the constant offset we will apply to the gaze angles
    if(subjectID ~= "FLIC_2002")
        gaze_angles = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/walkIndoor/temporalFrequency/%s_walkIndoor_pupilData_contrast-1x5.mat", subjectID, subjectID)).pupilData.radiusSmoothed.eyePoses.values;      
    else
        gaze_angles = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/walkIndoor/temporalFrequency/%s_walkIndoor_pupilData_contrast-1x5.mat", subjectID, subjectID)).pupilData.sceneConstrained.eyePoses.values ;
    end 

    offsets = load(sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_SceneGeometryMetadata.mat", subjectID, subjectID)).gazeOffset;

    % Define the output path where this video will write to 
    output_path = sprintf("/Users/zacharykelly/%s_walkIndoor_%s.avi", subjectID, video_type) % change this to testing when you do single frames

    % Load in the camera intrinscis of the world camera 
    path_to_intrinsics = "/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 

    % Load in the perspective projection object used to transform sensor positions to eye coordinates 
    % NOTE: If you do not have this, please consult calculate_perspective_transform_w2e.m 
    path_to_perspective_projection = sprintf("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_perspectiveProjection.mat", subjectID, subjectID)

    start_end = start_ends{ii}
    manual_offset = manual_offsets{ii};

    virtuallyFoveateVideo(world_video, gaze_angles, offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, "num_frames_to_process", start_end, "verbose", true, "manual_offset", manual_offset, "testing", true);

    
end 
