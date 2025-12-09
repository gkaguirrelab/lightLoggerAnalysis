function generateVirtuallyFoveatedVideos(subjectIDs, start_ends, options)

    arguments 
        subjectIDs 
        start_ends 
        options.manual_offsets = {};
        options.activity (1, 1) string = "activityUnspecified"
        options.video_type (1,1) string {mustBeMember(options.video_type, ["virtuallyFoveatedVideoUnspecified", "virtuallyFoveatedVideo", "virtuallyFoveatedVideoAprilTag", "projectionNoFoveation"])} = "virtuallyFoveatedVideoUnspecified"
        options.output_dir (1, 1) string = "./"
        options.verbose (1, 1) logical = false
        options.testing (1, 1) logical = false; 
        options.just_projection (1, 1) logical = false; 
    end 

    % Let's get the dropbox base dir for any references we make to dropbox 
    dropbox_base_dir = getpref("lightLoggerAnalysis", "dropboxBaseDir"); 
    
    % Let's take out activity from options to not have to keep accessing the struct 
    activity = options.activity; 

    for ii = 1:numel(subjectIDs)
        subjectID = "FLIC_" + subjectIDs{ii};
        
        % First, we will define a path to the playable video of the world camera we want to virtually foveate
        % and its original chunks, to get the respective timestamps of all the sensors 
        path_to_world_video = sprintf("/Volumes/GKA spare/scriptedIndoorOutdoorVideos/%s/%s/temporalFrequency/W.avi", subjectID, activity);
        path_to_recording_chunks = sprintf("/Volumes/EXTERNAL_1/scriptedIndoorOutdoorVideos/%s/%s/temporalFrequency", subjectID, activity);
        
        % Load in the gaze angles and the constant offset
        gaze_angles_folder = fullfile(dropbox_base_dir, sprintf("FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/%s/temporalFrequency", subjectID, activity));
        gaze_angles_file = find_gaze_angles_file(gaze_angles_folder); 
        gaze_angles_path = fullfile(gaze_angles_folder, gaze_angles_file); 
        gaze_angles_struct = load(gaze_angles_path); % .pupilData.smoothPupilTime_02.eyePoses.values; 
        gaze_angles_field = gaze_angles_struct.pupilData.currentField; 
        gaze_angles = gaze_angles_struct.pupilData.(gaze_angles_field).eyePoses.values; 
        if(options.just_projection)
            gaze_angles(:, :, :, :) = 0; 
            if(any(gaze_angles(:)) ~= 0)
                error("Projection only mode was selected but non zero gaze angles detected");
            end 
        end 

        offsets_path = fullfile(dropbox_base_dir, sprintf("/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_SceneGeometryMetadata.mat", subjectID, subjectID)); 
        offsets = load(offsets_path).gazeOffset;

        % Define the output path where this video will write to 
        output_path = fullfile(options.output_dir, sprintf("/%s_%s_%s.avi", subjectID, activity, options.video_type));

        % Load in the camera intrinscis of the world camera 
        path_to_intrinsics = "~/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 

        % Load in the perspective projection object used to transform sensor positions to eye coordinates 
        % NOTE: If you do not have this, please consult calculate_perspective_transform_w2e.m 
        path_to_perspective_projection = fullfile(dropbox_base_dir, sprintf("/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_perspectiveProjection.mat", subjectID, subjectID));
        
        % Determine which the range of frames to virtually foveate
        start_end = start_ends{ii};
        
        % Retrieve the manual offsets for this video
        if(numel(options.manual_offsets) == 0)
            manual_offset = [0, 0]; 
        else 
            manual_offset = options.manual_offsets{ii};
        end 

        % Output information before we process if we would like 
        if(options.verbose)
            fprintf("Virtually foveating subject: %s\n", subjectID);
            fprintf("\tusing projection only (no gaze angle) mode: %d\n", options.just_projection); 
            fprintf("\tusing world video: %s\n", path_to_world_video);
            fprintf("\tusing recording chunks: %s\n", path_to_recording_chunks); 
            fprintf("\tusing gaze angles from field (%s): %s\n", gaze_angles_field, gaze_angles_path)
            fprintf("\tusing offsets: %s\n", offsets_path); 
            fprintf("\tusing intrinsics: %s\n", path_to_intrinsics); 
            fprintf("\tusing projection: %s\n", path_to_perspective_projection); 
            fprintf("\tusing start/end = [%d %d]\n", start_end(1), start_end(2));
            fprintf("\tusing manual offset = [%d %d]\n", manual_offset(1), manual_offset(2)); 
            fprintf("\twith testing output: %d\n", options.testing)
            fprintf("\toutputting to: %s\n", output_path);
        end     

        if(options.testing)
            warning("Testing has been enabled. This means that for each frame a figure will appear. This will make any sort of long video unable to finish due to amount of figures");
        end 

        % Virtually foveate and output the video 
        virtuallyFoveateVideo(path_to_world_video, gaze_angles, offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection,... 
                              "num_frames_to_process", start_end,...
                              "verbose", options.verbose,...
                              "manual_offset", manual_offset,...
                              "testing", options.testing...
                             );
    end 
end 

function gaze_angles_file = find_gaze_angles_file(folder)
    % Get all files in the folder
    files = dir(folder);

    % Keep only files (ignore folders)
    files = files(~[files.isdir]);

    % Extract filenames
    names = {files.name};

    % Select those containing the keyword
    matches = names(contains(names, "pupilData_contrast", 'IgnoreCase', true));

    % Output the list
    gaze_angles_file = matches;
    if((numel(gaze_angles_file)) == 0)
        error(sprintf("No gaze angle files found in %s", folder));
    end 

    if(numel(gaze_angles_file) > 1)
        error(sprintf("Multiple gaze angle files detected in %s", folder)); 
    end 

end 


% 2001 lunch april tag frames 
% 10548   11595 1.2340e+04  1.3100e+04      13863


% 2003 lunch april tag frames
% 9339 9814 10289 10810 11388

%{ 


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




%}