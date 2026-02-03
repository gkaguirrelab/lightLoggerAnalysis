function generateVirtuallyFoveatedVideos(subjectIDs, options)

    arguments 
        subjectIDs 
        options.start_ends  = {[1, inf]};
        options.manual_offsets = {};
        options.activity (1, 1) string = "activityUnspecified"
        options.video_type (1,1) string {mustBeMember(options.video_type, ["virtuallyFoveatedVideoUnspecified", "virtuallyFoveatedVideo", "virtuallyFoveatedVideoAprilTag", "projectionNoFoveation"])} = "virtuallyFoveatedVideoUnspecified"
        options.output_dir (1, 1) string = "./"
        options.verbose (1, 1) logical = false
        options.testing (1, 1) logical = false; 
        options.just_projection (1, 1) logical = false; 
        options.non_contiguous_target_frames = {};
        options.video_read_cache_size = 500;
    end 

    % Let's get the dropbox base dir for any references we make to dropbox 
    dropbox_base_dir = getpref("lightLoggerAnalysis", "dropboxBaseDir"); 
    
    % Let's take out activity from options to not have to keep accessing the struct 
    activity = options.activity; 

    % Retrieve the start/ending frame 
    start_ends = options.start_ends;

    for ii = 1:numel(subjectIDs)
        subjectID = "FLIC_" + subjectIDs{ii};
        
        % First, we will define a path to the playable video of the world camera we want to virtually foveate
        % and its original chunks, to get the respective timestamps of all the sensors 
        path_to_world_video = "/Users/zacharykelly/Desktop/generated_videos/FLIC_2001_work_temporalSensitivity_1/W.avi" 
        world_t = load_world_timestamps("/Users/zacharykelly/Desktop/NeonWorkResult/alternative_camera_timestamps.csv");

        % Load in the camera intrinscis of the world camera 
        path_to_intrinsics = "~/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 

        % Next, we will load in the gaze angles (originally in px form, but we will convert)
        % and we will also load in the BLNK events 
        [pupil_t, gaze_angles] = load_gaze_angles("/Users/zacharykelly/Desktop/NeonWorkRESULT/alternative_camera_gaze.csv", path_to_intrinsics); %gaze_angles_struct.pupilData.(gaze_angles_field).eyePoses.values; 
        if(options.just_projection)
            gaze_angles(:, :) = 0; 
            if(any(gaze_angles(:)) ~= 0)
                error("Projection only mode was selected but non zero gaze angles detected");
            end 
        end 
        blnk_events = load_blnk_events("/Volumes/EXTERNAL_1/PilotWorkNeon/2026-01-29_16-16-52-59eb727c/blinks.csv");


        % Offset just set to 0, 0 since we assume Neon calculates this
        offsets = [0, 0]; 

        % Define the output path where this video will write to 
        output_path = fullfile(options.output_dir, sprintf("/%s_%s_%s.avi", subjectID, activity, options.video_type));

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

        % Retrieve the non contiguous target frames if provided 
        if(numel(options.non_contiguous_target_frames) == 0)
            non_contiguous_target_frames = [];
        else 
            non_contiguous_target_frames = options.non_contiguous_target_frames{ii};
        end 

        % Output information before we process if we would like 
        if(options.verbose)
            fprintf("Virtually foveating subject: %s\n", subjectID);
            fprintf("\tusing projection only (no gaze angle) mode: %d\n", options.just_projection); 
            fprintf("\tusing world video: %s\n", path_to_world_video);
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
        sensor_t_matrix = {world_t, pupil_t};
        virtuallyFoveateVideo(path_to_world_video, sensor_t_matrix, gaze_angles, offsets, blnk_events, output_path, path_to_intrinsics, path_to_perspective_projection,... 
                              "frames_to_process", start_end,...
                              "verbose", options.verbose,...
                              "manual_offset", manual_offset,...
                              "testing", options.testing,...
                              "non_contiguous_target_frames",non_contiguous_target_frames,...
                              "video_read_cache_size", options.video_read_cache_size);
    end 
end 


% Local function to load in the gaze angles from a given path 
% Gaze angles are stored in a .csv in px form, so we need to load 
% them in and also convert to deg 
function [pupil_t, gaze_angles] = load_gaze_angles(path, intrinsics_path)

    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');

    % Force timestamp column to int64 (preserves ns precision)
    opts = setvartype(opts, 'timestamp [ns]', 'int64');

    gaze_table = readtable(path, opts);

    pupil_t = gaze_table.('timestamp [ns]');   % seconds (safe in double))
    gx = double(gaze_table.('gaze x [px]'));
    gy = double(gaze_table.('gaze y [px]'));

    gaze_angles_px = [gx, gy];

    intr = load(intrinsics_path).camera_intrinsics_calibration.results.Intrinsics;
    gaze_angles = anglesFromIntrinsics(gaze_angles_px, intr);
end


% Load in the world timestamps 
function world_t = load_world_timestamps(path)
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');

    % Force timestamp column to int64 (preserves ns precision)
    opts = setvartype(opts, 'timestamp [ns]', 'int64');

    world_timestamps_table = readtable(path, opts); 
    world_t = int64(world_timestamps_table.('timestamp [ns]')); 
    
    return 

end 


% Local function to load in the timestamp ranges of blinks
function blnk_events = load_blnk_events(path)
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');

    % Force timestamp column to int64 (preserves ns precision)
    opts = setvartype(opts, {'start timestamp [ns]', 'end timestamp [ns]'}, 'int64');

    blnk_events_table = readtable(path, opts); 
    blnk_events = int64(blnk_events_table{:, {'start timestamp [ns]', 'end timestamp [ns]'}}); 

    assert(size(blnk_events, 2) == 2);

    return 
end 