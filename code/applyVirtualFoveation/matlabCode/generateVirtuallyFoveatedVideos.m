function generateVirtuallyFoveatedVideos(subjectIDs, options)

    arguments
    subjectIDs
    options.output_dir (1,1) string = ""
    options.activity (1,1) string = "unspecified"
    options.video_type (1,1) string ...
        {mustBeMember(options.video_type, ["april","task","unspecified"])} = "unspecified"
    options.start_ends = {[1, inf]}
    options.manual_offsets = {}
    options.verbose (1,1) logical = true
    options.testing (1,1) logical = false
    options.just_projection (1,1) logical = false
    options.non_contiguous_target_frames = {}
    options.video_read_cache_size (1,1) double = 1000
end

    % Let's get the dropbox base dir for any references we make to dropbox 
    dropbox_base_dir = getpref("lightLoggerAnalysis", "dropboxBaseDir"); 
    NAS_base_dir = "/Volumes/";

    % Let's take out activity from options to not have to keep accessing the struct 
    activity = options.activity; 

    % Retrieve the start/ending frame 
    start_ends = options.start_ends;

    % Iterate over the subjectIDs enteterd 
    for ii = 1:numel(subjectIDs)
        % Format the subject IDs as they are on dropbox 
        subjectID = "FLIC_" + string(subjectIDs{ii});

        % Now, let's make the path to this subject's files 
        subject_nas_path_raw = fullfile(NAS_base_dir, "FLIC_raw", "scriptedIndoorVideos", subjectID, activity); 
        subject_nas_path_processing = replace(subject_nas_path_raw, "FLIC_raw", "FLIC_processing");

        % Assert these folders exist 
        assert(isfolder(subject_nas_path_raw) && isfolder(subject_nas_path_processing));

        % Define shortcuts to neon/gka for raw and processing 
        gka_dir_raw = fullfile(subject_nas_path_raw, "GKA"); 
        neon_dir_raw = fullfile(subject_nas_path_raw, "Neon"); 

        gka_dir_processing = fullfile(subject_nas_path_processing, "GKA"); 
        neon_dir_processing = fullfile(subject_nas_path_processing, "Neon"); 

        % Get the name of the neon recording (this is the data when it was recorded)
        neon_raw_recording_name = getSingleSubfolder(neon_dir_raw); 

        if(options.verbose)
            fprintf("Processing subject: %s\n", subjectID); 
            fprintf("\tWith directories:\n"); 
            fprintf("\t\tFLIC_raw: \n");
            fprintf("\t\t\tGKA: %s\n", gka_dir_raw); 
            fprintf("\t\t\tNeon: %s\n", neon_dir_raw); 
            fprintf("\t\tFLIC_processing:\n");
            fprintf("\t\t\tGKA: %s\n", gka_dir_processing);
            fprintf("\t\t\tNeon: %s\n", neon_dir_processing) 
        end 

        assert(isfolder(gka_dir_raw) && isfolder(neon_dir_raw) && isfolder(gka_dir_processing) && isfolder(neon_dir_processing)); 

        % Also save the directory where the output of the egocentric video mapper lives 
        egocentric_video_mapper_output_dir = fullfile(neon_dir_processing, "egocentricVideoMapperResults"); 

        % First, we will define a path to the playable video of the world camera we want to virtually foveate
        % and its timestamp vector output by the neon
        path_to_world_video = fullfile(gka_dir_processing, "W.avi"); 
        path_to_world_t = fullfile(egocentric_video_mapper_output_dir, "alternative_camera_timestamps.csv");
        
        if(options.verbose)
            fprintf("\twith world video paths:\n");
            fprintf("\t\tt: %s\n", path_to_world_t);
            fprintf("\t\tv: %s\n", path_to_world_video);
        end     
        assert(isfile(path_to_world_video) && isfile(path_to_world_t));
        world_t = load_world_timestamps(path_to_world_t);

        % Load in the camera intrinscis of the world camera 
        path_to_intrinsics = "~/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 

        if(options.verbose)
            fprintf("\twith intrinsics:\n");
            fprintf("\t\tpath: %s\n", path_to_intrinsics);
        end 
        assert(isfile(path_to_intrinsics));

        % Next, we will load in the gaze angles (originally in px form, but we will convert)
        % and we will also load in the BLNK events 
        path_to_pupil_data = fullfile(egocentric_video_mapper_output_dir, "alternative_camera_gaze.csv")
        
        if(options.verbose)
            fprintf("\twith pupil data:\n");
            fprintf("\t\tpath: %s\n", path_to_pupil_data);
        end 
        assert(isfile(path_to_pupil_data))
        
        [pupil_t, gaze_angles] = load_gaze_angles(path_to_pupil_data, path_to_intrinsics); 

        % If we are just doing project and not full virtual foveation, set all gaze angles to 0 
        if(options.just_projection)
            gaze_angles(:, :) = 0; 
            if(any(gaze_angles(:)) ~= 0)
                error("Projection only mode was selected but non zero gaze angles detected");
            end 
        end 

        % Next we will load in the blink events 
        path_to_blnk_data = fullfile(neon_dir_raw, neon_raw_recording_name, "blinks.csv"); 
        if(options.verbose)
            fprintf("\twith blink events:\n");
            fprintf("\t\tpath: %s\n", path_to_blnk_data);
        end 
        
        assert(isfile(path_to_blnk_data));
        blnk_events = load_blnk_events(path_to_blnk_data);


        % Constant Offset just set to 0, 0 since we assume Neon calculates this
        offsets = [0, 0]; 

        % Define the output path where this video will write to 
        output_filename = sprintf("/%s_%s_%s_virtuallyFoveated.avi", subjectID, activity, options.video_type); 
        if(options.output_dir == "")
            output_path = fullfile(subject_nas_path_processing, output_filename); 
        else 
            output_path = fullfile(options.output_dir, output_filename); 

        end 

        if(options.verbose)
            fprintf("\twith output path:\n");
            fprintf("\t\tpath: %s\n", output_path);
        end 

        % Load in the perspective projection object used to transform sensor positions to eye coordinates 
        % NOTE: If you do not have this, please consult calculate_perspective_transform_w2e.m 
        path_to_perspective_projection = fullfile(dropbox_base_dir, sprintf("/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/%s/gazeCalibration/temporalFrequency/%s_gazeCal_perspectiveProjection.mat", subjectID, subjectID));


        % Determine which the range of frames to virtually foveate
        if(options.video_type == "unspecified") % If type unspecified, then we simply manual enter 
            start_end = start_ends{ii};    

        else 
            % Otherwise, gather the file with the start ends 
            % for this task 
            start_end_struct_path = fullfile(subject_nas_path_processing, "tag_task_start_end.mat");
            assert(isfile(start_end_struct_path));
            start_end = load(start_end_struct).(options.video_type); 
        end 

        if(options.verbose)
            fprintf("\twith start_ends: \n"); 
            fprintf("\t\tstart: %d\n", start_end(1));
            fprintf("\t\tend: %d\n", start_end(end));
        end 

        


        
        % Retrieve the manual offsets for this video (further manual adjsutments per video beyond 
        % individual offset per participant)
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


        if(options.testing)
            warning("Testing has been enabled. This means that for each frame a figure will appear. This will make any sort of long video unable to finish due to amount of figures");
        end 

        % Virtually foveate and output the video 
        sensor_t_matrix = {world_t, pupil_t};
        
        %{
        virtuallyFoveateVideo(path_to_world_video, sensor_t_matrix, gaze_angles, offsets, blnk_events, output_path, path_to_intrinsics, path_to_perspective_projection,... 
                              "frames_to_process", start_end,...
                              "verbose", options.verbose,...
                              "manual_offset", manual_offset,...
                              "testing", options.testing,...
                              "non_contiguous_target_frames",non_contiguous_target_frames,...
                              "video_read_cache_size", options.video_read_cache_size);

        %}
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

% Local function to get the name of the Neon recording given 
% its output folder 
function folderName = getSingleSubfolder(parentDir)
%GETSINGLESUBFOLDER Returns the single subfolder name inside a directory.
%
%   folderName = getSingleSubfolder(parentDir)
%
%   - parentDir must exist
%   - Exactly one subfolder (excluding . and ..) must exist
%   - Returns the name (not full path)
    % Ensure directory exists
    assert(isfolder(parentDir), ...
        "Input path does not exist or is not a folder.")

    % Get directory listing
    listing = dir(parentDir);

    % Keep only folders and remove "." and ".."
    isSubfolder = [listing.isdir] & ...
                  ~ismember({listing.name}, {'.','..'});

    subfolders = listing(isSubfolder);

    % Assert exactly one subfolder exists
    assert(numel(subfolders) == 1, ...
        "Directory must contain exactly one subfolder. Found %d.", ...
        numel(subfolders));

    % Return folder name
    folderName = string(subfolders.name);
end