function generateVirtuallyFoveatedVideos(subjectIDs, options)
% Generate virtually foveated world-camera videos for one or more subjects
%
% Syntax:
%   generateVirtuallyFoveatedVideos(subjectIDs)
%   generateVirtuallyFoveatedVideos(subjectIDs, options)
%
% Description:
%   This function generates virtually foveated versions of world-camera
%   videos for one or more FLIC subjects. For each subject, the function
%   locates the corresponding raw and processed recording directories,
%   loads the processed world video and its timestamps, loads gaze samples
%   and blink events from the associated Neon outputs, and then calls
%   virtuallyFoveateVideo to render a new video in which the sampling of
%   the world video is centered on gaze position.
%
%   The function supports processing either user-specified frame ranges or
%   frame ranges loaded from a saved tag/task start-end structure. It also
%   supports optional per-subject manual gaze offsets, a projection-only
%   mode in which all gaze angles are set to zero, and optional output to
%   either the subject's processing folder or a user-specified directory.
%
% Inputs:
%   subjectIDs            - Cell array / array-like container. Subject IDs
%                           to process. Each entry should correspond to the
%                           numeric part of a subject folder name, e.g.
%                           {2001, 2003}, which will be converted to
%                           "FLIC_2001", "FLIC_2003", etc.
%
% Optional key/value pairs:
%   input_dir             - Char/string. Base directory containing the
%                           FLIC_raw and FLIC_processing folders. Defaults
%                           to "/Volumes/".
%   output_dir            - String. Output directory for the virtually
%                           foveated videos. If empty, videos are written
%                           to the subject's FLIC_processing folder.
%   activity              - String. Name of the activity folder to process.
%   video_type            - String. Type of video range to process. Must
%                           be one of:
%                               "tag"
%                               "task"
%                               "unspecified"
%                           If "unspecified", frame ranges are taken from
%                           options.start_ends. Otherwise, frame ranges are
%                           loaded from tag_task_start_end.mat.
%   start_ends            - Cell array. Per-subject start/end frame ranges
%                           used when video_type is "unspecified".
%   manual_offsets        - Cell array. Optional per-subject manual gaze
%                           offsets to apply in addition to the default
%                           [0, 0] offset.
%   verbose               - Logical. If true, print status information and
%                           relevant file paths while processing.
%   just_projection       - Logical. If true, set all gaze angles to zero
%                           so that the output reflects projection without
%                           gaze-dependent virtual foveation.
%   overwrite_existing    - Logical. If false, skip subjects whose output
%                           video already exists at the target path.
%   video_read_cache_size - Scalar. Cache size passed to
%                           virtuallyFoveateVideo for reading video data.
%
% Outputs:
%   none
%
% Examples:
%{
    generateVirtuallyFoveatedVideos( ...
        {2001, 2003}, ...
        "activity", "walkIndoorFoveate", ...
        "video_type", "task", ...
        "verbose", true, ...
        "overwrite_existing", false ...
    );

    generateVirtuallyFoveatedVideos( ...
        {2005}, ...
        "input_dir", "/Volumes/", ...
        "activity", "taskOutdoor", ...
        "video_type", "unspecified", ...
        "start_ends", {[1000 2500]}, ...
        "manual_offsets", {[1.5 -0.5]}, ...
        "output_dir", "/path/to/output" ...
    );
%}
    arguments
        subjectIDs
        options.input_dir = "/Volumes/"
        options.output_dir (1,1) string = ""
        options.activity (1,1) string = "unspecified"
        options.video_type (1,1) string ...
            {mustBeMember(options.video_type, ["tag","task","unspecified"])} = "unspecified"
        options.start_ends = {[1, inf]}
        options.manual_offsets = {}
        options.verbose (1,1) logical = true
        options.just_projection (1,1) logical = false
        options.overwrite_existing = false 
        options.video_read_cache_size (1,1) double = 1000
    end

    % Let's take out activity from options to not have to keep accessing the struct 
    activity = options.activity; 

    % Retrieve the start/ending frame 
    start_ends = options.start_ends;

    % Iterate over the subjectIDs enteterd 
    for ii = 1:numel(subjectIDs)
        % Format the subject IDs as they are on dropbox / NAS
        subjectID = "FLIC_" + string(subjectIDs{ii});

        % Now, let's make the path to this subject's files 
        subject_nas_path_raw = fullfile(options.input_dir, "FLIC_raw", "NEWscriptedIndoorOutdoorVideos2026", subjectID, activity); 
        subject_nas_path_processing = replace(subject_nas_path_raw, "FLIC_raw", "FLIC_processing");

        % Assert these folders exist 
        assert(isfolder(subject_nas_path_raw) && isfolder(subject_nas_path_processing), sprintf("Problem with %s or %s", subject_nas_path_raw, subject_nas_path_processing));

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
        egocentric_video_mapper_output_dir = fullfile(neon_dir_processing, "egocentric_mapper_results"); 

        % First, we will define a path to the playable video of the world camera we want to virtually foveate
        % and its timestamp vector output by the neon
        path_to_world_video = fullfile(gka_dir_processing, "W.avi"); 
        path_to_world_t = fullfile(egocentric_video_mapper_output_dir, "alternative_camera_timestamps.csv");
        
        if(options.verbose)
            fprintf("\twith world video paths:\n");
            fprintf("\t\tt: %s\n", path_to_world_t);
            fprintf("\t\tv: %s\n", path_to_world_video);
        end     
        assert(isfile(path_to_world_video) && isfile(path_to_world_t), sprintf("Problem with %s or %s", path_to_world_video, path_to_world_t));
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
        path_to_pupil_data = fullfile(egocentric_video_mapper_output_dir, "alternative_camera_gaze.csv");
        
        if(options.verbose)
            fprintf("\twith pupil data:\n");
            fprintf("\t\tpath: %s\n", path_to_pupil_data);
        end 
        assert(isfile(path_to_pupil_data))
        
        [pupil_t, gaze_angles] = load_gaze_angles(path_to_pupil_data, path_to_intrinsics); 

        % If we are just doing project and not full virtual foveation, set all gaze angles to the median
        if(options.just_projection)
            size_before = size(gaze_angles);
            median_gaze_angle = median(gaze_angles, 1);
            gaze_angles = repmat(median_gaze_angle, size(gaze_angles, 1), 1);
            size_after = size(gaze_angles);
            
            assert(isequal(size_before, size_after));

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
        projection_type = "virtuallyFoveated"; 
        if(options.just_projection)
            projection_type = "justProjection";
        end 

        output_filename = sprintf("%s_%s_%s_%s.avi", subjectID, activity, options.video_type, projection_type); 
        if(options.output_dir == "")
            output_path = fullfile(subject_nas_path_processing, output_filename); 
        else 
            output_path = fullfile(options.output_dir, output_filename); 

        end 

        % Skip existing videos if we do not want to overwrite them 
        if(~options.overwrite_existing && isfile(output_path))
            continue; 
        end 

        if(options.verbose)
            fprintf("\twith output path:\n");
            fprintf("\t\tpath: %s\n", output_path);
        end 

        % Determine which the range of frames to virtually foveate
        if(options.video_type == "unspecified") % If type unspecified, then we simply manual enter 
            start_end = start_ends{ii};    

        else 
            % Otherwise, gather the file with the start ends 
            % for this task 
            start_end_struct_path = fullfile(subject_nas_path_processing, "tag_task_start_end.mat");
            assert(isfile(start_end_struct_path));
            start_end = load(start_end_struct_path).tag_task_start_end.(options.video_type); 
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

        % Virtually foveate and output the video 
        sensor_t_matrix = {world_t, pupil_t};
        
        virtuallyFoveateVideo(path_to_world_video, sensor_t_matrix, gaze_angles, offsets, blnk_events, output_path, path_to_intrinsics,... 
                              "frames_to_process", start_end,...
                              "verbose", options.verbose,...
                              "manual_offset", manual_offset,...
                              "video_read_cache_size", options.video_read_cache_size...
                            );

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