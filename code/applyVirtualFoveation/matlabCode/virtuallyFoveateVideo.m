function virtuallyFoveateVideo(world_video, gaze_angles, gaze_offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, options)
% Virtually foveate desired frames of a video with given gaze angles
%
% Syntax:
%   virtuallyFoveateVideo(world_video, gaze_angles, offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, options)
%
% Description:
%   TODO 
%       
%   
% Inputs:
%   TODO                          
%
% Outputs:
%
%   NONE
%
% Examples:
%{
    world_video = "/Volumes/T7 Shield/scriptedIndoorOutdoorVideos/FLIC_2001/gazeCalibration/temporalFrequency/W.avi"; 
    path_to_recording_chunks = "/Volumes/EXTERNAL_1/FLIC_2001/gazeCalibration/temporalFrequency";
    output_path = "./testingVirtualFoveation.avi"; 
    gaze_angles = load("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2001/gazeCalibration/temporalFrequency/FLIC_2001_gazeCal_pupilData.mat").pupilData.radiusSmoothed.eyePoses.values; 
    offsets = load("/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2001/gazeCalibration/temporalFrequency/FLIC_2001_gazeCal_SceneGeometryMetadata.mat").gazeOffset;
    path_to_intrinsics = "/Users/zacharykelly/Documents/MATLAB/projects/lightLoggerAnalysis/data/intrinsics_calibration.mat"; 
    path_to_perspective_projection = "/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2001/gazeCalibration/temporalFrequency/FLIC_2001_gazeCal_perspectiveProjection.mat";
    virtuallyFoveateVideo(world_video, gaze_angles, offsets, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection)
    
%}

    arguments 
        world_video {mustBeText}; 
        gaze_angles {mustBeMatrix}; 
        gaze_offsets {mustBeNumeric}; 
        output_path {mustBeText}; 
        path_to_recording_chunks {mustBeText};
        path_to_intrinsics {mustBeText};
        path_to_perspective_projection {mustBeText}; 
        options.num_frames_to_process = [1, inf]; 
        options.pupil_fps {mustBeNumeric} = 120; 
        options.pupil_world_phase_offset {mustBeNumeric} = 0.005; 
    end     

    % Import the Python util library 
    virutal_foveation_util = import_pyfile(getpref("lightLoggerAnalysis", "virtual_foveation_util_path"));
    video_io_util = import_pyfile(getpref("lightLoggerAnalysis", "video_io_util_path"));

    % Create a video IO reader wrapper we will use to read into the original video
    world_frame_reader = videoIOWrapper(world_video); 

    % Now we will retrieve the start and end time of all of the sensors 
    start_ends = find_sensor_start_ends(virutal_foveation_util, path_to_recording_chunks); 
    world_start_end = start_ends.("world");
    pupil_start_end = start_ends.("pupil");

    % Create the T vectors that will be used to do mapping of gaze angles to frames, given 
    % that the sensors may sometimes be off on FPS 
    world_t = linspace(world_start_end(1), world_start_end(2), world_frame_reader.NumFrames);
    pupil_t = linspace(pupil_start_end(1), pupil_start_end(2), size(gaze_angles, 1));

    if(numel(world_t) ~= world_frame_reader.NumFrames)
        error("Miscalculation of world timestamps");
    end 

    % Next, add the slight offset that we measured in the calibration procedure. That is, the pupil 
    % is actually 0.005 seconds phase advanced
    pupil_t = pupil_t + options.pupil_world_phase_offset; 

    % Make a tempdir for the output so we can reconstruct to an .avi video later (MALTAB does not support this)
    temp_dir = random_string(5);
    if(~exist(temp_dir, "dir"))
        mkdir(temp_dir)
    end 

    % Initialize a blank frame we will use to pad frames that have nan gaze angles 
    blank_frame = zeros(world_frame_reader.Height, world_frame_reader.Width, 3, 'uint8'); 

    % Iterate over the world frames 
    start_frame = options.num_frames_to_process(1); 
    end_frame = world_frame_reader.NumFrames; 
    if(options.num_frames_to_process(2) ~= inf)
        end_frame = options.num_frames_to_process(2);
    end 

    tic; 
    for ii = start_frame:end_frame
        % Retrieve the world frame and its timestamp 
        world_frame = world_frame_reader.read(ii, 'grayscale', true); 
        world_timestamp = world_t(ii); 
        
        % Find the gaze angle that corresponds to this frame 
        [~, gaze_angle_idx] = min(abs(pupil_t - world_timestamp));
        gaze_angle = (gaze_angles(gaze_angle_idx, 1:2) .* [1, 1]) + ([gaze_offsets(1), gaze_offsets(2)] .* [-1, -1] );
        gaze_angle = [0, 0]; 

        % Virtually foveat the frame 
        virtually_foveated_frame = []; 
        if(any(isnan(gaze_angle)))
            virtually_foveated_frame = blank_frame; 
        else    
            virtually_foveated_frame = uint8(virtuallyFoveateFrame(world_frame, gaze_angle, path_to_intrinsics, path_to_perspective_projection)); 
        end 

        figure; 
        imshow(world_frame)
        hold on; 
        title(sprintf("Gaze angle: %f %f", gaze_angle(1), gaze_angle(2)));

        figure; 
        imshow(virtually_foveated_frame); 
        hold on; 

        % Write the resulting image out as a frame 
        imwrite(virtually_foveated_frame, fullfile(temp_dir, sprintf('frame_%d.png', ii)));
    end     

    % Convert the temp dir back into a playable video at the desired output location 
    video_io_util.dir_to_video(temp_dir, output_path, world_frame_reader.FrameRate); 

    % Remove the temp dir 
    rmdir(temp_dir, 's'); 

    elapsed_seconds = toc; 
    fprintf("Elapsed Time: %f seconds\n", elapsed_seconds); 

    return; 

end 

function start_ends = find_sensor_start_ends(virutal_foveation_util, path_to_recording_chunks)
    % Find the start ends 
    start_ends = struct(virutal_foveation_util.find_sensor_start_end_times(path_to_recording_chunks));
    field_names = fieldnames(start_ends);
    for ff = 1:numel(field_names)
        start_ends.(field_names{ff}) = double(start_ends.(field_names{ff})); 
    end 

end 

function s = random_string(n)
    chars = ['A':'Z' 'a':'z' '0':'9'];  % character set
    s = chars(randi(numel(chars), [1, n])); 
end



