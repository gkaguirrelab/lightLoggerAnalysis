function virtuallyFoveateVideo(world_video, gaze_angles, output_path, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, options)

    arguments 
        world_video {mustBeText}; 
        gaze_angles {mustBeMatrix}; 
        output_path {mustBeText}; 
        path_to_recording_chunks {mustBeText};
        path_to_intrinsics {mustBeText};
        path_to_perspective_projection {mustBeText}; 
        options.pupil_fps {mustBeNumeric} = 120; 

    end     

    % Import the Python util library 
    virutal_foveation_util = import_pyfile(getpref("lightLoggerAnalysis", "virtual_foveation_util_path"));

    % Create a video IO reader wrapper we will use to read into the original video
    world_frame_reader = videoIOWrapper(world_video); 

    % Now we will retrieve the start and end time of all of the sensors 
    start_ends = find_sensor_start_ends(virutal_foveation_util, path_to_recording_chunks); 
    world_start_end = start_ends.("world");
    pupil_start_end = start_ends.("pupil");

    % Find the time difference in the start of the world and the pupil. This will
    % inform us how to map between associated world frames and gaze angles 
    start_time_difference = pupil_start_end(1) - world_start_end(1); % TODO: Need to add offset from paper  
    start_frame_difference = start_time_difference / 1; 

    % Iterate over the world frames 
    for ii = 1:world_frame_reader.NumFrames
        % Retrieve the world frame and the pupil gaze angle 
        world_frame = world_frame_reader.read(ii); 
        gaze_angle = gaze_angles(ii + start_frame_difference, 1:2); 
        
        % Virtually foveat the frame 
        virtually_foveated_frame = virtuallyFoveateFrame(world_frame, gaze_angle, path_to_intrinsics, path_to_perspective_projection); 
        
        % Write this image out as bytes to a text file (since .avi video writing does not work in MATLAB )
    end 

    % Convert the temp text file back to the playable video 





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

function reconstruct_video_from_bytes(virutal_foveation_util, path_to_file, output_path)
    


end 



