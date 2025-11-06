function virtuallyFoveateVideo(world_video, gaze_angles, path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection, options)

    arguments 
        world_video {mustBeText}; 
        gaze_angles {mustBeMatrix}; 
        output_path {mustBeText}; 
        path_to_recording_chunks {mustBeText};
        path_to_intrinsics {mustBeText};
        path_to_perspective_projection {mustBeText}; 
        options.pupil_fps {mustBeDouble} = 120; 

    end 

    % Create a video IO reader wrapper we will use to read into the original video
    world_frame_reader = videoIOWrapper(world_video); 

    % Now we will retrieve the start and end time of all of the sensors 
    

    %% Iterate over the frames of the video, reading them in 
    %parfor ii = 1:world_frame_reader.NumFrames


    %   end 


    return; 

end 





