function close(obj)
    arguments 
        obj 
    end 

    % Construct the output path of the video 
    temp_dir_path = fullfile(obj.Path, obj.filename); 
    video_output_path = fullfile(obj.full_video_path); 

    % Construct the output video 
    obj.utility_library.dir_to_video(temp_dir_path, video_output_path, obj.FrameRate); 

    % Once the video is made, remove the temp dir 
    rmdir(temp_dir_path, 's');

end 