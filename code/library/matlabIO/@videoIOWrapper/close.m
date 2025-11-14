function close(obj)
    arguments 
        obj 
    end 

    % Construct the output path of the video 
    temp_dir_path = fullfile(obj.Path, obj.filename); 
    video_output_path = fullfile(obj.full_video_path); 

    % Check if the temp dir existed, if it did, we remove it, if it did not, we do nothing 
    if(exist(temp_dir_path, "dir"))
        % List everything in the directory
        files = dir(temp_dir_path);

        % Remove "." and ".."
        files = files(~ismember({files.name}, {'.','..'}));

        % If there is a non-zero amount of frames, make the video
        if(~isempty(files))
            obj.utility_library.dir_to_video(temp_dir_path, video_output_path, obj.FrameRate); 
        end 


        % Once the video is made, remove the temp dir 
        rmdir(temp_dir_path, 's');
    end 

end 