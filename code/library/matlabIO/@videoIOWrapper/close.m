function close(obj)
    arguments 
        obj 
    end 

    % If obj.Path was not defined, this obj is being destroyed before properly created, 
    % so let's just close 
    if(isempty(obj.Path))
        return;
    end 

    % Construct the temporary input/output paths of the video 
    read_temp_file_path = obj.temporary_reading_hdf5_filepath; 
    write_temp_dir_path = fullfile(obj.Path, obj.filename); 
    video_output_path = fullfile(obj.full_video_path); 

    % If no video input/output was generated, just skip 
    if(isempty(read_temp_file_path) && isempty(write_temp_dir_path))
        return; 
    end 

    % Delete the temporary reading file 
    if(~isempty(read_temp_file_path) && exist(read_temp_file_path, "file"))
        delete(read_temp_file_path);
    end     

    % Check if the temp dir existed for if we were using, 
    % this in the writing mode, if it did, we remove it, if it did not, we do nothing 
    if(~isempty(write_temp_dir_path) && exist(write_temp_dir_path, "dir"))
        % List everything in the directory
        files = dir(write_temp_dir_path);

        % Remove "." and ".."
        files = files(~ismember({files.name}, {'.','..'}));

        % If there is a non-zero amount of frames, make the video
        if(~isempty(files))
            obj.utility_library.dir_to_video(write_temp_dir_path, video_output_path, double(obj.FrameRate)); 
        end 

        % Once the video is made, remove the temp dir 
        rmdir(write_temp_dir_path, 's');
    end

end 