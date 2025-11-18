function open(obj)
    arguments 
        obj 
    end     

    % If we are reading and not writing, just do nothing 
    if(~strcmp(obj.mode, 'write'))
        return; 
    end     

    % Construct the path to the temporary file where we will write 
    % out frames to disk 
    temp_dir_path = fullfile(obj.Path, obj.filename); 
    
    if ~exist(temp_dir_path, 'dir')
        mkdir(temp_dir_path);
    end
end