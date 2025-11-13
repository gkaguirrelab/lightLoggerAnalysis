function frame = write(obj, frame)
    arguments 
        obj 
        frame
    end 

    % Save the frame dimensions for writing 
    if(isempty(obj.Height))
        [height, width] = size(frame); 
        obj.Height = height; 
        obj.Width = width; 
    end 

    % Because we are overloading this function to write out frames to a directory 
    % we need to keep track of a persistent (static) variable frame number 
    persistent frame_num 
    if(isempty(frame_num))
        frame_num = 1; 
    end     

    % Construct the path to this frame to be written 
    frame_path = fullfile(obj.Path, obj.filename, sprintf("%d.jpg", frame_num)); 

    % Write the frame 
    imwrite(frame, frame_path);

    % Increment the frame number 
    frame_num = frame_num + 1; 

    return ; 

end