function writeVideo(obj, frame)
    arguments 
        obj 
        frame
    end 

    % If frame is a struct (such as in transparentTrack, retrieve the field with the data)
    if(isstruct(frame))
        frame = frame.cdata; 
    end 

    % Save the frame dimensions for writing 
    if(isempty(obj.Height))
        [height, width] = size(frame); 
        obj.Height = height; 
        obj.Width = width; 
    end 
    
    % Retrieve the number to name this frame 
    frame_num = obj.last_frame_written + 1; 

    % Construct the path to this frame to be written 
    frame_path = fullfile(obj.Path, obj.filename, sprintf("%d.jpg", frame_num)); 

    % Write the frame 
    imwrite(frame, frame_path);

    % Update the last written frame 
    obj.last_frame_written = frame_num; 

    return ; 

end