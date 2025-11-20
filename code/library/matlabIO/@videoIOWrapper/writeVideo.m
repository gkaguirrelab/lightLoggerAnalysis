function writeVideo(obj, frame)
    arguments 
        obj 
        frame
    end 

    % Ensure the FPS is set by this point 
    if(isempty(obj.FrameRate))
        error("The FrameRate proprerity of the writer was not set before writing.")
    end 

    % Ensure the frame object is not empty 
    if(isempty(frame))
        error("Frame argument is empty"); 
    end 

    % If frame is a struct (such as in transparentTrack, retrieve the field with the data)
    if(isstruct(frame))
        frame = frame.cdata; 
    end 
    
    % Retrieve the size of the frame to write 
    [height, width] = size(frame);

    % Save the frame dimensions for writing 
    % if not already set 
    if(isempty(obj.Height))
        obj.Height = height; 
        obj.Width = width; 
    end 

    % Ensure the frame to write has matching shape as others 
    if(height ~= obj.Height || width ~= obj.Width)
        error("Frames of video have become inhomogenously shaped")
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