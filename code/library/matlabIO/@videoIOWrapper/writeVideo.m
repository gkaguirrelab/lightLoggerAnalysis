function writeVideo(obj, frame, options)
    arguments 
        obj 
        frame
        options.apply_floor_ceiling = false; 
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
    [height, width, channels] = size(frame);

    % Save the frame dimensions for writing 
    % if not already set 
    if(isempty(obj.Height))
        obj.Height = height; 
        obj.Width = width; 
    end 

    % Ensure the frame to write has matching shape as others 
    if(height ~= obj.Height || width ~= obj.Width)
        fprintf("Video Dimensions: (%d, %d)\n", obj.Height, obj.Width);
        fprintf("Frame Shape: (%d, %d)\n", height, width);

        error("Frames of video have become inhomogenously shaped")
    end 


    
    % Retrieve the number to name this frame 
    frame_num = obj.last_frame_written + 1; 

    % Construct the path to this frame to be written 
    frame_path = fullfile(obj.Path, obj.filename, sprintf("%d.png", frame_num)); 

    % If we apply floor ceiling, 
    % in a 2D image all pixels are clamped to 0, 255 
    % in a 3D image, if any channel of a given pixel is <= 0 
    % all channels for that pixel become 0 
    % likewise, if any channel value is >= 255, all channel values 
    % for that pixel become 255
    if(options.apply_floor_ceiling)
        if(channels == 1)
            frame(frame <= 0) = 0;
            frame(frame >= 255) = 255;

        elseif(channels == 3)
            floor_mask = any(frame <= 0, 3);
            ceiling_mask = any(frame >= 255, 3);

            for cc = 1:channels
                current_channel = frame(:, :, cc);
                current_channel(floor_mask) = 0;
                current_channel(ceiling_mask) = 255;
                frame(:, :, cc) = current_channel;
            end

        else
            error("Unsupported number of channels: %d", channels);
        end
    end

    % Write the frame 
    imwrite(frame, frame_path);

    % Update the last written frame 
    obj.last_frame_written = frame_num; 

    return ; 

end
