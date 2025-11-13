function frame = readFrame(obj, options)
    arguments 
        obj 
        options.frameNum {mustBeNumeric} = []; 
        options.grayscale {mustBeNumericOrLogical} = false; 
        options.color {mustBeText} = "rgb"; 
    end 


    % If no frame num was provided, we simply read the last frame read + 1 
    was_frame_passed = ~isempty(options.frameNum); 
    if(was_frame_passed)
        frameNum = options.frameNum; 
    else
        frameNum = obj.last_frame_read + 1; 
    end 

    % Retrieve the frame from the video
    frame = squeeze(uint8(obj.utility_library.extract_frames_from_video( py.str(obj.full_video_path), {frameNum-1}, py.bool(options.grayscale))));
    py.gc.collect(); 
    
    % Convert to desired color if not grayscale 
    if (~options.grayscale)
        switch(options.color)
            % If user wants, just return bgr
            case "bgr"
                
            % Otherwise, need to flip channels 
            case "rgb"  
                frame = frame(:, :, [3 2 1]);
            
            otherwise 
                error("Unsupported color mode: %s", colorMode);

        end
    end 

    % Update the last read frame 
    obj.last_frame_read = frameNum; 


    return ; 

end
