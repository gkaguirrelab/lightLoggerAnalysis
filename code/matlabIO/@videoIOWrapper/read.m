function frame = read(obj, frameNum, options)
    arguments 
        obj 
        frameNum {mustBeNumeric}; 
        options.grayscale {mustBeNumericOrLogical} = false; 
        options.color {} = "rgb"; 
    end 

    % Retrieve the frame from the video
    frame = squeeze(uint8(obj.utility_library.extract_frames_from_video( py.str(obj.full_video_path), {frameNum-1}, py.bool(options.grayscale))));
    
    % Convert to desired color if not grayscale 
    if (~options.grayscale)
        switch(options.color)
            % If user wants, just return bgr
            case "bgr"
                return ; 

            % Otherwise, need to flip channels 
            case "rgb"
                frame = frame(:, :, [3 2 1]);
                return ; 
            otherwise 
                error("Unsupported color mode: %s", colorMode);

        end
    end 

    return ; 

end