function frame = readAviFrame(video_path, idx, options)
    arguments 
        video_path {mustBeText}; 
        idx {mustBeNumeric}; 
        options.grayscale {mustBeNumericOrLogical} = false; 
        options.color {mustBeText} = "rgb"
    end 

    % Only import Pi_util once
    persistent Pi_util
    if(isempty(Pi_util))
        Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    end 

    % Retrieve the frame from the video
    frame = squeeze(uint8(Pi_util.extract_frames_from_video(video_path, {idx-1}, options.grayscale)));
    
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
        end
    end 

    return ; 
end 