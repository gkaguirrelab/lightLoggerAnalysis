function frame = readAviFrame(video_path, idx, options)
    arguments 
        video_path {mustBeText}; 
        idx {mustBeNumeric}; 
        options.grayscale {mustBeNumericOrLogical} = false; 
    end 

    % Only import Pi_util once
    persistent Pi_util
    if(isempty(Pi_util))
        Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    end 

    % Retrieve the frame from the video
    frame = squeeze(uint8(Pi_util.extract_frames_from_video(video_path, {idx}, options.grayscale)));

    return ; 
end 