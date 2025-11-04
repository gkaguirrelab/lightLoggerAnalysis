function frames = avi2frames(video_path, options)
    arguments
        video_path {mustBeText}
        options.grayscale {mustBeNumericOrLogical} = false; 
        options.verbose {mustBeNumericOrLogical} = false; 
        options.start {mustBeNumeric} = 1; 
        options.end {mustBeNumeric} = inf; 
    end 

    % Next, let's get the number of frames in the video 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    num_frames = double(Pi_util.inspect_video_frame_count(video_path)); 
    frame_size = cell(Pi_util.inspect_video_framesize(video_path)); 
    rows = double(frame_size{1}); 
    cols = double(frame_size{2}); 

    start = options.start; 
    if(options.end == inf)
        endpoint = num_frames; 
    else
        endpoint = options.end
    end

    % Allocate return variable 
    if(options.grayscale)
        frames = zeros((endpoint - start) + 1, rows, cols, 'uint8');
    else
        frames = zeros( (endpoint - start) + 1, rows, cols, 3, 'uint8');
    end

    % Iterate over the frames 
    for ii = start:endpoint 
        tic; 
        % Retrieve the frame from the video
        if(options.grayscale)
            frames(ii, :, :) = readAviFrame(video_path, ii, "grayscale", options.grayscale);
        else 
            frames(ii, :, :, :) = readAviFrame(video_path, ii, "grayscale", options.grayscale); 
        end        
        elapsed_seconds = toc; 

        if(options.verbose)
            fprintf("Frame %d/%d | Elapsed time: %.3f\n", ii, endpoint, elapsed_seconds); 
        end 

    end 

    return ; 

end 