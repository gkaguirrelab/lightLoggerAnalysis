function global_means = calculate_global_channels_mean(video_path, options)
    arguments 
        video_path string;         
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY", "LMS", "L+M+S", "L-M", "a", "c_lm", "c_s"])}; 
        options.sum_of {mustBeMember(options.sum_of, ["value", "loge"])}; 
        options.channels = [1, 2, 3]; 
        options.verbose = false 
        options.start_end = false; 
    end 

    % First, let's open a video reader to stream frames from the video 
    video_reader = videoIOWrapper(video_path, 'ioAction', 'read'); 

    % Initialize the sum variable. It should be a 
    % vector consistening of the number of channels 
    % we would like to mean over 
    global_sums = zeros(numel(options.channels), 1); 

    % Define the transformation that will handle any modifications 
    % to the value that we are summing (suc as by taking the natural log)
    switch(options.sum_of)
        
        % If we just want the value, 
        % transformation is the identity transformation
        case "value"
            transformation = @(x) x; 

        % If we want the natural log 
        % it is as follows 
        case "loge"
            transformation = @(x) log(x + 10e-9);

        otherwise
            error("Unsupported transformation %s", options.sum_of);

    end 



    % Next, let's iterate over the frames of the video 
    if options.verbose
        h_wait = waitbar(0, 'Calculating Channel Means | Processing video frames...');
    end

    n_frames = video_reader.NumFrames; 
    if(~islogical(options.start_end))
        start = options.start_end(1); 
        end_ = options.start_end(2);
        assert(end_ <= n_frames); 

    else 
        start = 1; 
        end_ = n_frames; 
    end 

    for ii = start:end_
        % Read the target frame from the video 
        frame = video_reader.readFrame('frameNum', ii, ...
                                        'color', options.color,... 
                                        'zeros_as_nans', false...
                                       ); 
        
        % Sum the target channels 
        % applying the desired transformation
        for ch = 1:numel(options.channels)
            target_channel = options.channels(ch);
            channel_data = frame(:, :, target_channel);
            global_sums(ch) = global_sums(ch) + sum( transformation(channel_data(:)), 'omitnan');                
        end 
        assert(~any(isnan(global_sums(:) ))); 

        % ---- UPDATE PROGRESS BAR ----
        if options.verbose
            waitbar(ii / n_frames, h_wait, ...
                sprintf('Calculating Channel Means | Processing frame %d / %d', ii, n_frames));
        end
    
    end

    % ---- CLOSE PROGRESS BAR ----
    if options.verbose
        close(h_wait);
    end

    % The number of total values is num_frames * (IMAGE_SIZE)
    global_means = global_sums; 
    total_num_values = n_frames * video_reader.Height * video_reader.Width; 
    for ch = 1:numel(options.channels)
        global_means(ch) = global_sums(ch) / total_num_values;
    end     
    assert(~any(isnan(global_means(:) ))); 

    return 

end 