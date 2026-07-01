function global_means = calculate_global_channels_mean(video_path, options)
% Compute channel means across an entire video or a selected frame range.
%
% Syntax:
%   global_means = calculate_global_channels_mean(video_path, options)
%
% Description:
%   This function streams frames from a video through `videoIOWrapper`,
%   applies an optional per-pixel transformation such as the natural log,
%   and accumulates channelwise sums and valid-sample counts across the
%   requested frame interval. It is primarily used to estimate the global
%   normalization constants needed by the LMS-derived opponent-color
%   representations, but it also supports ordinary mean estimation in a
%   variety of color spaces.
%
% Inputs:
%   video_path               - String. Path to the video file to analyze.
%   options                  - Name/value options controlling the color
%                              space, transformation (`value` or `loge`),
%                              channels to average, frame subset, and
%                              handling of floor/ceiling values.
%
% Outputs:
%   global_means             - Column vector containing the mean of each
%                              requested channel after the specified
%                              transformation and NaN handling.
%
% Examples:
%{
    global_means = calculate_global_channels_mean("W.avi", ...
        "color", "LMS", "sum_of", "loge", "channels", [1 2 3]);
%}

    arguments 
        video_path string;         
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY", "LMS", "L+M+S", "L-M", "a", "c_lm", "c_s"])} = "RGB"; 
        options.sum_of {mustBeMember(options.sum_of, ["value", "loge"])} = "value"; 
        options.channels = [1, 2, 3]; 
        options.verbose = false 
        options.start_end = false; 
        options.zeros_as_nans = false; 
        options.ceiling_as_nans = false; 
        options.apply_floor_ceiling = false; 
    end 

    % First, let's open a video reader to stream frames from the video 
    video_reader = videoIOWrapper(video_path, 'ioAction', 'read'); 

    % Initialize the sum variable. It should be a 
    % vector consistening of the number of channels 
    % we would like to mean over 
    global_sums = zeros(numel(options.channels), 1); 

    % Initialuze the count of the number of nun NaN values for each channel. 
    total_non_nan = zeros(numel(options.channels), 1); 
    
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
            if(~options.zeros_as_nans)
                transformation = @(x) log(x + 10e-9); % epsiolon to fight inf
            
            else
                transformation = @(x) log(x); 
            end 

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
                                        'zeros_as_nans', options.zeros_as_nans,...
                                        'ceiling_as_nans', options.ceiling_as_nans,...
                                        'apply_floor_ceiling', options.apply_floor_ceiling...
                                       ); 
        
        % Sum the target channels 
        % applying the desired transformation
        for ch = 1:numel(options.channels)
            target_channel = options.channels(ch);
            channel_data = frame(:, :, target_channel);
        
            if(~options.zeros_as_nans && options.sum_of == "loge")
                assert(~any(frame(:) <= 0)); 
            end     

            global_sums(ch) = global_sums(ch) + sum( transformation(channel_data(:)), 'omitnan');                

            % Keep adding to the count of non NaN values by channel 
            total_non_nan(ch) = total_non_nan(ch) + sum(~isnan(channel_data(:))); 
        end 

        % If at ANY point these go to NaN, this is bad 
        assert(~any(isnan(global_sums))); 
        assert(~any(isnan(total_non_nan)));

        % If ANY point these go to INF this is also bad 
        assert(~any(isinf(global_sums)));
        assert(~any(isinf(total_non_nan)));

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

    % calculate the means 
    global_means = global_sums; 
    for ch = 1:numel(options.channels)
        total_num_values = total_non_nan(ch); 

        if(total_num_values == 0)
            global_means(ch) = 0; 
            continue; 
        end 

        global_means(ch) = global_sums(ch) / total_num_values;
    end     

    % Assert the final results are not nan and not INF 
    assert(~any(isnan(global_means(:)))); 
    assert(~any(isinf(global_means(:))));

    return 

end 
