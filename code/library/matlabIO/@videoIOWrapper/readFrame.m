function frame = readFrame(obj, options)
    arguments 
        obj 
        options.frameNum {mustBeNumeric} = []; 
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY", "LMS", "L+M+S", "L-M", "a", "c_lm", "c_s"])} = "RGB";
        options.zeros_as_nans = false; 
        options.verbose = false; 
        options.force_rebuffer = false; 
    end 
    verbose = options.verbose; 

    % Save a list of the LMS colorspace types  use later 
    persistent LMS_colorspace_types two_dimensional_colorspace_types normalized_color_types
    if(isempty(LMS_colorspace_types))
        LMS_colorspace_types = ["LMS", "L+M+S", "L-M", "a", "c_lm", "c_s"];
        two_dimensional_colorspace_types = ["GRAY", "L+M+S", "L-M", "a", "c_lm", "c_s"];
        normalized_color_types = ["a", "c_lm", "c_s"]; 
    end 

    if(options.verbose)
        tic; 
    end     

    % If no frame num was provided, we simply read the last frame read + 1 
    was_frame_passed = ~isempty(options.frameNum); 
    if(was_frame_passed)
        frameNum = options.frameNum; 
    else
        frameNum = obj.last_frame_read + 1; 
    end 

    % If somehow there are no more frames in the dataset, then raise an error 
    remaining = obj.NumFrames - frameNum + 1;
    if(remaining <= 0)
        error(sprintf("Frame number %d is out of bounds for video with NumFrames %d\n", frameNum, obj.NumFrames));
        return ;  
    end

    % If no color has been saved, save this current color 
    if(isempty(obj.current_reading_color_mode))
        obj.current_reading_color_mode = options.color; 
    end  

    % Determine if we need to do the conversion to HDF5 file 
    % when the first frame is read, or when desired color mode changes, or when the buffer needs to be re-read
    % by going under or over the range of the buffer 
    no_buffer_exists = ~exist(obj.temporary_reading_hdf5_filepath, "file"); 

    % If no buffer exists, set the buffer start frame to the current frameNum 
    if(no_buffer_exists)
        obj.buffer_start_frame = frameNum; 
    end 

    color_changed = string(obj.current_reading_color_mode) ~= string(options.color); 
    overran_buffer = frameNum > ( (obj.buffer_start_frame + obj.read_ahead_buffer_size) - 1); 
    underran_buffer = frameNum < obj.buffer_start_frame; 

    % Persistent normalization factors used with normalized color types. 
    % This should be recalculated when they are used 
    persistent normalization_factors

    should_rebuffer = any([no_buffer_exists, color_changed, overran_buffer, underran_buffer, options.force_rebuffer]);
    if(should_rebuffer)
        if(verbose)
            disp("REBUFFERING");
        end 

        % let's define the start and end indices of the read buffer 
        block_size = double(min(remaining, obj.read_ahead_buffer_size));

        % Because we will do MATLAB conversion of RGB -> LMS, 
        % if the desired color mode is LMS, we need to send RGB to Python and then convert the buffer in MATLAB 
        python_color = options.color; 
        if(ismember(options.color, LMS_colorspace_types))
            python_color = "RGB";
        end     

        % Convert a buffer of the video to HDF5 to be able to be read in 
        if(verbose)
            disp("CONVERTING TO HDF5"); 
        end     
        obj.utility_library.video_to_hdf5(obj.full_video_path, obj.temporary_reading_hdf5_filepath,...
                                          python_color,...
                                          py.int(frameNum - 1), py.int((frameNum - 1) + block_size),...
                                          options.zeros_as_nans,...
                                          options.verbose...
                                         )


        % If the read ahead buffer was not iniitalized, initialize it 
        % This let's us maintain a buffer of N frames in memory so we 
        % do not need to access the file for each individual frame 
        obj.buffer_start_frame = frameNum; 

        if(verbose)
            disp("READING FROM HDF5")
        end 
        info = h5info(obj.temporary_reading_hdf5_filepath, "/video");
        % Read gray frames from the temp file 
        read_ahead_buffer = []; 
        if(options.color == "GRAY")
            start = [1, 1, 1];
            count = [Inf, Inf, block_size];

            read_ahead_buffer = (h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count));
            read_ahead_buffer = permute(read_ahead_buffer, [3 2 1]);

        % Read color channel frames from the temp file
        else
            start = [1, 1, 1, 1];
            count = [3, Inf, Inf, block_size];  
            
            read_ahead_buffer = h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count);
            read_ahead_buffer = permute(read_ahead_buffer, [4 3 2 1]);
        end 

        % The buffer is now in memory. If we want to convert to LMS, we should do that now
        % so the whole buffer is in LMS space
        if(ismember(options.color, LMS_colorspace_types))
            if(verbose)
                disp("converting to LMS"); 
            end 

            % Retrieve the transformation matrices used to convert to LMS 
            % and the camera choice
            T_camera = obj.T_camera; 
            T_receptors = obj.T_receptors;  
            camera_used = obj.camera_used; 

            % Iterate over the buffer 
            parfor pp = 1:size(read_ahead_buffer, 1)
                read_ahead_buffer(pp, :, :, :) = rgb2lms(squeeze(read_ahead_buffer(pp, :, :, :)), T_receptors, T_camera,...
                                                         "camera", camera_used...
                                                        );
            end     
            
            
            % We first check if we are not in the normalized space 
            % This is a more simple calculation 
            % that does not require knowing certain normalization values
            if(~ismember(options.color, normalized_color_types))
                % If L+M+S, sum the channels together 
                % Leaving us a 2D image 
                if(options.color == "L+M+S")
                    read_ahead_buffer = sum(read_ahead_buffer, 4); 
                

                % If L-M, we need to subtract 2 channels and remove the third 
                % This leaves us a 2D image
                elseif (options.color == "L-M")
                    L = read_ahead_buffer(:,:,:,1);
                    M = read_ahead_buffer(:,:,:,2);
                    read_ahead_buffer = squeeze(L - M);
                
                end 
            
            else

                % The below are sourced from paper: Processing of Natural Temporal Stimuli by Macaque Retinal Ganglion Cells (J. H. van Hateren,1 L. Ru¨ ttiger,2,3 H. Sun,4 and B. B. Lee2,)

                % they first perform the following transformation 

                %% ------------------------------------------------------------------------
                % Cone signal preprocessing (contrast normalization)
                %
                % Step 1: Log transform + mean subtraction
                %
                % For each cone channel (L, M, S):
                %
                %   i_hat = log(i) - mean(log(i))
                %          
                %
                %
                % Interpretation:
                % - Log transform converts multiplicative changes → additive (intensity → contrast)
                % - Subtracting mean centers the signal over time
                %
                % Result:
                % - i_hat > 0  → brighter than average
                % - i_hat < 0  → darker than average
                %
                % Same applied to:
                %   l_hat, m_hat, s_hat
                %
                % i_hat represents deviation from average in contrast units
                %% ------------------------------------------------------------------------
                
                % If we are calculating any of the normalized 
                % color types for the first time, we need 
                % to calculate the normalization factors 
                % for each of l, m, s
                if(isempty(normalization_factors))
                    if(verbose)
                        disp("Calculating normalization factors")
                    end 

                    % Calculate the normalization factors 
                    % using the Python helper function 
                    normalization_factors = calculate_global_channels_mean(obj.full_video_path, ...
                                                                            "color", "LMS",...
                                                                            "sum_of", "loge",...
                                                                            "channels", [1, 2, 3],...
                                                                            "verbose", verbose...
                                                                            ); 

                end     


                % Let's calculate the normalized buffer
                l_hat = log(read_ahead_buffer(:, :, :, 1)) - normalization_factors(1); 
                m_hat = log(read_ahead_buffer(:, :, :, 2)) - normalization_factors(2);
                s_hat = log(read_ahead_buffer(:, :, :, 3)) - normalization_factors(3);

                % %% ------------------------------------------------------------------------
                % Achromatic signal (overall brightness)
                %
                % Combine L and M cone contrast signals:
                %
                %   a = (l_hat + m_hat) / sqrt(2)
                %
                % Interpretation:
                % - Represents overall brightness (luminance-like signal)
                % - S cone is ignored (minor contribution to luminance in primates)
                %
                % Rough meaning:
                %   "How bright is the scene relative to its average?"
                %% ------------------------------------------------------------------------
                if(options.color == "a")
                    read_ahead_buffer = (l_hat + m_hat) / sqrt(2); 

                %% ------------------------------------------------------------------------
                % Chromatic signal 1 (red–green opponent channel)
                %
                % Difference between L and M cones:
                %
                %   c_lm = (l_hat - m_hat) / sqrt(2)
                %
                % Interpretation:
                % - Measures red vs green balance
                %
                % Sign:
                %   c_lm > 0  → more L → reddish
                %   c_lm < 0  → more M → greenish
                %% ------------------------------------------------------------------------
                elseif(options.color == "c_lm")
                    read_ahead_buffer = (l_hat - m_hat) / sqrt(2); 
                
                %% ------------------------------------------------------------------------
                % Chromatic signal 2 (blue–yellow opponent channel)
                %
                % Compare S cone vs combined L+M:
                %
                %   c_s = (2*s_hat - (l_hat + m_hat)) / sqrt(6)
                %
                % Interpretation:
                % - Compares blue (S) vs yellow (L+M)
                %
                % Sign:
                %   c_s > 0  → blue-ish
                %   c_s < 0  → yellow-ish
                %% ------------------------------------------------------------------------
                elseif(options.color == "c_s")
                    read_ahead_buffer = (2 * s_hat - (l_hat + m_hat)) / sqrt(6); 

                end
            
            end

        end 
        
        % Save the read ahead buffer to the obj 
        obj.read_ahead_buffer = read_ahead_buffer;

        % Garbage collect from Python and save the current color/buffer position
        py.gc.collect();                                          
        obj.current_reading_color_mode = options.color; 
    end 

    % Retrieve the frame from the buffer of frames we have read 
    buffer_position = (frameNum - obj.buffer_start_frame) + 1; 
    frame = squeeze(obj.read_ahead_buffer(buffer_position, :, :, :)); 
    if(ismember(obj.current_reading_color_mode, two_dimensional_colorspace_types))
        frame = reshape(frame, obj.Height, obj.Width); 
    else
        frame = reshape(frame, obj.Height, obj.Width, 3); 

    end 

    % Ensure frame size somehow has not chnanged
    frame_size = size(frame);
    if(frame_size(1) ~= obj.Height || frame_size(2) ~= obj.Width)
        disp("Metadata: "); 
        disp([obj.Height, obj.Width]);

        disp("Read Frame");
        disp(frame_size);
        error("Frame read of inhomogenous shape to metadata");
    end     

    % If we are converting zeros to NaNs
    if (options.zeros_as_nans)

        if (ndims(frame) == 2 || size(frame,3) == 1)
            % 2 Channel: replace all 0s with NaN
            frame(frame == 0) = nan;

        elseif (ndims(frame) == 3 && size(frame,3) == 3)
            % 3 channel: find pixels where ALL channels are zero
            zeroMask = (frame(:,:,1) == 0) & ...
                    (frame(:,:,2) == 0) & ...
                    (frame(:,:,3) == 0);

            % Apply mask to all channels
            for c = 1:3
                tmp = frame(:,:,c);
                tmp(zeroMask) = nan;
                frame(:,:,c) = tmp;
            end
        end

    end
    % Ensure we are not returning an empty variable 
    if(isempty(frame))
        error("readFrame is about to return an empty variable");
    end 

    % Update the last read frame 
    obj.last_frame_read = frameNum; 

    if(options.verbose)
        fprintf("Frame read took: %f seconds\n", toc);
    end 

    return ; 

end
