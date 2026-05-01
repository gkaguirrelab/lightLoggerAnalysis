function frame = readFrame(obj, options)
    arguments 
        obj 
        options.frameNum {mustBeNumeric} = []; 
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY", "LMS", "L+M+S", "L-M", "a", "c_lm", "c_s"])} = "RGB";
        options.zeros_as_nans = false; 
        options.ceiling_as_nans = false;
        options.verbose = false; 
        options.force_rebuffer = false; 
        options.parallel_buffer = true; 
    end 
    
    verbose = options.verbose; 

    % Save a list of the LMS colorspace types  use later 
    persistent LMS_colorspace_types two_dimensional_colorspace_types normalized_color_types
    if(isempty(LMS_colorspace_types))
        LMS_colorspace_types = ["LMS", "L+M", "L+M+S", "L-M", "S", "a", "c_lm", "c_s"];
        two_dimensional_colorspace_types = ["GRAY", "L+M+S", "L-M", "L+M", "S", "a", "c_lm", "c_s"];
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

    % When converting to LMS we must ensure ceiling as NaN 
    % and zeros_as_nans are true 
    if(ismember(obj.current_reading_color_mode, LMS_colorspace_types))
        assert(options.zeros_as_nans && options.ceiling_as_nans, "When converting to LMS types, both ceiling and floor must be NaN"); 
    end 


    % If any rebuffer condition has been met, let's rebuffer
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

        % Python reading ONLY supports format of GRAY, RGB, BGR. Therefore, values should be 0-255
        % once they come out of here at MAXIMUM.
        % if floor and ceiling is desired to be NaNs, then 0 and 255 should come out as NaNs
        obj.utility_library.video_to_hdf5(obj.full_video_path, obj.temporary_reading_hdf5_filepath,...
                                          python_color,...
                                          py.int(frameNum - 1), py.int((frameNum - 1) + block_size),...
                                          options.zeros_as_nans,...
                                          options.ceiling_as_nans,...
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

        % Guard against bad values immediately after reading in the buffer 
        if(options.zeros_as_nans)
            assert(~any(read_ahead_buffer(:) <= 0), "0 was requested to be NaNs but values <= 0 were detected in the buffer after it has been read in"); 
        end 

        if(options.ceiling_as_nans)
            assert(~any(read_ahead_buffer(:) >= 255), "255 was requested to be NaNs but values >= 255 were detected in the buffer after it has been read in");
        end 

        % Display a progress bar if desired 
        % for in-MATLAB processing
        dq = [];
        if(options.verbose)
            % Create a waitbar for the parallel loop
            hWait = waitbar(0, 'Converting RGB frames to LMS...');

            % Progress counter
            nFrames = size(read_ahead_buffer, 1);
            progressCount = 0;

            % DataQueue lets parfor workers notify the client process
            dq = parallel.pool.DataQueue;

            % Update waitbar each time a worker finishes one frame
            afterEach(dq, @updateWaitbar);
            
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

            % Iterate over the buffer either in parallel or sequentially 

            % Parallel buffering 
            if(options.parallel_buffer)

                parfor pp = 1:size(read_ahead_buffer, 1)
                    % Read the RGB frame
                    rgb_frame = squeeze(read_ahead_buffer(pp, :, :, :)); 
                    assert(~any(rgb_frame(:) <= 0), "RGB Frame has value of 0 instead of NaN when converting to LMS"); 
                    assert(~any(rgb_frame(:) >= 255), "RGB Frame has value of 255 instead of NaN when converting to LMS"); 
                    
                    % Convert to LMS
                    converted = rgb2lms(rgb_frame, T_receptors, T_camera, "camera", camera_used); 
                    read_ahead_buffer(pp, :, :, :) = converted;
                    
                    if(options.verbose)
                        % Send progress update back to client
                        send(dq, pp);
                    end 

                end     

                % Close waitbar when complete
                if(options.verbose)
                    close(hWait);
                end 


            % Sequential buffering 
            else 
                for pp = 1:size(read_ahead_buffer, 1)
                    % Read the RGB frame
                    rgb_frame = squeeze(read_ahead_buffer(pp, :, :, :)); 
                    assert(~any(rgb_frame(:) <= 0), "RGB Frame has value of 0 instead of NaN when converting to LMS"); 
                    assert(~any(rgb_frame(:) >= 255), "RGB Frame has value of 255 instead of NaN when converting to LMS"); 
                    
                    % Convert to LMS
                    converted = rgb2lms(rgb_frame, T_receptors, T_camera, "camera", camera_used); 
                    read_ahead_buffer(pp, :, :, :) = converted;

                    if(options.verbose)
                        % Send progress update back to client
                        send(dq, pp);
                    end 

                end 

                % Close waitbar when complete
                if(options.verbose)
                    close(hWait);
                end 
            end

            % We first check if we are not in the normalized space 
            % This is a more simple calculation 
            % that does not require knowing certain normalization values
            if(~ismember(options.color, normalized_color_types))
                % If L+M+S, sum the channels together 
                % Leaving us a 2D image 
                if(options.color == "L+M+S")
                    read_ahead_buffer = sum(read_ahead_buffer, 4); 
                
                % If L-M, we need to add 2 channels and remove the third 
                % This leaves us a 2D image
                elseif(options.color == "L+M")
                    L = read_ahead_buffer(:,:,:,1);
                    M = read_ahead_buffer(:,:,:,2);
                    read_ahead_buffer = squeeze(L + M);

                % If L-M, we need to subtract 2 channels and remove the third 
                % This leaves us a 2D image
                elseif (options.color == "L-M")
                    L = read_ahead_buffer(:,:,:,1);
                    M = read_ahead_buffer(:,:,:,2);
                    read_ahead_buffer = squeeze(L - M);
                
                % If S, then we need to simply extract the last channel of the image 
                elseif (options.color == "S")
                    S = read_ahead_buffer(:, :, :, 3);
                    read_ahead_buffer = squeeze(S);
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
                if(isempty(obj.normalization_factors))
                    if(verbose)
                        disp("Calculating normalization factors")
                    end 

                    % Calculate the normalization factors 
                    % using the Python helper function 
                    obj.normalization_factors = calculate_global_channels_mean(obj.full_video_path, ...
                                                                            "color", "LMS",...
                                                                            "sum_of", "loge",...
                                                                            "channels", [1, 2, 3],...
                                                                            "verbose", verbose,...
                                                                            "zeros_as_nans", options.zeros_as_nans...
                                                                            "ceiling_as_nans", options.ceiling_as_nans... 
                                                                            ); 

                end     


                % Let's calculate the normalized buffer
                l_hat = log(L) - obj.normalization_factors(1); 
                m_hat = log(M) - obj.normalization_factors(2);
                s_hat = log(S) - obj.normalization_factors(3);

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

    % Local function used to update progress waitbars
    function updateWaitbar(~)
        progressCount = progressCount + 1;

        if isvalid(hWait)
            waitbar(progressCount / nFrames, hWait, ...
                sprintf('Converting RGB frames to LMS... %d / %d', ...
                progressCount, nFrames));
        end
    end

end
