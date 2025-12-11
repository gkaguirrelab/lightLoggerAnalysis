function frame = readFrame(obj, options)
    arguments 
        obj 
        options.frameNum {mustBeNumeric} = []; 
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY"])} = "RGB";
        options.zeros_as_nans = false; 
        options.verbose = true; 
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

    % Determine if we need to do the conversion to HDF5 file 
    % when the first frame is read, or when desired color mode changes, or when the buffer needs to be re-read
    % by going under or over the range of the buffer 
    if(~exist(obj.temporary_reading_hdf5_filepath, "file")...
       || ~strcmp(obj.current_reading_color_mode, options.color)...
       || frameNum > ( (obj.buffer_start_frame + obj.read_ahead_buffer_size) - 1)...
       || frameNum < obj.buffer_start_frame...
       )

        % Convert to HDF5 
        obj.utility_library.video_to_hdf5(obj.full_video_path, obj.temporary_reading_hdf5_filepath,...
                                          options.color,...
                                          py.int(0), py.float("inf"),...
                                          options.zeros_as_nans,...
                                          options.verbose...
                                         )

        % If the read ahead buffer was not iniitalized, initialize it 
        % This let's us maintain a buffer of N frames in memory so we 
        % do not need to access the file for each individual frame 
        obj.buffer_start_frame = frameNum; 
        remaining = obj.NumFrames - frameNum + 1;
        
        % If somehow there are no more frames in the dataset, then simply return nans
        if(remaining == 0)
            if(options.color == "GRAY")
                frame = nan(obj.Height, obj.Width);
            else 
                frame = nan(obj.Height, obj.Width, 3);
            end 

            return ;  
        end

        if(options.color == "GRAY")
            start = [1, 1, frameNum];

            block_size = double(min(remaining, obj.read_ahead_buffer_size));
            count = [Inf, Inf, block_size];

            obj.read_ahead_buffer = (h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count));
            obj.read_ahead_buffer = permute(obj.read_ahead_buffer, [3 2 1]);
        else
            start = [1, 1, frameNum, 1];
            block_size = double(min(remaining, obj.read_ahead_buffer_size));
            count = [Inf, Inf, block_size, 3]; 
            
            obj.read_ahead_buffer = (h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count));
            obj.read_ahead_buffer = permute(obj.read_ahead_buffer, [3 2 1 4]);
        end 

        if(options.zeros_as_nans)
            obj.read_ahead_buffer(obj.read_ahead_buffer == 0) = nan; 
        end 

        % Garbage collect from Python and save the current color/buffer position
        py.gc.collect();                                          
        obj.current_reading_color_mode = options.color; 
    end 

    % Retrieve the frame from the buffer of frames we have read 
    buffer_position = (frameNum - obj.buffer_start_frame) + 1; 
    frame = squeeze(obj.read_ahead_buffer(buffer_position, :, :, :)); 

    % Ensure frame size somehow has not chnanged
    frame_size = size(frame);
    if(frame_size(1) ~= obj.Height || frame_size(2) ~= obj.Width)
        disp("Metadata: "); 
        disp([obj.Height, obj.Width]);

        disp("Read Frame");
        disp(frame_size);
        error("Frame read of inhomogenous shape to metadata");
    end     

    % If we are in grayscale and we want to convert zeros to nans, do it now 
    if(options.zeros_as_nans)
        frame(frame == 0) = nan; 
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
