function frame = readFrame(obj, options)
    arguments 
        obj 
        options.frameNum {mustBeNumeric} = []; 
        options.color {mustBeMember(options.color, ["RGB","BGR","GRAY"])} = "RGB";
        options.zeros_as_nans = false; 
        options.verbose = false; 
        options.force_rebuffer = false; 
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

    should_rebuffer = any([no_buffer_exists, color_changed, overran_buffer, underran_buffer, options.force_rebuffer]);
    if(should_rebuffer)
        % let's define the start and end indices of the read buffer 
        block_size = double(min(remaining, obj.read_ahead_buffer_size));

        % Convert a buffer of the video to HDF5 to be able to be read in 
        obj.utility_library.video_to_hdf5(obj.full_video_path, obj.temporary_reading_hdf5_filepath,...
                                          options.color,...
                                          py.int(frameNum - 1), py.int((frameNum - 1) + block_size),...
                                          options.zeros_as_nans,...
                                          options.verbose...
                                         )


        % If the read ahead buffer was not iniitalized, initialize it 
        % This let's us maintain a buffer of N frames in memory so we 
        % do not need to access the file for each individual frame 
        obj.buffer_start_frame = frameNum; 

        info = h5info(obj.temporary_reading_hdf5_filepath, "/video");
        if(options.color == "GRAY")
            start = [1, 1, 1];
            count = [Inf, Inf, block_size];

            obj.read_ahead_buffer = (h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count));
            obj.read_ahead_buffer = permute(obj.read_ahead_buffer, [3 2 1]);
        else
            start = [1, 1, 1, 1];
            count = [3, Inf, Inf, block_size];  
            
            obj.read_ahead_buffer = h5read(obj.temporary_reading_hdf5_filepath, "/video", start, count);
            obj.read_ahead_buffer = permute(obj.read_ahead_buffer, [4 3 2 1]);
        end 

        % Garbage collect from Python and save the current color/buffer position
        py.gc.collect();                                          
        obj.current_reading_color_mode = options.color; 
    end 

    % Retrieve the frame from the buffer of frames we have read 
    buffer_position = (frameNum - obj.buffer_start_frame) + 1; 
    if(obj.current_reading_color_mode == "GRAY")
        frame = reshape(obj.read_ahead_buffer(buffer_position, :, :, :), obj.Height, obj.Width); 
    else
        frame = reshape(obj.read_ahead_buffer(buffer_position, :, :, :), obj.Height, obj.Width, 3); 
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
