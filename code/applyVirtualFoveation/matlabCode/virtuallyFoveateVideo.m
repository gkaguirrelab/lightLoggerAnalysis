function virtuallyFoveateVideo(world_video, sensor_t_cell, gaze_angles, gaze_offsets, blnk_events, output_path, path_to_intrinsics, options)
% Generate a virtually foveated video from a world-camera recording
%
% Syntax:
%   virtuallyFoveateVideo(world_video, sensor_t_cell, gaze_angles, gaze_offsets, ...
%       blnk_events, output_path, path_to_intrinsics)
%   virtuallyFoveateVideo(world_video, sensor_t_cell, gaze_angles, gaze_offsets, ...
%       blnk_events, output_path, path_to_intrinsics, options)
%
% Description:
%   This function generates a virtually foveated output video from a
%   world-camera recording by aligning each world frame to the nearest
%   available gaze sample in time and resampling the frame into a
%   gaze-centered retinal representation. Gaze samples are adjusted by a
%   participant-specific constant offset and an optional manual offset
%   before being used for virtual foveation.
%
%   The function processes the video in chunks to reduce memory pressure.
%   For each frame in a chunk, it finds the nearest pupil timestamp,
%   retrieves the corresponding gaze angle, checks whether the gaze is
%   invalid or falls during a blink event, and if valid, calls
%   virtuallyFoveateFrame to produce a gaze-centered output frame. Frames
%   with invalid gaze, blink events, or empty source data are written as
%   blank frames in the output video.
%
%   Frame processing within each chunk is parallelized with parfor. The
%   output video is written to disk incrementally after each chunk is
%   processed.
%
% Inputs:
%   world_video            - Char/string. Path to the playable world-camera
%                            video to be processed.
%   sensor_t_cell          - Cell array. Two-element cell array containing:
%                                {world_t, pupil_t}
%                            where world_t gives the timestamps for the
%                            world frames and pupil_t gives the timestamps
%                            for the gaze samples.
%   gaze_angles            - Numeric matrix. Matrix of gaze angles with one
%                            row per pupil sample. The first two columns
%                            are interpreted as gaze azimuth and elevation.
%   gaze_offsets           - Numeric vector. Constant participant-specific
%                            gaze offset to subtract from the gaze angles.
%   blnk_events            - Numeric matrix. Blink event table with one row
%                            per blink and two columns giving blink start
%                            and end timestamps.
%   output_path            - Char/string. Path where the virtually
%                            foveated output video will be written.
%   path_to_intrinsics     - Char/string. Path to the camera intrinsics
%                            calibration file used by virtuallyFoveateFrame.
%
% Optional key/value pairs:
%   world_fps              - Scalar. Frame rate assigned to the output
%                            world video writer.
%   pupil_fps              - Scalar. Sampling rate of the pupil camera.
%                            Included for metadata/compatibility.
%   nan_deg_threshold      - Scalar. Absolute gaze-angle threshold in
%                            degrees above which gaze samples are treated
%                            as invalid and replaced with blank frames.
%   frames_to_process      - Two-element numeric vector. Start and end
%                            world-frame indices to process. Use [1, inf]
%                            to process the full video.
%   verbose                - Logical. If true, print progress information
%                            during video reading, processing, and writing.
%   manual_offset          - Two-element numeric vector. Additional manual
%                            offset applied after subtracting the
%                            participant-specific gaze offset.
%   video_read_cache_size  - Scalar. Number of frames to process per chunk
%                            and read-ahead buffer size for video I/O.
%
% Outputs:
%   none
%
% Examples:
%{
    sensor_t_cell = {world_t, pupil_t};

    virtuallyFoveateVideo( ...
        "/path/to/W.avi", ...
        sensor_t_cell, ...
        gaze_angles, ...
        [0 0], ...
        blnk_events, ...
        "/path/to/output.avi", ...
        "/path/to/intrinsics_calibration.mat", ...
        "frames_to_process", [1 inf], ...
        "manual_offset", [0 0], ...
        "verbose", true ...
    );
%}

    arguments
        world_video {mustBeText}
        sensor_t_cell
        gaze_angles {mustBeMatrix}
        gaze_offsets {mustBeNumeric}
        blnk_events
        output_path {mustBeText}
        path_to_intrinsics {mustBeText}
        
        options.world_fps = 120
        options.pupil_fps {mustBeNumeric} = 120
        options.nan_deg_threshold = 45
        options.frames_to_process = [1, inf]
        options.verbose = false
        options.manual_offset = [0, 0]
        options.video_read_cache_size = 1000
    end
    % Import the Python util library 
    virutal_foveation_util = import_pyfile(getpref("lightLoggerAnalysis", "virtual_foveation_util_path"));

    % Create a video IO reader wrapper we will use to read/write into the original video
    if(options.verbose)
        disp("Opening video reader/writer")
    end 
    world_frame_reader = videoIOWrapper(world_video,... 
                                        "ioAction", 'read', ...
                                        "readAheadBufferSize", options.video_read_cache_size...
                                        ); 
    world_frame_writer = videoIOWrapper(output_path, "ioAction", 'write'); 
    world_frame_writer.FrameRate = options.world_fps; 

    % Apply the gaze offsets to the gaze angles, and adjust their coordinate system 
    gaze_angles_original = gaze_angles(:, 1:2) - gaze_offsets; 

    % We first subtract the constant gaze offset from the gaze angles (measured once per pariticpant)
    % Then, we flip the signs to be upsidedown and left handed. Then, 
    % we apply a manual offset from the April Tag. The sign of this corresponds to the follow
    % +azi = move right, +ele = move up
    gaze_angles(:, 1:2) = ( ( gaze_angles(:, 1:2) - gaze_offsets ) .* [1, -1    ] ) + options.manual_offset;

    % Extract world and pupil t 
    world_t = sensor_t_cell{1};
    pupil_t = sensor_t_cell{2};

    % Choose the bounds for our virtual foveation 
    start_frame = options.frames_to_process(1); 
    end_frame = world_frame_reader.NumFrames; 
    if(options.frames_to_process(2) ~= inf)
        end_frame = options.frames_to_process(2);
    end 
    assert(start_frame > 0 && end_frame <= world_frame_reader.NumFrames);

    % Open the writer to start writing frames 
    open(world_frame_writer); 

    % Iterate over the world frames 
    if(options.verbose)
        disp("Beginning frame processing")
    end 
    
    % First, let's find the number of chunks in the portion of the
    % video selected 
    num_chunks = ceil( ( (end_frame - start_frame) + 1) / options.video_read_cache_size);
    current_start = start_frame; 

    % Next, we will iterate over the chunks of the video 
    total_time = tic; 
    for cc = 1:num_chunks
        chunk_time = tic; 
        
        % If verbose, print what chunk we are on
        if(options.verbose)
            fprintf("Processing chunk: %d/%d\n", cc, num_chunks);
        end  


        % Find the end frame for this chunk 
        current_end = min((current_start + options.video_read_cache_size) - 1, end_frame); 

        % Initialize container for the frames to be read in 
        chunk_size = (current_end - current_start) + 1; 
        chunk_frames = zeros(chunk_size, world_frame_reader.Height, world_frame_reader.Width, 3, 'uint8'); 

        % Initialize output buffer that will write to the disk 
        output_buffer = zeros(chunk_size, 480, 480, 3, 'uint8');
        desiredN = size(output_buffer, 2); 

        dq = parallel.pool.DataQueue;
        h = waitbar(0, 'Processing frames...');

        counter = java.util.concurrent.atomic.AtomicInteger(0);
        afterEach(dq, @(~) waitbar(counter.incrementAndGet() / chunk_size, h));

        % Next, let's read in the frames for the chunk 
        input_buffer_insertion_idx = 1; 
        for ff = current_start : current_end; 
            % Force rebuffer every time a chunk starts (should not need to do this, but just for safety)
            force_rebuffer = input_buffer_insertion_idx == 1; 
            
            % Read the frame and insert it into the chunk frames buffer 
            chunk_frames(input_buffer_insertion_idx, :, :, :) = world_frame_reader.readFrame('frameNum', ff,...
                                                                  'color', 'RGB',...
                                                                  'verbose', false,...
                                                                  'force_rebuffer', force_rebuffer...
                                                                ); 

            % Update the insertion index 
            input_buffer_insertion_idx = input_buffer_insertion_idx + 1; 
        end  

        % Now, we read in all the frames. We will now process them in parallel
        parfor pp = 1:size(chunk_frames, 1)
            % Find the global world frame number 
            % that corresponds to this frame 
            global_world_frame_number = (pp + current_start) - 1; 

            % Retrieve the world frame timestamp and frame 
            world_timestamp = world_t(global_world_frame_number); 
            world_frame = squeeze(chunk_frames(pp, :, :, :)); 
            
            % Find the gaze angle that corresponds to this frame 
            [~, gaze_angle_idx] = min(abs(pupil_t - world_timestamp));
            gaze_angle = gaze_angles(gaze_angle_idx, 1:2); 

            % NaN the gaze angle if it's above a large threshold
            if( any(abs(gaze_angle) > options.nan_deg_threshold) )
                gaze_angle(:) = nan;
            end
            
            % If any of the cases for which we should skip a frame are true, just skip and we 
            % will output a blank frame already pre-allocated in the output buffer 
            if(~any(world_frame(:)) || any(isnan(gaze_angle)) || is_blnk_event(pupil_t(gaze_angle_idx), blnk_events))
                send(dq, 1);
                continue
            end 
            
            % Otherwise, virtually foveate the frame 
            virtually_foveated_frame = uint8(virtuallyFoveateFrame(world_frame, gaze_angle, ...
                                                                   path_to_intrinsics, 'desiredN', desiredN...
                                                                   )...
                                            );
            output_buffer(pp, :, :, :) = virtually_foveated_frame; 

            % Update the wait bar 
            send(dq,1);
        end 

        close(h);

       % Now that all these frames have been processed in parallel, let's write them out to the disk
        for oo = 1:size(output_buffer, 1)
            frame_to_write = squeeze(output_buffer(oo, :, :, :)); 
            world_frame_writer.writeVideo(frame_to_write); 
        end 

        % Update the start for the next chunk 
        current_start = current_start + chunk_size;
        
        chunk_elapsed_seconds = toc(chunk_time);
        fprintf("Chunk elpased time: %f seconds\n", chunk_elapsed_seconds);

    end 

    % Close the world video writer 
    close(world_frame_writer); 

    elapsed_seconds = toc(total_time); 
    fprintf("Total Elapsed Time: %f seconds\n", elapsed_seconds); 

    return; 

end 


% Local function to determine if a given pupil frame timestamp 
% is in a range of blink events 
function is_blink = is_blnk_event(timestamp, blnk_events)
    is_blink = false; 

    % Iterate over the blnk_events
    for rr = 1:size(blnk_events, 1)
        row = blnk_events(rr, :);
        event_start = row(1);
        event_end = row(2);
        
        % If the current event end time is before the current event, we can just skip 
        if(event_end < timestamp)
            continue; 
        end 

        % If the current event start time is after the timestamp 
        % we are searching for, just return early. No event after this 
        % could be a range where this timestamp lies 
        if(event_start > timestamp)
            return; 
        end 

        % If the timestamp is in this range, it is a BLNK event, so 
        % return true 
        if((event_start <= timestamp) && (event_end >= timestamp))
            is_blink = true;
            return ; 
        end 
    end 


    return ; 

end