function [x_offset, y_offset] = calculateManualOffset(subjectID, activity, options)
    arguments
        subjectID; 
        activity;
        options.april_tag_frames {mustBeNumeric} = [];
    end     

    % If not provided april tag frames, let's calculate them 
    if(numel(options.april_tag_frames) == 0)
        error("Must pass april tag frames. Auto finding is not yet implemented")
    end 

    % Generate the video that has just these april tag frames 
    generateVirtuallyFoveatedVideos({subjectID}, {[1, inf]}, "activity", activity, ...
                                    "non_contiguous_target_frames", {options.april_tag_frames},...
                                    "video_read_cache_size", 10000 ...
                                   );

    % Define the path to this video and read its frames  
    temp_output_path = fullfile(".", sprintf("/FLIC_%s_%s_%s.avi", subjectID, activity, "virtuallyFoveatedVideoUnspecified"));
    video_reader = videoIOWrapper(temp_output_path, "ioAction", 'read');
    frames = nan(numel(options.april_tag_frames), video_reader.Height, video_reader.Width);
    if(video_reader.NumFrames ~= size(frames, 1))
        fprintf("Num frames in video: %d\n", video_reader.NumFrames);
        fprintf("Expected; %d\n", size(frames, 1));
        error("Unequal number of frames in video and april tag frames");
    end 

    for frame_num = 1:size(frames, 1)   
        frame = video_reader.readFrame("frameNum", frame_num, "color", "GRAY", "zeros_as_nans", true);
        if(size(frame) ~= [video_reader.Height, video_reader.Width])
            disp("Frame size")
            disp(size(frame))
            disp("Video reader metadata")
            disp([video_reader.Height, video_reader.Width])
            error("Frame size does not equal video reader metadata")
        end 
    end 

    % Pass the frames to Python to calculate the offset from the frames 
    virtual_foveation_util = import_pyfile(getpref("lightLoggerAnalysis", "virtual_foveation_util_path"));
    manual_offsets = double(virtual_foveation_util.calculate_manual_offsets(frames)); 
    x_offset_raw = manual_offsets(1);
    y_offset_raw = manual_offsets(2); 
    
    % Next, we need to convert the offsets to the correct coordinate space. 
    % + is left for x and + is up for y 
    x_offset = -x_offset_raw; 
    y_offset = y_offset_raw;

    % Remove the temp file 
    delete(temp_output_path);

    % Garbage collect from Python and return 
    py.gc.collect(); 

    return; 
end 


