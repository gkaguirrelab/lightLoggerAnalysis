function [x_offset, y_offset] = calculateManualOffset(subjectID,  activity,... 
                                                      path_to_world_video, ...
                                                      gaze_angles, participant_offsets,...
                                                      path_to_recording_chunks, path_to_intrinsics, path_to_perspective_projection,...
                                                      options
                                                     )
    arguments
        video; 
        options.april_tag_frames {mustBeNumeric} = [];
    end     

    % If not provided april tag frames, let's calculate them 
    if(numel(april_tag_frames) == 0)
        error("Must pass april tag frames. Auto finding is not yet implemented")
    end 

    % Generate the video that has just these april tag frames 
    generateVirtuallyFoveatedVideos({subjectID}, {[1, inf]}, "activity", activity, "non_contiguous_target_frames", {options.april_tag_frames});

    % Define the path to this video and read its frames  
    temp_output_path = fullfile(".", sprintf("/%s_%s_%s.avi", subjectID, activity, "virtuallyFoveatedVideoUnspecified"));
    video_reader = videoIOWrapper(temp_output_path, "ioAction", 'read');
    frames = nan(numel(options.april_tag_frames, video_reader.Height, video_reader.Width));
    for frame_num = 1:size(frames, 1)   
        frames(frame_num, :, :) = video_reader.readFrame("frameNum", frame_num, "color", "GRAY", "zeros_as_nans", true);
    end 

    % Pass the frames to Python to calculate the offset from the frames 
    virtual_foveation_util = import_pyfile(getpref("lightLoggerAnalysis", "virtual_foveation_util_path"));
    offsets = virtual_foveation_util.calculate_manual_offsets(frames); 
    [x_offset_raw, y_offset_raw] = offsets; 
    
    % Next, we need to convert the offsets to the correct coordinate space. 
    % + is left for x and + is up for y 
    x_offset = -x_offset_raw; 
    y_offset = y_offset_raw;

    % Garbage collect from Python and return 
    py.gc.collect(); 

    return; 
end 


