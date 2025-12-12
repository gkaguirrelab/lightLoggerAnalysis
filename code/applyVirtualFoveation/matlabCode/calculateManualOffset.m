function [x_offset, y_offset] = calculateManualOffset(subjectID, activity, intrinsics, options)
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
        frames(frame_num, :, :) = frame; 
    end 

    % Click the targets to get the manual offsets
    manual_offsets_screen = click_points(frames)

    % Next, we need to convert it to degrees
    manual_offsets_deg = anglesFromIntrinsics(manual_offsets, intrinsics)
    x_offset_raw = manual_offsets_deg(1);
    y_offset_raw = manual_offsets_deg(2);

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

% Local function to do the clicking in screen coordinates for manual offsets
function mean_manual_offsets = click_points(frames)

    % Initialize output array for our offsets 
    nFrames = size(frames, 1);
    manual_offsets_per_target = nan(nFrames, 2);

    % Determine subplot grid size based on number of frames
    [rows, cols] = find_min_figsize(nFrames);

    figure;
    for ii = 1:nFrames
        % Show the image 
        subplot(rows, cols, ii);
        imshow(squeeze(frames(ii, :, :)), []);
        title(sprintf("Target: %d", ii))
        axis image off;
        
        % Wait for one click on the image 
        [x, y] = ginput(1);  
        manual_offsets_per_target(ii, :) = [x, y];
        
        % Plot the results of that click
        hold on;
        plot(x, y, 'r+', 'MarkerSize', 6, 'LineWidth', 1);
    end

    % Return the mean manual offset 
    mean_manual_offsets = mean(manual_offsets_per_target, 1, 'omitnan');

end