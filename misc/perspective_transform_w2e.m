function transformed_frames = perspective_transform_w2e(frames, fisheyeIntrinsicsPath, transformationPath, center_offsets)

    % Retrieve the number of frames to transform in frames 
    n_frames = size(frames,  1); 

    % Allocate output buffer 
    transformed_frames = zeros(size(frames)); 

    % Parallelize the transformation of frames
    disp("starting parallel")
    parfor ii = 1:n_frames
        % Extract the current frame 
        % and its associated offset 
        frame = squeeze(frames(ii, :, :)); 
        center_offset = squeeze(center_offsets(ii, :)); 

        % Transform the frame 
        transformed_frames(ii, :, :) = coordinateTransformFinal(frame, fisheyeIntrinsicsPath, transformationPath, center_offset)
    end 

    % Output the result 
    return ; 

end