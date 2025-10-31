% Iterate over the extracted frames 


for ii = 1:size(extracted_frames, 1)
    % Retrieve the frame and the gaze angle 
    frame = double(squeeze(extracted_frames(ii, :, :))); 
    gaze_angle = gaze_angles(ii, :); 
    

    % Calculate the transformed frame 
    %perspective_transform.fit.geometric_transform 
    transformed_frame = coordinateTransformFinal(background_frame, world_camera_intrinsics, perspective_transform.fit.geometric_transform, [-1*((gaze_angle(2)+5.4)-10),(gaze_angle(1)-0.4)]); 

    size(transformed_frame)
        
    % Show the resulting frame 
    %{
    figure; 
    imagesc(xq, yq, Vq);
    axis image;
    set(gca,'YDir','normal');
    colormap gray;
    colorbar;
    xlim([-80, 80]); 
     ylim([-80, 80])
    title(sprintf("Frame Num %d", ii));
    xlabel('Visual angle [deg]');
    ylabel('Visual angle [deg]');
    %}
end 