function perspective_transform = calculate_perspective_transform_w2e(world_camera_intrinsics, target_pos_ang_intended, target_pos_screen)

    arguments 
        world_camera_intrinsics; % Struct containing the world camera intrinsics
        target_pos_ang_intended; % Position of gaze targets in degrees of visual angle 
        target_pos_screen; % Position of gaze targets on the sensor screen 
    end 

    % First, if we do not already know the target positions on screen, we will calculate them 
    if(isempty(target_pos_screen))
        error("Not yet implemented"); 
    end 

    % Ensure we have the same number of pos in angles and screen 
    assert(size(target_pos_ang_intended, 1) == size(target_pos_screen, 1));  

    % Convert the screen pixel coordinates into [N x azimuth x elevation]
    pixel_azimuth_elevation = anglesFromIntrinsics(target_pos_screen, world_camera_intrinsics);
    geometric_transform = fitgeotform2d( pixel_azimuth_elevation, target_pos_ang_intended, 'projective');

    % Transform the gaze targets as seen by the camera into the eye rotation
    % coordinate space
    target_pos_angles_measured = transformPointsForward( geometric_transform, pixel_azimuth_elevation );

    % initialize return struct
    perspective_transform.data.world_camera_intrinsics = world_camera_intrinsics; 
    perspective_transform.data.target_pos_ang_intended = target_pos_ang_intended; 
    perspective_transform.data.target_pos_screen = target_pos_screen; 
    perspective_transform.fit.target_pos_angles_measured = target_pos_angles_measured; 
    perspective_transform.fit.geometric_transform = geometric_transform; 


end     