function [target_pos_angles_measured] = perspective_transform_w2e(world_camera_intrinsics, target_pos_ang, target_pos_screen, calibration_video_path)

    arguments 
        world_camera_intrinsics {mustBeStruct}; % Struct containing the world camera intrinsics
        target_pos_ang_intended {mustBeDouble}; % Position of gaze targets in degrees of visual angle 
        target_pos_screen {mustBeDouble}; % Position of gaze targets on the sensor screen 
        calibration_video_path {mustBeText}; % Path to the .avi calibration video used to find the screen target pos 
    end 

    % First, if we do not already know the target positions on screen, we will calculate them 
    if(isempty(target_pos_screen))
        error("Not yet implemented"); 
    end 

    % Ensure we have the same number of pos in angles and screen 
    assert(size(target_pos_ang_intended, 1) == size(target_pos_screen, 1));  

    % Convert the screen pixel coordinates into [N x azimuth x elevation]
    pixel_azimuth_elevation = anglesFromIntrinsics(target_position_screen, world_camera_intrinsics);
    geometric_transform = fitgeotform2d( pixel_azimuth_elevation, target_pos_ang_intended, 'projective');

    % Transform the gaze targets as seen by the camera into the eye rotation
    % coordinate space
    target_pos_angles_measured = transformPointsForward( geometric_transform, pixel_azimuth_elevation );


    % Get the set of gaze targets in units of visual angle as seen by the eye of the observer
    % If one assumes a distance (e.g. 1 meter) you can then convert any world camera image into a 3D curved surface using this final projection mapping.


end     