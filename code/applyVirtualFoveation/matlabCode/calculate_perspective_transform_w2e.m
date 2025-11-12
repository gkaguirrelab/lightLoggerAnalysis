function perspective_projection = calculate_perspective_transform_w2e(world_camera_intrinsics, target_pos_ang_intended, target_pos_screen)

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

    % initialize return struct
    perspective_projection.data.world_camera_intrinsics = world_camera_intrinsics; 
    perspective_projection.data.target_pos_ang_intended = target_pos_ang_intended; 
    perspective_projection.data.target_pos_screen = target_pos_screen; 
    perspective_projection.fit.target_pos_angles_measured = pixel_azimuth_elevation; 
    perspective_projection.fit.geometric_transform = geometric_transform; 

    figure; 
    t = tiledlayout(2, 2);
    
    nexttile; 
    title("Intended Target Angles");
    xlabel("Azimuth [deg]"); 
    ylabel("Elevation [deg]");
    hold on; 
    
    % Plot the intended gaze target positions and their number 
    for ii = 1:size(target_pos_ang_intended, 1)
        % Plot the target gaze angle
        plot(target_pos_ang_intended(ii, 1), target_pos_ang_intended(ii, 2), 'o', 'DisplayName', sprintf("%d", ii)); 

        % Add text label just above the marker
        text(target_pos_ang_intended(ii, 1), target_pos_ang_intended(ii, 2) + 0.8,... 
             sprintf('%d', ii), 'FontSize', 10, 'HorizontalAlignment', 'center'...
            );

    end     
    xlim([-30, 30]); 
    ylim([-30, 30]); 
    legend show; 

    nexttile; 
    title("Screen Coordinate Positions");
    xlabel("Col [px]"); 
    ylabel("Row [px]");
    hold on; 
    
    % Plot the intended gaze target positions and their number 
    for ii = 1:size(target_pos_screen, 1)
        % Plot the target gaze angle
        plot(target_pos_screen(ii, 1), target_pos_screen(ii, 2), 'o', 'DisplayName', sprintf("%d", ii)); 

        text(target_pos_screen(ii, 1), target_pos_screen(ii, 2) - 10,... 
             sprintf('%d', ii), 'FontSize', 10, 'HorizontalAlignment', 'center'...
            );

    end     
    set(gca, 'YDir', 'reverse');
    xlim([0, 700]); 
    ylim([0, 500]); 
    legend show; 



    nexttile; 
    title("Screen Coordinate as Degrees of Visual Angle");
    xlabel("Azimuth [deg]"); 
    ylabel("Elevation [deg]");
    hold on; 
    
    % Plot the intended gaze target positions and their number 
    for ii = 1:size(pixel_azimuth_elevation, 1)
        % Plot the target gaze angle
        plot(pixel_azimuth_elevation(ii, 1), pixel_azimuth_elevation(ii, 2), 'o', 'DisplayName', sprintf("%d", ii)); 

        text(pixel_azimuth_elevation(ii, 1), pixel_azimuth_elevation(ii, 2) + 0.8,... 
             sprintf('%d', ii), 'FontSize', 10, 'HorizontalAlignment', 'center'...
            );

    end     
    xlim([-30, 30]); 
    ylim([-30, 30]); 
    legend show; 


    eye_visual_angle_coords = transformPointsForward(geometric_transform, pixel_azimuth_elevation);
    nexttile; 
    title("Screen Coordinate as Eye Visual Angle");
    xlabel("Azimuth [deg]"); 
    ylabel("Elevation [deg]");
    hold on; 
    
    % Plot the intended gaze target positions and their number 
    for ii = 1:size(eye_visual_angle_coords, 1)
        % Plot the target gaze angle
        plot(target_pos_ang_intended(ii, 1), target_pos_ang_intended(ii, 2), 'o', 'DisplayName', sprintf("%d", ii)); 
        plot(eye_visual_angle_coords(ii, 1), eye_visual_angle_coords(ii, 2), 'x', 'DisplayName', sprintf("%d", ii)); 

    end     
    xlim([-30, 30]); 
    ylim([-30, 30]); 
    legend show; 




end     