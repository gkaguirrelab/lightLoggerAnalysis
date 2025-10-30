function virtually_foveated_frame = coordinateTransformFinal(I, gaze_angle, fisheyeIntrinsicsPath, transformationPath)
    % TODO: 
    %    center_offset is in degrees but not what you expect on the plot 
    %    first number is somehow y and then x. 
    %    also, left and up on the image is positive. 

    %


    % Need to make available to this function the fisheye intrinsics, and the
    % set of "imgPts" and "worldPts".
    % 
    % Example usage:
    %{
        %  ---- GET INTRINSICS ----
        data = load('/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/intrinsics_calibration/intrinsics_calibration.mat');
        fisheyeIntrinsics = camera_intrinsics_calibration.results.Intrinsics;  
        %  ---- GET VIDEO DATA ----
        final_AS
        %  ---- USE (example) ----
        myChoice = 'interceptMap';
        lightLevel = 'highAS';
        [SensorFigure, CameraFigure, EyeFigure] = coordinateTransformFinal(fisheyeIntrinsics, myChoice, lightLevel)
    %}
        
    
    figure; 
    plot(gaze_angle(1), gaze_angle(2), '+', 'MarkerSize', 12, 'DisplayName', 'Gaze Angle'); 
    hold on;  


    % --- 1. Load or Define Your Data ---
    % Load in the intrinsics and the transformation calculated we previously calculated to map TODO: XX to YY 
    fisheyeIntrinsics = load(fisheyeIntrinsicsPath).camera_intrinsics_calibration.results.Intrinsics; 
    transformation = load(transformationPath).perspective_transform.fit.geometric_transform; 
    
    % Get the camera visual field positions corresponding to positions of all
    % locations on the camera sensor
    [nRows, nCols]    = size(I);
    [xg, yg]          = meshgrid(1:nCols, 1:nRows);
    sensorPoints      = [xg(:),yg(:)];
    visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    
    % Transform the camera visual field points to eye rotation coordinate space
    % using the previously calculated tform
    eyeRotationCoordinates = transformPointsForward(transformation, visualFieldPoints);

    % START GEOFF CODE BLOCK

    % Extract x, y, v 
    x = eyeRotationCoordinates(:, 1);
    y = eyeRotationCoordinates(:, 2);
    v = I(:); 

    % --- 2. Define the Regular Grid for Interpolation ---
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);

    num_x_points = 600;
    num_y_points = 600;

    xi = linspace(xmin, xmax, num_x_points);
    yi = linspace(ymin, ymax, num_y_points);
    [XX, YY] = meshgrid(xi, yi);
    
    % --- 3. Interpolate the Scattered Data onto the Grid ---
    % The 'griddata' function performs the interpolation.
    % Method can be 'linear', 'cubic', 'nearest', or 'v4'.
    VI = griddata(x, y, v, XX, YY, 'linear');
    
        % --- 4. Apply circular mask around gaze angle ---
    radius_deg = 30; % Adjust this depending on your desired radius in degrees

    % Compute distance from each grid point to gaze angle
    dist = sqrt((XX - gaze_angle(1)).^2 + (YY - gaze_angle(2)).^2);

    % Create mask: keep pixels inside the radius
    mask = dist <= radius_deg;

    xi = xi - gaze_angle(1); 
    yi = yi - gaze_angle(2);

    % Apply mask (set outside region to NaN or 0)
    VI_masked = VI;
    VI_masked(~mask) = NaN;  % NaN for transparency in plotting (or 0 for black fill)

    % --- 5. Display the Result as an Image ---
    imagesc(xi, yi, VI_masked);
    hold on;
    
    axis xy;          % make y increase upward (important!)
    axis equal;       % keep aspect ratio 1:1
    xlim([-35, 35]);  % now works
    ylim([-35, 35]);
    colormap gray;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    %title('Masked Around Gaze Angle');
    axis off;            % removes axes, ticks, label


    exportgraphics(gcf, '/Users/zacharykelly/Desktop/virtually_foveated.png', 'Resolution', 300);
    close all; 

    virtually_foveated_frame = imread('/Users/zacharykelly/Desktop/virtually_foveated.png');
 

    return ; 
end

% local function to adjust plot ...
function plotCircle3d(center,normal,radius)
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,3))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'g-');
end