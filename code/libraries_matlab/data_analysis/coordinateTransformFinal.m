function virtually_foveated_frame = coordinateTransformFinal(I, fisheyeIntrinsics, transformation, center_offset)
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
    

    % Obtain the set of gaze calibration targets as seen by the world camera,
    % expressed in sensor coordinate locations
    % Obtain the set of eye rotations that correspond to these gaze target
    % locations


    
    % Transform the gaze targets as seen by the camera into the eye rotation
    % coordinate space
    %gazeTargetEyeRotation = transformPointsForward( tform, gazeTargetCameraFieldCoord );
    
    % Obtain the data map in the sensor (u,v) space.
    I = mean(I,3);
    myMap = 'gray'; 
    barRange = [0,255]; 
    gazePlotFlag = false;

    
    % Get the camera visual field positions corresponding to positions of all
    % locations on the camera sensor
    [nRows, nCols]    = size(I);
    [xg, yg]          = meshgrid(1:nCols, 1:nRows);
    sensorPoints      = [xg(:),yg(:)];
    visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics, center_offset);

    % Show what I looks like in the sensor grid coordinate space, and add the
    % gaze targets
    %{
    SensorFigure = figure;
    surf(reshape(sensorPoints(:,1),480,640), reshape(sensorPoints(:,2),nRows,nCols), I, 'edgeColor','none');
    view([0,-90]); % This command rotates the plot so we are looking straight at it
    colormap(myMap)
    hold on
    %if gazePlotFlag
    %    plot3(gazeTargetSensorCoord(:,1),gazeTargetSensorCoord(:,2),repmat(0,size(gazeTargetSensorCoord(:,1))),'xr');
    %end
    title('Camera Image in Sensor Coordinates')
    xlabel('sensor position [pixels]');
    ylabel('sensor position [pixels]');
    %colorbar
    %clim(barRange);
    %}

    % Now show what I looks like in the camera visual field coordinate space
    %{
    CameraFigure = figure;
    surf(reshape(visualFieldPoints(:,1),nRows,nCols),reshape(visualFieldPoints(:,2),nRows,nCols),I,'edgeColor','none');
    view([90,90]); % This command rotates the plot so we are looking straight at it
    colormap(myMap)
    hold on
    %if gazePlotFlag
    %    plot3(gazeTargetCameraFieldCoord(:,1),gazeTargetCameraFieldCoord(:,2),repmat(255,size(gazeTargetCameraFieldCoord(:,1))),'xr');
    %end
    title('Camera Image in Camera Visual Angle Coordinates')
    %colorbar
    %clim(barRange);
    %}

    % Transform the camera visual field points to eye rotation coordinate space
    % using the previously calculated tform
    eyeRotationCoordinates = transformPointsForward(transformation, visualFieldPoints);
    
    % Now show what I looks like in eye rotation coordinates
    %{
    EyeFigure = figure;
    surf(reshape(eyeRotationCoordinates(:,1),480,640),reshape(eyeRotationCoordinates(:,2),nRows,nCols),I,'edgeColor','none');
    view([0,-90]);
    axis ij;    
    %colormap(myMap)
    hold on
    %if gazePlotFlag
    %    plot3(gazeTargetEyeRotation(:,1),gazeTargetEyeRotation(:,2),repmat(0,size(gazeTargetEyeRotation(:,1))),'xr');
    %    plot3(veridicalEyeRotations(:,1),veridicalEyeRotations(:,2),repmat(0,size(veridicalEyeRotations(:,1))),'xb');
    %end
    title('camera image in eye rotation coords')
    %colorbar
    %clim(barRange);
    %}

    % Plot this for 60 degree of eccentricity
    idx = vecnorm(eyeRotationCoordinates,2,2) > 60;
    subI = I; subI(idx)=nan;

    virtually_foveated_X = reshape(eyeRotationCoordinates(:,1),nRows,nCols); 
    virtually_foveated_Y = reshape(eyeRotationCoordinates(:,2),nRows,nCols); 
    
    figure;     
    surf(virtually_foveated_X, virtually_foveated_Y, subI,'edgeColor','none');    
    view([0,-90]);
    axis ij;    % NEED TO DO THIS OTHERWISE THE THING IS ROTATED (FIGURE OUT WHY BETTER)
    %colormap(myMap)
    hold on
    %if gazePlotFlag
    %    plot3(gazeTargetEyeRotation(:,1),gazeTargetEyeRotation(:,2),repmat(-5,size(gazeTargetEyeRotation(:,1))),'xr');
    %    plot3(veridicalEyeRotations(:,1),veridicalEyeRotations(:,2),repmat(-5,size(veridicalEyeRotations(:,1))),'xb');
    %end
    axis square
    grid off
    % Add some polar angle coordinate grids
    plot3([-60,60],[0,0],[-5,-5],'-g');
    plot3([0,0],[-60,60],[-5,-5],'-g');
    for r = [15,30,60]
        plotCircle3d([0 0 -5],[0 0 1],r)
    end
    title('camera image in eye rotation coords')
    xlabel('Visual angle [deg]');
    ylabel('Visual angle [deg]');
    colorbar
    %clim(barRange);

    % 
    % Regular query grid (choose resolution = image size)
    xq = linspace(min(virtually_foveated_X(:)), max(virtually_foveated_X(:)), nCols);
    yq = linspace(min(virtually_foveated_Y(:)), max(virtually_foveated_Y(:)), nRows);
    [Xq, Yq] = meshgrid(xq, yq);

    % Interpolate onto regular grid (vectorize X, Y, Z)
    Vq = griddata(virtually_foveated_X(:), virtually_foveated_Y(:), subI(:), Xq, Yq);

    % Show rasterized result
    figure;
    imagesc(xq, yq, Vq);
    axis image;
    set(gca,'YDir','normal');
    colormap gray;
    colorbar;
    title('Rasterized camera image in eye rotation coords');
    xlabel('Visual angle [deg]');
    ylabel('Visual angle [deg]');

    %figure; 
    M = Vq; 
    M(isnan(M)) = 0;               % replace NaNs with 0 if needed
    M = mat2gray(M);               % scale to [0,1]
    
   % imshow(M)

    % Assign the transformed variable 
    virtually_foveated_frame = M; 
end

% local function to adjust plot ...
function plotCircle3d(center,normal,radius)
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,3))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'g-');
end