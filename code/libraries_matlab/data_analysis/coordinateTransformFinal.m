function [SensorFigure, CameraFigure, EyeFigure] = coordinateTransformFinal(fisheyeIntrinsics, maps, calibFile, myChoice, lightLevel)

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
    
    arguments
        fisheyeIntrinsics
        maps 
        calibFile
        myChoice {mustBeMember(myChoice, {'worldImage', 'interceptMap', 'slopeMap'})} = 'interceptMap'
        lightLevel {mustBeMember(lightLevel, {'highAS', 'lowAS', 'allAS'})} = 'allAS'
    end

    % Obtain the set of gaze calibration targets as seen by the world camera,
    % expressed in sensor coordinate locations
    % Obtain the set of eye rotations that correspond to these gaze target
    % locations
    [gazeTargetSensorCoord, veridicalEyeRotations] = get_calibration_dots(calibFile);
    
    %% AFFINE TRANSFORM
    
    % Transform these points to visual field locations as seen by the world camera,
    % and then calculate the affine transform between these locations and the
    % eye angles of rotation
    gazeTargetCameraFieldCoord = anglesFromIntrinsics( gazeTargetSensorCoord, fisheyeIntrinsics );
    tform = fitgeotform2d( gazeTargetCameraFieldCoord, veridicalEyeRotations, 'projective' );
    
    % Transform the gaze targets as seen by the camera into the eye rotation
    % coordinate space
    gazeTargetEyeRotation = transformPointsForward( tform, gazeTargetCameraFieldCoord );
    
    % Obtain the data map in the sensor (u,v) space.
    switch myChoice
        case 'worldImage'
            I = imread('/Users/zacharykelly/Downloads/exampleWorldImage.png');
            I = imresize(I,[480,640]);
            I = mean(I,3);
            myMap = 'gray'; barRange = [0,255]; gazePlotFlag = false;
        case 'interceptMap'
            if ~isfield(maps,'intercept') || ~isfield(maps.intercept,lightLevel)
                error('maps.intercept.%s is missing or empty.', lightLevel);
            end
            I = maps.intercept.(lightLevel);             % numeric matrix
            myMap = 'hot';  barRange = [-3.5,-2]; gazePlotFlag = false;

        case 'slopeMap'
            if ~isfield(maps,'slope') || ~isfield(maps.slope,lightLevel)
                error('maps.slope.%s is missing or empty.', lightLevel);
            end
            I = maps.slope.(lightLevel);                 % numeric matrix
            myMap = 'hot';  barRange = [-2.5,-2]; gazePlotFlag = false;
    end
    
    % Get the camera visual field positions corresponding to positions of all
    % locations on the camera sensor
    [nRows, nCols]    = size(I);
    [xg, yg]          = meshgrid(1:nCols, 1:nRows);
    sensorPoints      = [xg(:),yg(:)];
    visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);

    % Show what I looks like in the sensor grid coordinate space, and add the
    % gaze targets
    SensorFigure = figure;
    surf(reshape(sensorPoints(:,1),480,640), reshape(sensorPoints(:,2),nRows,nCols), I, 'edgeColor','none');
    view([0,-90]); % This command rotates the plot so we are looking straight at it
    colormap(myMap)
    hold on
    if gazePlotFlag
        plot3(gazeTargetSensorCoord(:,1),gazeTargetSensorCoord(:,2),repmat(0,size(gazeTargetSensorCoord(:,1))),'xr');
    end
    title('Camera Image in Sensor Coordinates')
    xlabel('sensor position [pixels]');
    ylabel('sensor position [pixels]');
    colorbar
    %clim(barRange);
    
    % Now show what I looks like in the camera visual field coordinate space
    CameraFigure = figure;
    surf(reshape(visualFieldPoints(:,1),480,640),reshape(visualFieldPoints(:,2),nRows,nCols),I,'edgeColor','none');
    view([90,90]); % This command rotates the plot so we are looking straight at it
    colormap(myMap)
    hold on
    if gazePlotFlag
        plot3(gazeTargetCameraFieldCoord(:,1),gazeTargetCameraFieldCoord(:,2),repmat(255,size(gazeTargetCameraFieldCoord(:,1))),'xr');
    end
    title('Camera Image in Camera Visual Angle Coordinates')
    colorbar
    %clim(barRange);
    
    % Transform the camera visual field points to eye rotation coordinate space
    % using the previously calculated tform
    eyeRotationCoordinates = transformPointsForward(tform,visualFieldPoints);
    
    % Now show what I looks like in eye rotation coordinates
    EyeFigure = figure;
    surf(reshape(eyeRotationCoordinates(:,1),480,640),reshape(eyeRotationCoordinates(:,2),nRows,nCols),I,'edgeColor','none');
    view([0,-90]);
    colormap(myMap)
    hold on
    if gazePlotFlag
        plot3(gazeTargetEyeRotation(:,1),gazeTargetEyeRotation(:,2),repmat(0,size(gazeTargetEyeRotation(:,1))),'xr');
        plot3(veridicalEyeRotations(:,1),veridicalEyeRotations(:,2),repmat(0,size(veridicalEyeRotations(:,1))),'xb');
    end
    title('camera image in eye rotation coords')
    colorbar
    %clim(barRange);
    
    % Plot this for 60 degree of eccentricity
    idx = vecnorm(eyeRotationCoordinates,2,2) > 60;
    figure
    subI = I; subI(idx)=nan;
    surf(reshape(eyeRotationCoordinates(:,1),480,640),reshape(eyeRotationCoordinates(:,2),nRows,nCols),subI,'edgeColor','none');
    view([0,-90]);
    colormap(myMap)
    hold on
    if gazePlotFlag
        plot3(gazeTargetEyeRotation(:,1),gazeTargetEyeRotation(:,2),repmat(-5,size(gazeTargetEyeRotation(:,1))),'xr');
        plot3(veridicalEyeRotations(:,1),veridicalEyeRotations(:,2),repmat(-5,size(veridicalEyeRotations(:,1))),'xb');
    end
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

end

% local function to obtain gaze calibration dots for affine transform
function [imgPts, worldPts] = get_calibration_dots(calibFile)
    arguments
        calibFile {mustBeText} = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/SophiaGazeCalib2/W_mjpeg.avi';
    end
    
    vr = VideoReader(calibFile);
    
    % video parameters
    fps = vr.FrameRate;
    startSec = 82;
    endSec   = 276;
    
    % move to start time
    vr.CurrentTime = startSec;
    
    % prepare max brightness mask
    maxBright = [];
    while hasFrame(vr) && vr.CurrentTime <= endSec
        frm = readFrame(vr);
        g   = rgb2gray(frm);
        % only keep pixels brighter than a threshold
        % normalize and threshold
        g = double(g);
        g(g < 200) = 0;
        % combine max
        if isempty(maxBright)
            maxBright = g;
        else
            maxBright = max(maxBright, g);
        end
    end
    
    % convert to uint8 for display
    maxBright = uint8(maxBright);
    
    % extract dot centroids
    BW = imbinarize(maxBright, graythresh(maxBright));  
    stats = regionprops(BW, 'Centroid');
    imgPts = vertcat(stats.Centroid); % k×2 [x y] in pixel coords
    x = imgPts(:,1);
    y = imgPts(:,2);
    r = sqrt(x.^2 + y.^2);
    
    % world point coords
    %% 
    worldPts = [ ...
            -20, -20;   -20, 0;   -20, 20;   -15, -15;  -15, 15; ...
            -10, -10;   -10, 0;   -10, 10;   -5, -5;    -5, 5;   ...
             0, -20;     0, -10;   0, 0;      0, 10;     0, 20;  ...
             5, -5;      5, 5;     10, -10;   10, 0;     10, 10; ...
             15, -15;    15, 15;   20, -20;   20, 0;     20, 20];
    end

% local function to obtain the theta and phi values from the camera
% intrinsics
function visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics)
    cx = fisheyeIntrinsics.DistortionCenter(1);
    cy = fisheyeIntrinsics.DistortionCenter(2);
    xC = sensorPoints(:,1) - cx; yC = sensorPoints(:,2) - cy;
    
    r = sqrt(xC.^2 + yC.^2);
    coeffs = fisheyeIntrinsics.MappingCoefficients;
    a0=coeffs(1); a2=coeffs(2); a3=coeffs(3); a4=coeffs(4);
    
    theta = nan(size(r));
    for k=1:numel(r)
        func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - r(k);
        theta(k) = fzero(func,[0,pi]);
    end
    
    phi = atan2(yC,xC); R = 1;
    X = R*sin(theta).*cos(phi);
    Y = R*sin(theta).*sin(phi);
    Z = R*cos(theta);
    
    azi = rad2deg(atan(Y./Z));
    ele = rad2deg(atan(X./Z));
    
    visualFieldPoints = [azi,ele];
end

% local function to adjust plot ...
function plotCircle3d(center,normal,radius)
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,3))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'g-');
end