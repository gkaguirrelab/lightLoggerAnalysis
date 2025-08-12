function coordinateTransform(imgPts, worldPts, fisheyeIntrinsics, myChoice)

% Need to make available to this function the fisheye intrinsics, and the
% set of "imgPts" and "worldPts".

% housekeeping
plotChoices = {'worldImage','slopeMap','interceptMap'};
myChoice = 'worldImage';

% Obtain the set of gaze calibration targets as seen by the world camera,
% expressed in sensor coordinate locations
gazeTargetSensorCoord = imgPts;

% Obtain the set of eye rotations that correspond to these gaze target
% locations
veridicalEyeRotations = worldPts;

% Transform these to visual field locations as seen by the world camera,
% and then calculate the affine transform between these locations and the
% eye angles of rotation
gazeTargetCameraFieldCoord = anglesFromIntrinsics(gazeTargetSensorCoord, fisheyeIntrinsics);
tform = fitgeotform2d(gazeTargetCameraFieldCoord,veridicalEyeRotations,'projective');

% Transform the gaze targets as seen by the camera into the eye rotation
% coordinate space
gazeTargetEyeRotation = transformPointsForward(tform,gazeTargetCameraFieldCoord);

% Obtain the data map in the sensor (u,v) space.
switch myChoice
    case 'worldImage'
        I = imread('/Users/zacharykelly/Downloads/exampleWorldImage.png');
        I = imresize(I,[480,640]);
        I = mean(I,3);
        myMap = 'gray';
        barRange = [0,255];
        gazePlotFlag = true;
    case 'interceptMap'
        I = hFigHighIntercept;
        myMap = 'hot';
        barRange = [-3.5,-2];
        gazePlotFlag = false;
    case 'slopeMap'
        I = hFigHighSlope;
        myMap = 'hot';
        barRange = [-2.5,-2];
        gazePlotFlag = false;
end

% Get the camera visual field positions corresponding to positions of all
% locations on the camera sensor
nCols = 640; nRows = 480;
[xg, yg] = meshgrid(1:nCols, 1:nRows);
sensorPoints = [xg(:),yg(:)];
visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);

% Show what I looks like in the sensor grid coordinate space, and add the
% gaze targets
figure
surf(reshape(sensorPoints(:,1),480,640),reshape(sensorPoints(:,2),480,640),I,'edgeColor','none');
view([0,-90]); % This command rotates the plot so we are looking straight at it
colormap(myMap)
hold on
if gazePlotFlag
    plot3(gazeTargetSensorCoord(:,1),gazeTargetSensorCoord(:,2),repmat(0,size(gazeTargetSensorCoord(:,1))),'xr');
end
title('Camera image in sensor coords')
xlabel('sensor position [pixels]');
ylabel('sensor position [pixels]');
colorbar
clim(barRange);

% Now show what I looks like in the camera visual field coordinate space
figure
surf(reshape(visualFieldPoints(:,1),480,640),reshape(visualFieldPoints(:,2),480,640),I,'edgeColor','none');
view([90,90]); % This command rotates the plot so we are looking straight at it
colormap(myMap)
hold on
if gazePlotFlag
    plot3(gazeTargetCameraFieldCoord(:,1),gazeTargetCameraFieldCoord(:,2),repmat(255,size(gazeTargetCameraFieldCoord(:,1))),'xr');
end
title('Camera image in camera visual angle coords')
colorbar
clim(barRange);

% Transform the camera visual field points to eye rotation coordinate space
% using the previously calculated tform
eyeRotationCoordinates = transformPointsForward(tform,visualFieldPoints);

% Now show what I looks like in eye rotation coordinates
figure
surf(reshape(eyeRotationCoordinates(:,1),480,640),reshape(eyeRotationCoordinates(:,2),480,640),I,'edgeColor','none');
view([0,-90]);
colormap(myMap)
hold on
if gazePlotFlag
    plot3(gazeTargetEyeRotation(:,1),gazeTargetEyeRotation(:,2),repmat(0,size(gazeTargetEyeRotation(:,1))),'xr');
    plot3(veridicalEyeRotations(:,1),veridicalEyeRotations(:,2),repmat(0,size(veridicalEyeRotations(:,1))),'xb');
end
title('camera image in eye rotation coords')
colorbar
clim(barRange);

% Plot this for 60 degree of eccentricity
idx = vecnorm(eyeRotationCoordinates,2,2) > 60;
figure
subI = I; subI(idx)=nan;
surf(reshape(eyeRotationCoordinates(:,1),480,640),reshape(eyeRotationCoordinates(:,2),480,640),subI,'edgeColor','none');
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
clim(barRange);



function plotCircle3d(center,normal,radius)
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,3))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'g-');
end