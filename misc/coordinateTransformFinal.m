function retinalImage = coordinateTransformFinal(I, gaze_angle, fisheyeIntrinsicsPath, transformationPath, options)
% Transform world camera image in pixels to retinal image in degrees
%
% Syntax:
%   retinalImage = coordinateTransformFinal(I, gaze_angle, fisheyeIntrinsicsPath, transformationPath)
%
% Description:
%   This function takes in an image (I) and transforms that image to a
%   retinal representation through several stages:
%     - From world camera pixels to world camera degrees of visual angle,
%       making use of the emipirically measured world camera instrinsics.
%     - From world camera degrees to eye coordinate degrees. This transform
%       makes use of a projection matrix calculated for each participant
%       based upon a gaze calibration procedure.
%     - From eue coordinate degrees to retinal degrees. This final
%       transform corrects for rotation of the eye from frame-to-frame.
%
% Inputs:
%   I                     - Matrix, double. The world camera image. 
%   gaze_angle            - 2x1 vector specificying the azimuth and
%                           elevation rotation of the eye in degrees.
%                           Coordinate frame is that defined by the
%                           gkaModelEye.
%   fisheyeIntrinsicsPath - String or char vector full path to previously 
%                           measured world camera intrinsics.
%   transformationPath    - String or char vector full path to previously
%                           calculated camera projection object.
%
% Optional key/value pairs:
%  "FOVradius"            - Scalar. The field of view of the retinal image
%                           in degrees of visual angle.
%  "degPerSample"         - Scalar. The visual angle degrees per sample in
%                           the returned retinal image.
%  "forceRecalc"          - Logical. Forces reloading and recalculation
%                           based upon the camera intrinsics and the
%                           projection from world to eye coordinates.
%
% Outputs:
%   retinalImage          - Matrix. The retinal image in degrees of
%                           visual angle.
%
% Examples:
%{
	foo = 1;
    bar = myFunc(foo);
	fprintf('Bar = %d \n',bar);   
%}

arguments
    I double {mustBeMatrix}
    gaze_angle double {mustBeVector}
    fisheyeIntrinsicsPath
    transformationPath
    options.FOVradius double = 30
    options.degPerSample double = 0.25
    options.forceRecalc logical = true
end

% Get the camera visual field positions corresponding to positions of all
% locations on the camera sensor
[nRows, nCols] = size(I);

% Load the camera intrinsics and derive the conversion of world camera
% pixels to world camera angle degrees. Then, load the projection matrix
% from world to eye coordinates and calculate the conversion for each
% sample. We use persistent variables to avoid re-calculating all this
% across frames. There is a mechnanism as well to detect a change in the
% path to the file, or the size of the image, which triggers a
% re-calculation.
persistent fisheyeIntrinsics lastFisheyeIntrinsicsPath
persistent worldCoordsInDegrees nRowsLast nColsLast
persistent transformation lastTransformationPath eyeRotationCoordinates

% Check for conditions that would trigger a re-calculation of the angles
% from intrinsics
recalcFlag = false;
if isempty(fisheyeIntrinsics) || ...
        ~strcmp(lastFisheyeIntrinsicsPath,fisheyeIntrinsicsPath) || ...
        ~isequal(nRowsLast,nRows) || ~isequal(nColsLast,nCols) || ...
        ~strcmp(lastTransformationPath,transformationPath) || ...
        options.forceRecalc
    recalcFlag = true;
end

% We are going to recalculate the angles from intrinsics
if recalcFlag
    % Update the persistent variables
    lastFisheyeIntrinsicsPath = fisheyeIntrinsicsPath;
    nRowsLast = nRows; nColsLast = nCols;

    % Load the camera intrinsics
    fisheyeIntrinsics = load(fisheyeIntrinsicsPath).camera_intrinsics_calibration.results.Intrinsics;

    % Calculate the world pixel --> degrees conversion
    [xg, yg]          = meshgrid(1:nCols, 1:nRows);
    sensorPoints      = [xg(:),yg(:)];
    worldCoordsInDegrees = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);

    % Load the projection matrix
    lastTransformationPath = transformationPath;
    transformation = load(transformationPath).perspective_transform.fit.geometric_transform;

    % Calculate world degrees --> eye degrees
    eyeRotationCoordinates = transformPointsForward(transformation, worldCoordsInDegrees);

end

% Extract x, y, v from the eyeRotationCoordinates
x = eyeRotationCoordinates(:, 1);
y = eyeRotationCoordinates(:, 2);
v = I(:);

% Define the Regular Grid for Interpolation
xmin = gaze_angle(1)-options.FOVradius;
xmax = gaze_angle(1)+options.FOVradius;
ymin = gaze_angle(2)-options.FOVradius;
ymax = gaze_angle(2)+options.FOVradius;

num_x_points = (xmax-xmin) / options.degPerSample;
num_y_points = (ymax-ymin) / options.degPerSample;

xi = linspace(xmin, xmax, num_x_points);
yi = linspace(ymin, ymax, num_y_points);
[XX, YY] = meshgrid(xi, yi);

% Interpolate the irregularly scattered samples onto a grid
retinalImage = griddata(x, y, v, XX, YY, 'linear');

end
