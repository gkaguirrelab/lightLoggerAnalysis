function retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
% Transform world camera image in pixels to retinal image in degrees
%
% Syntax:
%   retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
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
%                           Can be grayscale (nRows x nCols) or RGB
%                           (nRows x nCols x 3).
%   gaze_angle            - 2x1 vector specificying the azimuth and
%                           elevation rotation of the eye in degrees.
%                           Coordinate frame is that defined by the
%                           gkaModelEye.
%   fisheyeIntrinsicsPath - String or char vector full path to previously
%                           measured world camera intrinsics.
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
%                           visual angle. If input is RGB, output is
%                           (num_y_points x num_x_points x 3).
%
% Examples:
%{
    % TODO
%}

    arguments
        I double
        gaze_angle double {mustBeVector}
        fisheyeIntrinsicsPath
        options.FOVradius double = 60
        options.degPerSample double = 0.25
        options.forceRecalc logical = false
    end

    % Get the camera visual field positions corresponding to positions of all
    % locations on the camera sensor
    [nRows, nCols, nChannels] = size(I);
    I = fliplr(imrotate(I, 180));

    % Load the camera intrinsics and derive the conversion of world camera
    % pixels to world camera angle degrees. We use persistent variables to
    % avoid re-calculating all this across frames. There is a mechnanism as
    % well to detect a change in the path to the file, or the size of the
    % image, which triggers a re-calculation.
    persistent fisheyeIntrinsics lastFisheyeIntrinsicsPath
    persistent worldCoordsInDegrees nRowsLast nColsLast

    % Check for conditions that would trigger a re-calculation of the angles
    % from intrinsics
    recalcFlag = false;
    no_fisheye_intrinsics = isempty(fisheyeIntrinsics);
    new_fisheye_intrinsics_path = ~strcmp(lastFisheyeIntrinsicsPath, fisheyeIntrinsicsPath);
    new_image_shape = ~isequal(nRowsLast, nRows) || ~isequal(nColsLast, nCols);
    recalc_conditions = [no_fisheye_intrinsics, new_fisheye_intrinsics_path, new_image_shape];

    % Preserve your exact logic/structure: recalc if ANY conditions true.
    if ( any(recalc_conditions) || options.forceRecalc)
        disp("RECALCULATING ANGLES FROM INTRINSICS DUE TO CONDITIONS");
        disp(recalc_conditions);
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
        sensorPoints      = [xg(:), yg(:)];
        worldCoordsInDegrees = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    end

    % Extract x, y from the worldCoordsInDegrees
    x = worldCoordsInDegrees(:, 1);
    y = worldCoordsInDegrees(:, 2);

    % Define the Regular Grid for Interpolation
    xmin = gaze_angle(1) - options.FOVradius;
    xmax = gaze_angle(1) + options.FOVradius;
    ymin = gaze_angle(2) - options.FOVradius;
    ymax = gaze_angle(2) + options.FOVradius;

    % IMPORTANT: linspace requires integer number of points.
    % Keep your structure but make it robust by rounding to an integer
    % (and +1 so endpoints are included consistently).
    num_x_points = round((xmax - xmin) / options.degPerSample) + 1;
    num_y_points = round((ymax - ymin) / options.degPerSample) + 1;

    xi = linspace(xmin, xmax, num_x_points);
    yi = linspace(ymin, ymax, num_y_points);
    [XX, YY] = meshgrid(xi, yi);

    % Interpolate the irregularly scattered samples onto a grid
    % - Grayscale: identical to your original logic
    % - RGB: do the exact same thing channel-by-channel
    if nChannels <= 1
        v = I(:);
        retinalImage = griddata(x, y, v, XX, YY, 'linear');
    else
        retinalImage = zeros(num_y_points, num_x_points, nChannels, 'like', I);
        for c = 1:nChannels
            v = I(:, :, c);
            retinalImage(:, :, c) = griddata(x, y, v(:), XX, YY, 'linear');
        end
    end
end