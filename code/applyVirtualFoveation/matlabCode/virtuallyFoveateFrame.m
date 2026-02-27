function retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
% Transform world camera image in pixels to retinal image in degrees
%
% Syntax:
%   retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
%
% Description:
%   This function takes in an image (I) and transforms that image to a
%   retinal representation by mapping each world-camera pixel to a (x,y)
%   location in degrees of visual angle (via empirically measured fisheye
%   intrinsics), then resampling onto a regular grid centered at gaze_angle.
%
%   NOTE: This version supports a FIXED output size via options.desiredN.
%   If desiredN is provided, the output is always desiredN-by-desiredN
%   (or desiredN-by-desiredN-by-3 for RGB), eliminating 480/481 jitter.
%
% Inputs:
%   I                     - Matrix, double. The world camera image.
%                           Can be grayscale (nRows x nCols) or RGB
%                           (nRows x nCols x 3).
%   gaze_angle            - 2x1 vector specifying azimuth and elevation
%                           (degrees) in the SAME coordinate frame as the
%                           output of anglesFromIntrinsics().
%   fisheyeIntrinsicsPath - String/char full path to intrinsics .mat file.
%
% Optional key/value pairs (name-value):
%   "FOVradius"            - Scalar. Half-width field of view radius (deg).
%   "degPerSample"         - Scalar. (Only used if desiredN is empty.)
%   "desiredN"             - Integer. Fixed output width/height in samples.
%   "forceRecalc"          - Logical. Force recomputing intrinsics mapping.
%
% Output:
%   retinalImage           - Resampled image on regular degree grid.
%
% Example:
%   R = virtuallyFoveateFrame(I, [0;0], intrPath, "desiredN", 480);

    arguments
        I double
        gaze_angle double {mustBeVector}
        fisheyeIntrinsicsPath
        options.FOVradius (1,1) double = 60
        options.degPerSample (1,1) double = 0.25
        options.desiredN (1,1) double {mustBeInteger, mustBePositive} = 480
        options.forceRecalc (1,1) logical = false
    end

    % Basic image handling
    [nRows, nCols, nChannels] = size(I);

    % Keep your original orientation transform
    I = fliplr(imrotate(I, 180));

    % Persistent caches to avoid recomputation across frames
    persistent fisheyeIntrinsics lastFisheyeIntrinsicsPath
    persistent worldCoordsInDegrees nRowsLast nColsLast

    % Recalc conditions
    recalcFlag = false;
    no_fisheye_intrinsics = isempty(fisheyeIntrinsics);
    new_fisheye_intrinsics_path = isempty(lastFisheyeIntrinsicsPath) || ~strcmp(lastFisheyeIntrinsicsPath, fisheyeIntrinsicsPath);
    new_image_shape = isempty(nRowsLast) || isempty(nColsLast) || ~isequal(nRowsLast, nRows) || ~isequal(nColsLast, nCols);

    if any([no_fisheye_intrinsics, new_fisheye_intrinsics_path, new_image_shape]) || options.forceRecalc
        recalcFlag = true;
    end

    if recalcFlag
        lastFisheyeIntrinsicsPath = fisheyeIntrinsicsPath;
        nRowsLast = nRows;
        nColsLast = nCols;

        % Load intrinsics (adjust this field path if your .mat differs)
        tmp = load(fisheyeIntrinsicsPath);
        fisheyeIntrinsics = tmp.camera_intrinsics_calibration.results.Intrinsics;

        % Compute mapping from pixel -> (deg_x, deg_y)
        [xg, yg] = meshgrid(1:nCols, 1:nRows);
        sensorPoints = [xg(:), yg(:)];

        % anglesFromIntrinsics must return Nx2 degrees
        worldCoordsInDegrees = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    end

    % Extract scattered (x,y) degree coords for each pixel
    x = worldCoordsInDegrees(:, 1);
    y = worldCoordsInDegrees(:, 2);

    % ----- FIXED OUTPUT GRID SIZE -----
    N = options.desiredN;  % the single source of truth for output size

    xi = gaze_angle(1) + linspace(-options.FOVradius, options.FOVradius, N);
    yi = gaze_angle(2) + linspace(-options.FOVradius, options.FOVradius, N);
    [XX, YY] = meshgrid(xi, yi);

    % Resample onto the grid (griddata output matches size(XX))
    if nChannels <= 1
        retinalImage = griddata(x, y, I(:), XX, YY, 'linear');
    else
        retinalImage = zeros(N, N, nChannels, 'like', I);
        for c = 1:nChannels
            v = I(:, :, c);
            retinalImage(:, :, c) = griddata(x, y, v(:), XX, YY, 'linear');
        end
    end
end