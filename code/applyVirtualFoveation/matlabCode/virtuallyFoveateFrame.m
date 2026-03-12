function retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
% Resample a world-camera image onto a gaze-centered retinal coordinate grid
%
% Syntax:
%   retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath)
%   retinalImage = virtuallyFoveateFrame(I, gaze_angle, fisheyeIntrinsicsPath, options)
%
% Description:
%   This function transforms a world-camera image from pixel coordinates
%   into a retinal image represented on a regular grid in degrees of
%   visual angle. It uses empirically measured fisheye camera intrinsics
%   to map each input pixel to an angular location, then resamples the
%   image onto a gaze-centered grid spanning a specified field of view.
%
%   The function supports both grayscale and RGB inputs. To avoid repeated
%   recomputation across frames, the fisheye intrinsics and the pixel-to-
%   degrees mapping are cached using persistent variables. The output size
%   can be fixed by specifying options.desiredN, which ensures a constant
%   output resolution across frames.
%
% Inputs:
%   I                     - Numeric matrix. Input world-camera image. May
%                           be grayscale (nRows x nCols) or RGB
%                           (nRows x nCols x 3).
%   gaze_angle            - Numeric vector. Two-element vector specifying
%                           azimuth and elevation in degrees, in the same
%                           coordinate frame returned by
%                           anglesFromIntrinsics().
%   fisheyeIntrinsicsPath - Char/string. Full path to a .mat file
%                           containing the fisheye camera intrinsics used
%                           for the pixel-to-angle mapping.
%
% Optional key/value pairs:
%   FOVradius             - Scalar. Half-width of the retinal field of
%                           view in degrees around gaze center.
%   degPerSample          - Scalar. Angular spacing per output sample in
%                           degrees. Included for compatibility with
%                           versions that determine output size from
%                           sampling density.
%   desiredN              - Positive integer. Fixed output width and height
%                           in samples. The returned image will be
%                           desiredN-by-desiredN for grayscale input or
%                           desiredN-by-desiredN-by-3 for RGB input.
%   forceRecalc           - Logical. If true, force recomputation of the
%                           fisheye intrinsics mapping even if cached
%                           values are available.
%
% Outputs:
%   retinalImage          - Numeric matrix. Image resampled onto a regular
%                           retinal grid centered at gaze_angle.
%
% Examples:
%{
    intrPath = "/path/to/intrinsics_calibration.mat";
    retinalImage = virtuallyFoveateFrame(I, [0; 0], intrPath);

    retinalImage = virtuallyFoveateFrame( ...
        I, ...
        [5; -2], ...
        intrPath, ...
        "FOVradius", 60, ...
        "desiredN", 480 ...
    );
%}R = virtuallyFoveateFrame(I, [0;0], intrPath, "desiredN", 480);

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