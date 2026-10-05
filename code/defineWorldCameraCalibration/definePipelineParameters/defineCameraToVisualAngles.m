% Define the conversions of pixel locations to visual angles.
%
% The calibrated fisheye model maps every 640-by-480 pixel center to an
% azimuth and elevation. Following world_util.py, those angles are converted
% to unit viewing directions. The magnitude of the cross product between the
% row and column direction derivatives is the local spherical-area Jacobian,
% expressed as steradians per pixel.
%
% Output:
%   deltaSteradians - 480-by-640 double array. Element (row, column) is the
%                     solid angle subtended by that world-camera pixel.
%   eccentricityMap - 480-by-640 double array. Element (row, column) is the
%                     visual field eccentricity in degrees from the optical center.
%   unitDirections  - 480-by-640-by-3 double array. Element (row, column, :)
%                     is the 3D Cartesian unit vector for that pixel's viewing direction.
%
% The arrays are saved to derived/cameraToVisualAngles.mat.

projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');

% Load the calibrated fisheye model using the canonical, correctly spelled
% filename and variable name.
intrinsicsPath = fullfile( ...
    projectRoot, 'derived', 'arducamB0392cameraIntrinsics.mat');
assert(isfile(intrinsicsPath), ...
    'defineDeltaSteradians:MissingIntrinsics', ...
    'Could not find the derived world-camera fisheye intrinsics file.');

intrinsicsData = load(intrinsicsPath, 'arducamB0392cameraIntrinsics');
assert(isfield(intrinsicsData, 'arducamB0392cameraIntrinsics'), ...
    'defineDeltaSteradians:MissingIntrinsicsVariable', ...
    ['The intrinsics MAT file must contain ' ...
    'arducamB0392cameraIntrinsics.']);
fisheyeIntrinsics = ...
    intrinsicsData.arducamB0392cameraIntrinsics.results.Intrinsics;

% The calibration and raw recordings both use 480 rows by 640 columns.
imageSize = double(fisheyeIntrinsics.ImageSize(:).');
assert(isequal(imageSize, [480 640]), ...
    'defineDeltaSteradians:UnexpectedImageSize', ...
    'Expected 480-by-640 intrinsics, but found %g-by-%g.', ...
    imageSize(1), imageSize(2));
rows = imageSize(1);
columns = imageSize(2);

% anglesFromIntrinsics expects MATLAB one-based [x, y] pixel-center
% coordinates. Flatten in MATLAB column-major order and reshape the returned
% angles in the same order so each result remains aligned to its source pixel.
[xCoordinates, yCoordinates] = meshgrid(1:columns, 1:rows);
sensorPoints = [xCoordinates(:), yCoordinates(:)];
visualAngles = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
visualAngles = reshape(visualAngles, rows, columns, 2);

% Extract native azimuth and elevation maps directly in degrees
nativeAzimuthMap = visualAngles(:, :, 1);
nativeElevationMap = visualAngles(:, :, 2);

% Match world_frame_visual_angle_to_steradians in world_util.py exactly.
azimuth = deg2rad(nativeAzimuthMap);
elevation = deg2rad(nativeElevationMap);
cosElevation = cos(elevation);
unitDirections = cat(3, ...
    cosElevation .* sin(azimuth), ...
    -sin(elevation), ...
    cosElevation .* cos(azimuth));

% Calculate eccentricity (in degrees) from the forward optical axis. 
% The forward axis is [0, 0, 1], so the dot product with unit directions
% is exactly the Z-component.
eccentricityMap = rad2deg(acos(unitDirections(:, :, 3)));

% NumPy's gradient(..., edge_order=2) uses centered differences internally
% and second-order one-sided differences at both image boundaries. The local
% helper below reproduces those formulas explicitly for MATLAB arrays.
directionChangePerRow = secondOrderFiniteDifference(unitDirections, 1);
directionChangePerColumn = secondOrderFiniteDifference(unitDirections, 2);

% The cross-product magnitude is the area of the local parallelogram on the
% unit sphere, which is the solid angle represented by the pixel.
areaVectors = cross( ...
    directionChangePerColumn, directionChangePerRow, 3);
deltaSteradians = sqrt(sum(areaVectors.^2, 3));

assert(isequal(size(deltaSteradians), [rows columns]), ...
    'defineCameraToVisualAngles:UnexpectedOutputSize', ...
    'The solid-angle map does not match the calibrated image size.');
assert(all(isfinite(deltaSteradians), 'all') && ...
    all(deltaSteradians > 0, 'all'), ...
    'defineCameraToVisualAngles:InvalidOutput', ...
    'Every pixel solid angle must be finite and positive.');

% Compare the summed per-pixel approximation with an independent integration
% over the calibrated rectangular sensor boundary. The finite-difference map
% samples area at pixel centers, so close agreement rather than exact equality
% is expected at the outer half-pixel boundary.
summedPixelSteradians = sum(deltaSteradians, 'all');
integratedFieldOfViewSteradians = ...
    calculateFisheyeSolidAngle(fisheyeIntrinsics);
relativeFieldOfViewError = abs( ...
    summedPixelSteradians - integratedFieldOfViewSteradians) ./ ...
    integratedFieldOfViewSteradians;
assert(relativeFieldOfViewError < 0.01, ...
    'defineCameraToVisualAngles:FieldOfViewMismatch', ...
    ['The pixel solid-angle sum differs from the independently integrated ' ...
    'field of view by %.3f%%.'], 100 * relativeFieldOfViewError);

% Show the solid angle map
figure
imagesc(deltaSteradians);
colorbar
axis equal
title('Steradians per pixel')

% Show the eccentricity map
figure
imagesc(eccentricityMap);
axis equal
colorbar
title('Eccentricity (degrees)')

% Save the reusable maps plus human-readable provenance.
saveFileName = fullfile(projectRoot, 'derived', 'cameraToVisualAngles.mat');
readme = sprintf([ ...
    'Created by defineCameraToVisualAngles.\n' ...
    'deltaSteradians -- 480-by-640 solid angle per world-camera pixel, ' ...
    'in steradians.\n' ...
    'eccentricityMap -- 480-by-640 visual field eccentricity from optical center, ' ...
    'in degrees.\n' ...
    'unitDirections  -- 480-by-640-by-3 unit viewing vectors.\n' ...
    'nativeAzimuthMap -- 480-by-640 azimuth angles, in degrees.\n' ...
    'nativeElevationMap -- 480-by-640 elevation angles, in degrees.\n' ...
    'The arrays are aligned with raw world-camera [row, column] coordinates.\n' ...
    'Summed pixel solid angle: %.12g sr.\n' ...
    'Independently integrated field of view: %.12g sr.\n'], ...
    summedPixelSteradians, integratedFieldOfViewSteradians);

% Added nativeAzimuthMap and nativeElevationMap to the saved variables
save(saveFileName,'readme','deltaSteradians','eccentricityMap','unitDirections','nativeAzimuthMap', 'nativeElevationMap');



% LOCAL FUNCTIONS

function derivative = secondOrderFiniteDifference(values, dimension)
% Reproduce numpy.gradient(..., edge_order=2) at unit sample spacing.

dimensionLength = size(values, dimension);
assert(dimensionLength >= 3, ...
    'defineDeltaSteradians:FiniteDifferenceSize', ...
    'Second-order finite differences require at least three samples.');

derivative = zeros(size(values), 'like', values);
allIndices = repmat({':'}, 1, ndims(values));

% Interior samples use the centered difference (next - previous) / 2.
target = allIndices;
previous = allIndices;
next = allIndices;
target{dimension} = 2:(dimensionLength - 1);
previous{dimension} = 1:(dimensionLength - 2);
next{dimension} = 3:dimensionLength;
derivative(target{:}) = ...
    (values(next{:}) - values(previous{:})) ./ 2;

% Boundary samples use NumPy's second-order one-sided coefficients.
first = allIndices;
second = allIndices;
third = allIndices;
first{dimension} = 1;
second{dimension} = 2;
third{dimension} = 3;
derivative(first{:}) = ...
    (-3 .* values(first{:}) + 4 .* values(second{:}) - ...
    values(third{:})) ./ 2;

last = allIndices;
penultimate = allIndices;
antepenultimate = allIndices;
last{dimension} = dimensionLength;
penultimate{dimension} = dimensionLength - 1;
antepenultimate{dimension} = dimensionLength - 2;
derivative(last{:}) = ...
    (3 .* values(last{:}) - 4 .* values(penultimate{:}) + ...
    values(antepenultimate{:})) ./ 2;

end