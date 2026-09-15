% Define the solid angle subtended by each world-camera pixel.
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
%
% The same array is saved to derived/deltaSteradians.mat.

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

% Match world_frame_visual_angle_to_steradians in world_util.py exactly.
% The final visual-angle dimension is [azimuth, elevation], in degrees.
azimuth = deg2rad(visualAngles(:, :, 1));
elevation = deg2rad(visualAngles(:, :, 2));
cosElevation = cos(elevation);
unitDirections = cat(3, ...
    cosElevation .* sin(azimuth), ...
    -sin(elevation), ...
    cosElevation .* cos(azimuth));

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
    'defineDeltaSteradians:UnexpectedOutputSize', ...
    'The solid-angle map does not match the calibrated image size.');
assert(all(isfinite(deltaSteradians), 'all') && ...
    all(deltaSteradians > 0, 'all'), ...
    'defineDeltaSteradians:InvalidOutput', ...
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
    'defineDeltaSteradians:FieldOfViewMismatch', ...
    ['The pixel solid-angle sum differs from the independently integrated ' ...
    'field of view by %.3f%%.'], 100 * relativeFieldOfViewError);

% Show the map
figure
imagesc(deltaSteradians);
colorbar
title('Steradians per pixel')

% Save only the reusable pixel map plus human-readable provenance. The scalar
% totals are validation diagnostics and can always be recomputed from the map
% and intrinsics.
saveFileName = fullfile(projectRoot, 'derived', 'deltaSteradians.mat');
readme = sprintf([ ...
    'Created by defineDeltaSteradians.\n' ...
    'deltaSteradians -- 480-by-640 solid angle per world-camera pixel, ' ...
    'in steradians.\n' ...
    'The array is aligned with raw world-camera [row, column] coordinates.\n' ...
    'Summed pixel solid angle: %.12g sr.\n' ...
    'Independently integrated field of view: %.12g sr.\n'], ...
    summedPixelSteradians, integratedFieldOfViewSteradians);
save(saveFileName,'readme','deltaSteradians');

fprintf('Summed pixel solid angle: %2.2f sr\n', summedPixelSteradians);

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
