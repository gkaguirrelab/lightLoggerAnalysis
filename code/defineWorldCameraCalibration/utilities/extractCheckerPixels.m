function [rawPixelIndices, corners] = extractCheckerPixels(I, camIntrinsicsData, inCorners, showPlot)
% EXTRACTCHECKERPIXELS Extracts raw pixel indices for a 6x4 checkerboard.
%
% Inputs:
%   I                 - The raw, optically distorted image (2D grayscale or 3D RGB).
%   camIntrinsicsData - The variable containing the fisheye intrinsics.
%   showPlot          - Optional boolean (default: true). Controls if a figure is created.
%   inCorners         - Optional 4x2 array of [x, y] coordinates for the grid
%                       corners (Order: TL, TR, BR, BL). If provided, skips
%                       manual selection.
%
% Outputs:
%   rawPixelIndices   - A 4x6 cell array where each cell contains the linear
%                       pixel indices from the central 75% of the corresponding
%                       check in the raw image I.
%   corners           - The 4x2 array of [x, y] coordinates used for the grid corners.

% --- 1. Parse Inputs ---
if nargin < 3
    inCorners = [];
end
if nargin < 4
    if isempty(inCorners)
        showPlot = true;
    else
        showPlot = false;
    end
end

if ~showPlot && isempty(inCorners)
    error('Cannot perform manual corner selection without a figure. Provide ''inCorners'' or set ''showPlot'' to true.');
end

% Safely extract the fisheyeIntrinsics object
if isa(camIntrinsicsData, 'fisheyeIntrinsics')
    camIntrinsics = camIntrinsicsData;
else
    try
        camIntrinsics = camIntrinsicsData.results.Intrinsics;
    catch
        error('Could not extract fisheyeIntrinsics. Please pass the intrinsics object directly.');
    end
end

if ~isa(camIntrinsics, 'fisheyeIntrinsics')
    error('The extracted camera intrinsics must be a fisheyeIntrinsics object.');
end

[imageRows, imageCols, ~] = size(I);

% --- 2. Undistort the Image and Get Virtual Intrinsics ---
disp('Undistorting image...');
[J, newIntrinsics] = undistortFisheyeImage(I, camIntrinsics);
K = newIntrinsics.IntrinsicMatrix;

if size(J, 3) == 3, J_disp = rgb2gray(J); else, J_disp = J; end
if size(I, 3) == 3, I_disp = rgb2gray(I); else, I_disp = I; end

% --- 3. Figure and Axes Setup ---
if showPlot
    fig = figure('Name', 'Checkerboard Extraction', 'NumberTitle', 'off');
    set(fig, 'Position', [100, 100, 1200, 500]);

    ax1 = subplot(1, 2, 1);
    imshow(J_disp, 'Parent', ax1);
    title(ax1, 'Undistorted (J)');
    hold(ax1, 'on');

    ax2 = subplot(1, 2, 2);
    imshow(I_disp, 'Parent', ax2);
    title(ax2, 'Raw Distorted (I): Re-projected Boundaries');
    hold(ax2, 'on');
end

% --- 4. Corner Selection ---
if isempty(inCorners)
    axes(ax1);
    subtitle(ax1, 'Click 4 Outer Corners of the COLORED PATCHES (Order: TL, TR, BR, BL)');

    corners = zeros(4, 2);

    % Use modern drawpoint ROI for explicit white crosshairs
    for i = 1:4
        hPt = drawpoint(ax1, 'Color', 'w');
        corners(i, :) = hPt.Position;
        hPt.InteractionsAllowed = 'none';
    end
else
    corners = inCorners;
    if showPlot
        axes(ax1);
        subtitle(ax1, 'Using provided corner coordinates');
    end
end

x = corners(:, 1);
y = corners(:, 2);

if showPlot
    plot(ax1, x, y, 'rO', 'MarkerSize', 8, 'LineWidth', 2);
    plot(ax1, [x; x(1)], [y; y(1)], 'r-', 'LineWidth', 1.5);
    drawnow;
end

% --- 5. Fit the Grid Geometry ---
cols = 6;
rows = 4;

% Assume the colored patches make up 95% of the grid spacing
patch_width = 0.95;

fixed_max_x = (cols - 1) + patch_width;
fixed_max_y = (rows - 1) + patch_width;

fixedPoints = [
    0,            0;
    fixed_max_x,  0;
    fixed_max_x,  fixed_max_y;
    0,            fixed_max_y
    ];
movingPoints = corners;

tform = fitgeotrans(fixedPoints, movingPoints, 'projective');

% --- 6. Project, Re-distort, and Mask ---
disp('Calculating boundaries and mapping back to raw space...');
rawPixelIndices = cell(rows, cols);

safe_margin = 0.02;
numPtsPerEdge = 25;

% The linear scaling factor to achieve a 50% central area
area_scale = sqrt(0.5);

for r = 1:rows
    for c = 1:cols
        % --- Full Boundary Definition (For definition and plotting) ---
        x_left   = (c - 1) + safe_margin;
        x_right  = (c - 1) + patch_width - safe_margin;
        y_top    = (r - 1) + safe_margin;
        y_bot    = (r - 1) + patch_width - safe_margin;

        canonical_bounds_full = [
            x_left,  y_top;
            x_right, y_top;
            x_right, y_bot;
            x_left,  y_bot
            ];

        if showPlot
            % Map full bounds to undistorted space and plot
            undist_bounds = transformPointsForward(tform, canonical_bounds_full);
            plot(ax1, [undist_bounds(:,1); undist_bounds(1,1)], ...
                [undist_bounds(:,2); undist_bounds(1,2)], 'g-', 'LineWidth', 1);

            % Get distorted full boundaries for the raw image plot
            distPtsFull = getDistortedBoundary(canonical_bounds_full, numPtsPerEdge, tform, K, camIntrinsics);
            plot(ax2, [distPtsFull(:,1); distPtsFull(1,1)], ...
                [distPtsFull(:,2); distPtsFull(1,2)], 'c-', 'LineWidth', 1);
        end

        % --- Central 75% Extraction Boundary (For masking) ---
        cx = (x_left + x_right) / 2;
        cy = (y_top + y_bot) / 2;

        s_x_left  = cx - (cx - x_left) * area_scale;
        s_x_right = cx + (x_right - cx) * area_scale;
        s_y_top   = cy - (cy - y_top) * area_scale;
        s_y_bot   = cy + (y_bot - cy) * area_scale;

        canonical_bounds_extract = [
            s_x_left,  s_y_top;
            s_x_right, s_y_top;
            s_x_right, s_y_bot;
            s_x_left,  s_y_bot
            ];

        % Get distorted extraction boundaries
        distPtsExt = getDistortedBoundary(canonical_bounds_extract, numPtsPerEdge, tform, K, camIntrinsics);

        if showPlot
            plot(ax2, [distPtsExt(:,1); distPtsExt(1,1)], ...
                [distPtsExt(:,2); distPtsExt(1,2)], 'm--', 'LineWidth', 1.5);

            center_x = mean(distPtsExt(:,1));
            center_y = mean(distPtsExt(:,2));
            text(ax2, center_x, center_y, sprintf('%d,%d', r, c), ...
                'Color', 'yellow', 'HorizontalAlignment', 'center', 'FontSize', 8);
        end

        % Create mask and extract linear pixel indices using the 75% area
        checkMask = poly2mask(distPtsExt(:,1), distPtsExt(:,2), imageRows, imageCols);
        rawPixelIndices{r, c} = find(checkMask);
    end
end

if showPlot
    hold(ax1, 'off');
    hold(ax2, 'off');
    disp('Extraction complete. Figure left open for inspection.');
else
    disp('Extraction complete.');
end
end

% =========================================================================
% Local Helper Function
% =========================================================================
function distPts = getDistortedBoundary(bounds, numPtsPerEdge, tform, K, camIntrinsics)
% Maps canonical bounds to undistorted space, interpolates edges, and
% projects them back into the raw, distorted fisheye space.

undist_bounds = transformPointsForward(tform, bounds);

denseBoundaryUndist = [];
for i = 1:4
    p1 = undist_bounds(i, :);
    p2 = undist_bounds(mod(i, 4) + 1, :);

    x_edge = linspace(p1(1), p2(1), numPtsPerEdge)';
    y_edge = linspace(p1(2), p2(2), numPtsPerEdge)';

    denseBoundaryUndist = [denseBoundaryUndist; x_edge(1:end-1), y_edge(1:end-1)];
end

numPoints = size(denseBoundaryUndist, 1);
P_undist = [denseBoundaryUndist, ones(numPoints, 1)];
rays3D = P_undist / K;

distPts = worldToImage(camIntrinsics, eye(3), [0 0 0], rays3D);
end