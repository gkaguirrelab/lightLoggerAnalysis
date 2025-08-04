function [X, Y, Z, slopeColor, frq] = mapSlopeSPD(v, fps, window, step, channel, fisheyeIntrinsics)
% Computes a slope map of temporal SPD across image regions and projects
% it onto a 1 m visual field surface consistent with the camera's fisheye mapping.
%
% Inputs:
%   v                   - [frames x rows x cols] video chunk (double)
%   fps                 - Sampling rate (Hz). Defaults to 120 Hz in
%                         calcTemporalSPD.m
%   window              - [height width] of square region (e.g., [40 40])
%   step                - Step size for moving window
%   channel             - 'LM', 'L-M', or 'S'
%   fisheyeIntrinsics   - camera.FisheyeIntrinsics object
%
% Outputs:
%   X, Y, Z             - 3D matrices (rows x cols) on hemisphere
%   slopeColor          - Slope values mapped as color on surface
%   frq                 - Frequency vector 
% 
% Usage:
% v = chunks{1}.W.v
% data = load('~/FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/camera_intrinsics_calibration.mat');
% fisheyeIntrinsics = data.camera_intrinsics_calibration.camera_intrinsics.Intrinsics; 
% [X, Y, Z, slopeColor, frq] = mapSlopeSPD(v, 120, [40, 40], 20, 'LM', fisheyeIntrinsics)

[~, nRows, nCols] = size(v);

% Prepare for window scanning
maxRows = floor((nRows - window(1)) / step) + 1;
maxCols = floor((nCols - window(2)) / step) + 1;
total_patches = maxRows * maxCols;
fprintf('Processing %d patches\n', total_patches);

% Preallocate slope map layers for each patch
slopeMap3D = nan(nRows, nCols, total_patches);

% Get frequency bins
[~, frq] = calcTemporalSPD(v, fps, false, 'postreceptoralChannel', channel);
frq_bins = numel(frq);

% Initialize accumulators for SPD
accumSPD = zeros(nRows, nCols, frq_bins);
countSPD = zeros(nRows, nCols);

layer = 0;
for row = 1:step:(nRows - window(1) + 1)
    for col = 1:step:(nCols - window(2) + 1)
        layer = layer + 1;

        regionMatrix = zeros(nRows, nCols);
        regionMatrix(row:row+window(1)-1, col:col+window(2)-1) = 1;

        try
            [spd, frq] = calcTemporalSPD(v, fps, false, 'postreceptoralChannel', channel, 'regionMatrix', regionMatrix);
        catch ME
            warning('Skipping patch at (%d, %d): %s', row, col, ME.message);
            continue
        end

        spd = spd(:);
        frq = frq(:);
        validIdx = frq > 0 & spd > 0;
        validIdx(frq >= 52 & frq <= 71) = false;

        if nnz(validIdx) >= 2
            coeffs = polyfit(log10(frq(validIdx)), log10(spd(validIdx)), 1);
            slope = coeffs(1);
        else
            slope = NaN;
        end

        slopeMap3D(row:row+window(1)-1, col:col+window(2)-1, layer) = slope;

        for f = 1:frq_bins
            accumSPD(row:row+window(1)-1, col:col+window(2)-1, f) = ...
                accumSPD(row:row+window(1)-1, col:col+window(2)-1, f) + spd(f);
        end
        countSPD(row:row+window(1)-1, col:col+window(2)-1) = ...
            countSPD(row:row+window(1)-1, col:col+window(2)-1) + 1;
    end
end

slopeMap = mean(slopeMap3D, 3, 'omitnan');

%% PROJECT TO VISUAL FIELD (1m) CONSISTENT WITH mapFisheyeVideoChunkToVisualField

[imgHeight, imgWidth] = deal(nRows, nCols);
[xGrid, yGrid] = meshgrid(1:imgWidth, 1:imgHeight);

cx = fisheyeIntrinsics.DistortionCenter(1);
cy = fisheyeIntrinsics.DistortionCenter(2);
xCentered = xGrid - cx;
yCentered = yGrid - cy;
rPixels = sqrt(xCentered.^2 + yCentered.^2);

mappingCoefficients = fisheyeIntrinsics.MappingCoefficients;
a0 = mappingCoefficients(1);
a2 = mappingCoefficients(2);
a3 = mappingCoefficients(3);
a4 = mappingCoefficients(4);

theta = zeros(size(rPixels));
for i = 1:numel(rPixels)
    r_val = rPixels(i);
    func = @(theta_val) a0*theta_val + a2*theta_val^3 + ...
                         a3*theta_val^4 + a4*theta_val^5 - r_val;
    try
        theta(i) = fzero(func, [0, pi]);
    catch
        theta(i) = NaN;
    end
end

phi = atan2(yCentered, xCentered);

R = 1; % 1 meter
X = R * sin(theta) .* cos(phi);
Y = R * sin(theta) .* sin(phi);
Z = R * cos(theta);

X = reshape(X, imgHeight, imgWidth);
Y = reshape(Y, imgHeight, imgWidth);
Z = reshape(Z, imgHeight, imgWidth);

%% Display

figure;
surf(X, Y, Z, slopeMap, 'EdgeColor', 'none');
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title(sprintf('Local 1/f SPD Slope Map (%s) on 1m Visual Field', channel));
colormap jet;
colorbar;
view(3);
camlight;
lighting gouraud;
axis off;

slopeColor = slopeMap;

end
