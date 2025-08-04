function [X, Y, Z, slopeColor, frq] = mapSlopeSPD(v, fps, window, step, fisheyeIntrinsics)
% Computes a slope map of temporal SPD across image regions and projects
% it onto a 1 m visual field surface consistent with the camera's fisheye mapping.
%
% Inputs:
%   v                   - [frames x rows x cols] video chunk (double)
%   fps                 - Sampling rate (Hz). Defaults to 120 Hz in
%                         calcTemporalSPD.m
%   window              - [height width] of square region (e.g., [40 40])
%   step                - Step size for moving window
%   fisheyeIntrinsics   - camera.FisheyeIntrinsics object
%
% Outputs:
%   X, Y, Z             - 3D matrices (rows x cols) on hemisphere
%   slopeColor          - Slope values mapped as color on surface
%   frq                 - Frequency vector 
% 
% Usage:
% [X,Y,Z,slopeColor,frq] = mapSlopeSPD(v, 120, [40,40], 20, fisheyeIntrinsics)

[~, nRows, nCols] = size(v);

% Prepare for window scanning
maxRows = floor((nRows - window(1)) / step) + 1;
maxCols = floor((nCols - window(2)) / step) + 1;
total_patches = maxRows * maxCols;
fprintf('Processing %d patches\n', total_patches);

% Preallocate slope map layers for each patch
slopeMap3D = nan(nRows, nCols, total_patches);

% Get frequency bins for entire chunk
[~, frq] = calcTemporalSPD(v, fps, 'lineResolution', false);
frq_bins = numel(frq);

% Loop over patches
layer = 0;
for row = 1:step:(nRows - window(1) + 1)
    for col = 1:step:(nCols - window(2) + 1)
        layer = layer + 1;
        regionMatrix = zeros(nRows, nCols);
        regionMatrix(row:row+window(1)-1, col:col+window(2)-1) = 1;
        try
            [spd, frq] = calcTemporalSPD(v, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
        catch ME
            warning('Skipping patch at (%d,%d): %s', row, col, ME.message);
            continue;
        end

        spd = spd(:);
        frq = frq(:);
        % Censor exactly 30 Hz
        spd(frq == 30) = NaN;

        % Mask out invalid frequencies and noise bands
        validIdx = frq>0 & spd>0;
        validIdx(frq>=52 & frq<=71) = false;

        if nnz(validIdx) >= 2
            coeffs = polyfit(log10(frq(validIdx)), log10(spd(validIdx)), 1);
            slope = coeffs(1);
        else
            slope = NaN;
        end
        slopeMap3D(row:row+window(1)-1, col:col+window(2)-1, layer) = slope;
    end
end

% Average slope across overlapping patches
slopeMap = mean(slopeMap3D, 3, 'omitnan');

% Project to visual field
[xGrid, yGrid] = meshgrid(1:nCols, 1:nRows);
cx = fisheyeIntrinsics.DistortionCenter(1);
cy = fisheyeIntrinsics.DistortionCenter(2);
xCentered = xGrid - cx;
yCentered = yGrid - cy;
rPixels = sqrt(xCentered.^2 + yCentered.^2);
mapCoeffs = fisheyeIntrinsics.MappingCoefficients;
a0 = mapCoeffs(1); a2 = mapCoeffs(2); a3 = mapCoeffs(3); a4 = mapCoeffs(4);
theta = nan(size(rPixels));
for k = 1:numel(rPixels)
    r_val = rPixels(k);
    func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - r_val;
    theta(k) = fzero(func, [0, pi]);
end
phi = atan2(yCentered, xCentered);
R = 1;
X = R * sin(theta) .* cos(phi);
Y = R * sin(theta) .* sin(phi);
Z = R * cos(theta);
X = reshape(X, nRows, nCols);
Y = reshape(Y, nRows, nCols);
Z = reshape(Z, nRows, nCols);

% Display
gcf;
surf(X, Y, Z, slopeMap, 'EdgeColor', 'none');
axis equal off; colormap jet; colorbar;
title('1/f SPD Slope Map'); view(3); camlight; lighting gouraud;

slopeColor = slopeMap;
frq = frq;
end