function [slopeMap, interceptMap, frq] = ...
    mapSlopeIntSPD(v, fps, window, step, theta, phi, R)
% Computes slope and intercept maps of temporal SPD across image regions
% and projects them onto a 1 m visual field surface, then plots both maps.
%
% Required Inputs:
%   v                 - [frames x rows x cols] video chunk (double)
%   fps               - Sampling rate (Hz)
%   window            - [height width] of square region (e.g., [40 40])
%   step              - Step size for moving window
%   theta             - (radians) [rows x cols] elevation-from-optical-axis  (optional)
%   phi               - (radians) [rows x cols] azimuth                      (optional)
%
%
% Optional:
%   affineMat         - 2x3 affine matrix to apply to [az, el] before XYZ
%
% Outputs:
%   slopeMap          - local 1/f slope per pixel
%   interceptMap      - local y-intercept of log-log fit per pixel
%   frq               - frequency vector used in SPD
%   hFigSlope         - figure handle for slope map
%   hFigIntercept     - figure handle for intercept map

[~, nRows, nCols] = size(v);

maxRows = floor((nRows - window(1)) / step) + 1;
maxCols = floor((nCols - window(2)) / step) + 1;
total_patches = maxRows * maxCols;
slope3D     = nan(nRows, nCols, total_patches);
intercept3D = nan(nRows, nCols, total_patches);

[~, frq] = calcTemporalSPD(v, fps, 'lineResolution', false);

layer = 0;
for row = 1:step:(nRows - window(1) + 1)
    for col = 1:step:(nCols - window(2) + 1)
        layer = layer + 1;
        regionMatrix = zeros(nRows, nCols);
        regionMatrix(row:row+window(1)-1, col:col+window(2)-1) = 1;
        try
            [spd, fLoc] = calcTemporalSPD(v, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
        catch
            continue;
        end
        spd = spd(:); fLoc = fLoc(:);
        spd(fLoc==30) = NaN;
        valid = fLoc>0 & spd>0;
        valid(fLoc>=52 & fLoc<=71) = false;
        if nnz(valid)>=2
            C = polyfit(log10(fLoc(valid)), log10(spd(valid)), 1);
            slope3D(row:row+window(1)-1, col:col+window(2)-1, layer)     = C(1);
            intercept3D(row:row+window(1)-1, col:col+window(2)-1, layer) = C(2);
        end
    end
end

slopeMap     = mean(slope3D,     3, 'omitnan');
interceptMap = mean(intercept3D, 3, 'omitnan');

X = reshape(R .*sin(theta) .*cos(phi), nRows, nCols);
Y = reshape(R .*sin(theta) .*sin(phi), nRows, nCols);
Z = reshape(R .*cos(theta),         nRows, nCols);

% Plot slope map
hFigSlope = figure;
surf(X, Y, Z, slopeMap, 'EdgeColor','none'); shading interp; lighting none;
axis equal; colormap jet; colorbar; title('1/f SPD Slope Map');

% Plot intercept map
hFigIntercept = figure;
surf(X, Y, Z, interceptMap, 'EdgeColor','none'); shading interp; lighting none;
axis equal; colormap jet; colorbar; title('1/f SPD Intercept Map');

end