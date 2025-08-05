function [X, Y, Z, slopeMap, interceptMap, frq, hFigSlope, hFigIntercept] = mapSlopeSPD(v, fps, window, step, fisheyeIntrinsics)
% Computes slope and intercept maps of temporal SPD across image regions
% and projects them onto a 1 m visual field surface, then plots both maps.
%
% Inputs:
%   v                   - [frames x rows x cols] video chunk (double)
%   fps                 - Sampling rate (Hz)
%   window              - [height width] of square region (e.g., [40 40])
%   step                - Step size for moving window
%   fisheyeIntrinsics   - camera.FisheyeIntrinsics object
%
% Outputs:
%   X, Y, Z             - 3D coordinates (rows x cols) on hemisphere
%   slopeMap            - local 1/f slope per pixel
%   interceptMap        - local y-intercept of log-log fit per pixel
%   frq                 - frequency vector used in SPD

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

[xGrid, yGrid] = meshgrid(1:nCols, 1:nRows);
cx = fisheyeIntrinsics.DistortionCenter(1);
cy = fisheyeIntrinsics.DistortionCenter(2);
xC = xGrid - cx; yC = yGrid - cy;

r = sqrt(xC.^2 + yC.^2);
coeffs = fisheyeIntrinsics.MappingCoefficients;
a0=coeffs(1); a2=coeffs(2); a3=coeffs(3); a4=coeffs(4);

theta = nan(size(r));
for k=1:numel(r)
    func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - r(k);
    theta(k) = fzero(func,[0,pi]);
end

phi = atan2(yC,xC); R = 1;
X = reshape(R*sin(theta).*cos(phi), nRows, nCols);
Y = reshape(R*sin(theta).*sin(phi), nRows, nCols);
Z = reshape(R*cos(theta),         nRows, nCols);

% Plot slope map
hFigSlope = figure;
surf(X, Y, Z, slopeMap, 'EdgeColor','none'); shading interp; lighting none;
axis equal off; colormap jet; colorbar; title('1/f SPD Slope Map');

% Plot intercept map
hFigIntercept = figure;
surf(X, Y, Z, interceptMap, 'EdgeColor','none'); shading interp; lighting none;
axis equal off; colormap jet; colorbar; title('1/f SPD Intercept Map');

end