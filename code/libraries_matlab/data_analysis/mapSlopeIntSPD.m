function [slopeMap, interceptMap, frq] = ...
    mapSlopeIntSPD(v, fps, window, step, options)
% Computes slope and intercept maps of temporal SPD across image regions
% and projects them onto a 1 m visual field surface, then plots both maps.
%
% Required Inputs:
%   v                 - [frames x rows x cols] video chunk (double)
%   fps               - Sampling rate (Hz)
%   window            - [height width] of square region. Defaults to [40 40]
%   step              - Step size for moving window. Defaults to 20
%   doPlot            - (boolean) Visualize the SPD maps or not
%   theta             - (radians) [rows x cols] elevation-from-optical-axis  (optional)
%   phi               - (radians) [rows x cols] azimuth                      (optional)
%   R                 - (radians) (scalar) radius                            (optional)
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
%
% Example Usage: 
%{
    [slopeMap, interceptMap, frq] = mapSlopeIntSPD(v, fps, [40 40], 20, True, theta, phi, R)
%}

arguments
    v               (:,480,640) {mustBeNumeric}
    fps             (1,1) {mustBeNumeric}   = 120
    window          (1,2) {mustBeNumeric}   = [40 40] 
    step            (1,1) {mustBeNumeric}   = 20
    options.doPlot  (1,1) logical           = false
    options.theta   {mustBeNumeric}         = []
    options.phi     {mustBeNumeric}         = []
    options.R       (1,1) {mustBeNumeric}   = 1
end


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

if nargin < 5 || ~option.doPlot
end

if nargin >= 5
    if options.doPlot && nargin < 6
        error("Must provide theta, phi, and R values to visualize.");
    else 
        X = reshape(R .*sin(options.theta) .*cos(options.phi), nRows, nCols);
        Y = reshape(R .*sin(options.theta) .*sin(options.phi), nRows, nCols);
        Z = reshape(R .*cos(options.theta),         nRows, nCols);
        
        % Plot slope map
        figure;
        surf(X, Y, Z, slopeMap, 'EdgeColor','none'); shading interp; lighting none;
        axis equal; colormap jet; colorbar; title('1/f SPD Slope Map');
        
        % Plot intercept map
        figure;
        surf(X, Y, Z, interceptMap, 'EdgeColor','none'); shading interp; lighting none;
        axis equal; colormap jet; colorbar; title('1/f SPD Intercept Map');    
    end
end

end

