function [slopeMap, aucMap, frq] = mapSPDs(v, fps, window, step, options)
% Computes slope and Area Under the Curve (AUC) maps of temporal SPD across 
% image regions and projects them onto a 1 m visual field surface, then 
% plots both maps.
%
% Required Inputs:
%   v                 - [frames x rows x cols] video chunk (double)
%   fps               - Sampling rate (Hz)
%   window            - [height width] of square region. Defaults to [40 40]
%   step              - Step size for moving window. Defaults to 20
%   doPlot            - (boolean) Visualize the SPD maps or not
%   theta             - (radians) [rows x cols] elevation-from-optical-axis
%   phi               - (radians) [rows x cols] azimuth                    
%   R                 - (radians) (scalar) radius                       
%
% Optional:
%   affineMat         - 2x3 affine matrix to apply to [az, el] before XYZ
%
% Outputs:
%   slopeMap          - local 1/f slope per pixel
%   aucMap            - local Area Under the Curve (integrated power) per pixel
%   frq               - frequency vector used in SPD
%   hFigSlope         - figure handle for slope map
%   hFigAUC           - figure handle for AUC map
%
% Example Usage: 
%{
    [slopeMap, aucMap, frq] = mapSlopeAUCSPD(v, fps, [40 40], 20, True, theta, phi, R)
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
    % Retrieve the size of the video chunk 
    [~, nRows, nCols] = size(v);
    % Compute how many window positions fit vertically and horizontally
    maxRows = floor((nRows - window(1)) / step) + 1;
    maxCols = floor((nCols - window(2)) / step) + 1;
    
    % Total number of patches (window positions across the image)
    total_patches = maxRows * maxCols;
    % Initialize 3D arrays for slope and AUC values per window layer
    slope3D     = nan(nRows, nCols, total_patches);
    auc3D       = nan(nRows, nCols, total_patches); % *** CHANGE 1: Renamed intercept3D to auc3D ***
    % Calculate the temporal SPD of the input chunk
    [~, frq] = calcTemporalSPD(v, fps, 'lineResolution', false);
    % Counter for layer index (each patch corresponds to one layer)
    layer = 0;
    % Slide the analysis window across the image in row and column directions
    for row = 1:step:(nRows - window(1) + 1)
        for col = 1:step:(nCols - window(2) + 1)
            % Increment patch counter
            layer = layer + 1;
            
            % Create a binary mask selecting the current region of interest
            regionMatrix = zeros(nRows, nCols);
            regionMatrix(row:row+window(1)-1, col:col+window(2)-1) = 1;
            try
                % Compute temporal SPD restricted to this region
                [spd, fLoc] = calcTemporalSPD(v, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
            catch
                % If computation failes, skip that patch
                continue;
            end
            % Flatten vectors for fitting/AUC calculation
            spd = spd(:); 
            fLoc = fLoc(:);
            
            % Exclude exactly 30 Hz (set to NaN, e.g. avoid mains noise artifact)
            spd(fLoc==30) = NaN;
            
            % Define valid data points: positive frequencies and positive power
            valid = fLoc>0 & spd>0;
            % Exclude a specific band of frequencies (52–71 Hz)
            valid(fLoc>=52 & fLoc<=71) = false;
            % Only proceed if at least 2 valid frequency bins remain
            if nnz(valid)>=2
                
                % Calculate AUC using trapezoidal integration (trapz)
                % AUC is the sum of the power spectrum (SPD) multiplied by the 
                % frequency bin width (df), which is handled by trapz.
                % The integration is performed only over the valid frequencies.
                auc = trapz(fLoc(valid), spd(valid)); 
                
                % Fit a straight line in log-log space for the slope
                C = polyfit(log10(fLoc(valid)), log10(spd(valid)), 1);
                
                % Assign slope and AUC values to current region layer
                slope3D(row:row+window(1)-1, col:col+window(2)-1, layer)     = C(1);
                auc3D(row:row+window(1)-1, col:col+window(2)-1, layer)       = auc; % *** CHANGE 3: Store AUC instead of intercept ***
            end
        end
    end % Added missing 'end' for the row loop closure
    
    % Average slope across all overlapping layers (ignoring NaNs)
    slopeMap     = mean(slope3D,     3, 'omitnan');
    % Average AUC across all overlapping layers (ignoring NaNs)
    aucMap       = mean(auc3D, 3, 'omitnan'); % *** CHANGE 4: Renamed interceptMap to aucMap ***
    
    % If plotting is requested
    if options.doPlot
        % Ensure angular coordinates are provided
        assert(~isempty(options.theta) && ~isempty(options.phi), ...
            'To plot, provide options.theta and options.phi.');
        
        % Extract radius for projection
        R = options.R;
        % Convert spherical coordinates (R, theta, phi) to Cartesian (X,Y,Z)
        X = reshape(R .* sin(options.theta) .* cos(options.phi), nRows, nCols);
        Y = reshape(R .* sin(options.theta) .* sin(options.phi), nRows, nCols);
        Z = reshape(R .* cos(options.theta),                     nRows, nCols);
        % Plot slope map on the visual field surface
        figure;
        surf(X, Y, Z, slopeMap, 'EdgeColor','none'); shading interp; lighting none;
        axis equal; colormap jet; colorbar; title('1/f SPD Slope Map');
        % Plot AUC map on the visual field surface
        figure;
        surf(X, Y, Z, aucMap, 'EdgeColor','none'); shading interp; lighting none;
        axis equal; colormap jet; colorbar; title('1/f SPD Area Under Curve Map'); % *** CHANGE 5: Updated plot title ***
    end
end
