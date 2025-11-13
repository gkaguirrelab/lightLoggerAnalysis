function [slopeMap, aucMap, frq] = mapSPDs(video, fps, window, step, options)
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
        video               {mustBeText}
        fps             (1,1) {mustBeNumeric}   = 120
        window          (1,2) {mustBeNumeric}   = [40 40] 
        step            (1,1) {mustBeNumeric}   = 20
        options.doPlot  (1,1) logical           = false
        options.theta   {mustBeNumeric}         = []
        options.phi     {mustBeNumeric}         = []
        options.R       (1,1) {mustBeNumeric}   = 1
        options.num_frames_to_process {mustBeNumeric} = [1, inf];
        options.chunk_size_seconds {mustBeNumeric} = 15; 
    end 
        % Open a reader to the video 
        world_reader = videoIOWrapper(v); 
        nRows = world_reader.Height; 
        nCols = world_reader.Width; 

        % Compute how many window positions fit vertically and horizontally
        maxRows = floor((nRows - window(1)) / step) + 1;
        maxCols = floor((nCols - window(2)) / step) + 1;
        
        % Total number of patches (window positions across the image)
        total_patches = maxRows * maxCols;
        % Initialize 3D arrays for slope and AUC values per window layer
        slope3D     = nan(nRows, nCols, total_patches);
        auc3D       = nan(nRows, nCols, total_patches); % *** CHANGE 1: Renamed intercept3D to auc3D ***
        
        % Define the start and endpoint that we want to analyze of the video 
        start = options.num_frames_to_process(1); 
        endpoint = options.num_frames_to_process(2);
        if(endpoint == inf)
            endpoint = world_frame_reader.NumFrames; 
        end 

        % Calculate the step size in frames 
        chunk_size_frames = floor(chunk_size_seconds * world_reader.FrameRate); 
        
        % Counter for layer index (each patch corresponds to one layer)
        layer = 0;

        % Move over the chunks of the video
        for frame_num = start:chunk_size_frames:endpoint
            % Find the local chunk size (e.g. the last frame may not be the ideal large) 
            local_chunk_size_frames = min(chunk_size_frames, world_reader.NumFrames - frame_num + 1); 
            
            % Initialize a frame chunk we wil populate
            frame_chunk = zeros(local_chunk_size_frames, nRows, nCols); 

            % Populate the frame chunk in 30 second chunks 
            for insertion_index = 1:local_chunk_size_frames
                frame_chunk(insertion_index, :, :) = world_reader.readFrame("frameNum", frame_num+insertion_index - 1, "grayscale", true); 
            end

            % Find the SPD of the full spatial resolution
            [~, frq] = calcTemporalSPD(frame_chunk, fps, 'lineResolution', false);

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
                        [spd, fLoc] = calcTemporalSPD(local_chunk_size_frames, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
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
                    
                    % Only proceed if at least 2 valid frequency bins remain
                    if nnz(valid)>=2
                            
                        f_min = min(fLoc(valid));

                        f_max = max(fLoc(valid));
                        
                        % Fit a straight line in log-log space for the slope
                        C = polyfit(log10(fLoc(valid)), log10(spd(valid)), 1);
                        % Calculate Area Under the *FITTED LINE* (Analytic Integral)
                        % The function is SPD(f) = 10^C(2) * f^C(1)
                        % Check for the special case where the exponent C(1) is -1
                        if abs(C(1) + 1) < 1e-6
                            % Integral of 1/f is ln(f)
                            auc = (10^C(2)) * (log(f_max) - log(f_min));
                        else
                            % General case: Integral[ A * f^B df ]
                            A = 10^C(2); % Scaling factor
                            B = C(1);    % Exponent/Slope

                            % Closed-form integral: (A / (B+1)) * (f_max^(B+1) - f_min^(B+1))
                            auc = (A / (B + 1)) * ( (f_max^(B + 1)) - (f_min^(B + 1)) );
                        end

                    % Assign slope and AUC values to current region layer
                    slope3D(row:row+window(1)-1, col:col+window(2)-1, layer) = C(1);
                    auc3D(row:row+window(1)-1, col:col+window(2)-1, layer) = auc;

                end
            end % Added missing 'end' for the row loop closure
        
            % Average slope across all overlapping layers (ignoring NaNs)
            slopeMap     = mean(slope3D,     3, 'omitnan');
            % Average AUC across all overlapping layers (ignoring NaNs)
            aucMap       = mean(auc3D, 3, 'omitnan'); % 
            
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
                axis equal; colormap jet; colorbar; title('1/f SPD Area Under Curve Map'); % 
            end

        end % end chunks
        
    end % end function
end 