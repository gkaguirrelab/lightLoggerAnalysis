function [whole_video_mean_slopeMap, whole_video_mean_aucMap, whole_video_frequencies] = mapSPDs(video_hdf5_path, fps, window, step, options)
% Computes slope and Area Under the Curve (AUC) maps of temporal SPD across
% image regions and projects them onto a 1 m visual field surface, then
% plots both maps.
%
% Required Inputs:
%   video             - String. path to video
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
    [slopeMap, aucMap, frq] = mapSlopeAUCSPD(video, fps, [40 40], 20, True, theta, phi, R)
%}
arguments
    video_hdf5_path {mustBeText}
    fps             (1,1) {mustBeNumeric}   = 120
    window          (1,2) {mustBeNumeric}   = [24 24]
    step            (1,1) {mustBeNumeric}   = 12
    options.doPlot  (1,1) logical           = false
    options.num_frames_to_process {mustBeNumeric} = [1, inf];
    options.chunk_size_seconds {mustBeNumeric} = 5;
    options.aucRange {mustBeNumeric} = [log10(0.1),log10(60)];
end

% Load in some info about the video to get us started
video_info = h5info(video_hdf5_path, "/video");
video_size = video_info.Dataspace.Size;
nRows = video_size(1);
nCols = video_size(2);
nFrames = video_size(3);

if(nFrames <= 0)
    error("Video has no frames");
end

% Define the start and endpoint that we want to analyze of the video
start_frame = options.num_frames_to_process(1);
end_frame = options.num_frames_to_process(2);

if(end_frame == inf)
    end_frame = nFrames;
end

% Compute how many window positions fit vertically and horizontally
maxRows = floor((nRows - window(1)) / step) + 1;
maxCols = floor((nCols - window(2)) / step) + 1;

% Total number of patches (window positions across the image)
total_patches = maxRows * maxCols;

% Calculate the step size in frames
chunk_size_frames = floor(options.chunk_size_seconds * fps);

% Find the chunk starts
chunkStarts = 1:chunk_size_frames/2:(end_frame-chunk_size_frames);

% Find the number of chunks in the video
nChunks = length(chunkStarts);

% Allocate average Slope3D and AUC3D variables across all chunks
whole_video_mean_slopeMap = zeros(nRows, nCols);
whole_video_mean_aucMap = zeros(nRows, nCols);

% Move over the chunks of the video
current_chunk = 1;
for frame_num = chunkStarts

    % Counter for layer index (each patch corresponds to one layer)
    layer = 0;

    % Initialize 3D arrays for slope and AUC values per window layer
    slope3D = nan(nRows, nCols, total_patches);
    auc3D = nan(nRows, nCols, total_patches);

    tic;
    fprintf("Processing chunk: %d/%d\n", current_chunk, nChunks);

    % Find the local chunk size (e.g. the last frame may not be the ideal large)
    local_chunk_size_frames = min(chunk_size_frames, nFrames - frame_num + 1);

    if(local_chunk_size_frames == 0)
        disp("Local chunk size is 0");
    end

    % Initialize a frame chunk we wil populate
    frame_chunk = h5read(video_hdf5_path, "/video", [1, 1, frame_num], [inf, inf, local_chunk_size_frames]);
    frame_chunk = permute(frame_chunk, [3 1 2]);  % flip back to nFrames x nRows x nCols

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
                [spd, fLoc] = calcTemporalSPD(frame_chunk, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
            catch
                % If computation failes, skip that patch
                continue;
            end

            % If spd or fLoc are nan, skip this patch
            if(any(isnan(spd)) || any(isnan(fLoc)))
                continue
            end

            % Discard the zeroeth frequency
            fLoc = fLoc(2:end);
            spd = spd(2:end);

            % Fit a straight line in log-log space for the slope
            C = polyfit(log10(fLoc'), log10(spd), 1);

            % Calculate the auc
            auc = (mean(polyval(C,options.aucRange))/2)*diff(options.aucRange);

            % Assign slope and AUC values to current region layer
            slope3D(row:row+window(1)-1, col:col+window(2)-1, layer) = C(1);
            auc3D(row:row+window(1)-1, col:col+window(2)-1, layer) = auc;

        end

    end % Added missing 'end' for the row loop closure

    % Average slope across all overlapping layers (ignoring NaNs)
    slopeMap = mean(slope3D,     3, 'omitnan');
    % Average AUC across all overlapping layers (ignoring NaNs)
    aucMap = mean(auc3D, 3, 'omitnan'); %

    % Add this to the growing average for the whole video
    whole_video_mean_slopeMap(:,:,current_chunk) = slopeMap;
    whole_video_mean_aucMap(:,:,current_chunk) =  aucMap;

    % If plotting is requested
    if options.doPlot
        figure;
        imagesc(slopeMap);
        title(sprintf("Chunk %d | Slope Map", current_chunk));
        axis image;
        colormap('hot');
        colorbar;

        figure;
        imagesc(aucMap);
        title(sprintf("Chunk %d | AUC Map", current_chunk));
        axis image;
        colormap('hot');
        colorbar;
        drawnow;

    end

    % Calculate the elapsed tiem per chunk
    elapsed_seconds = toc;

    fprintf("Chunk %d took: %f seconds\n", current_chunk, elapsed_seconds)

    % Increment the chunk number we are on
    current_chunk = current_chunk + 1;

end % end chunks

% Finally, finish the average slope map and AUC map calculations
whole_video_mean_slopeMap = mean(whole_video_mean_slopeMap,3,'omitmissing');
whole_video_mean_aucMap = mean(whole_video_mean_aucMap,3,'omitmissing');

figure;
imagesc(whole_video_mean_slopeMap);
title(sprintf("Whole Video Mean | Slope Map"));
axis image;
colormap('hot');
colorbar;

figure;
imagesc(whole_video_mean_aucMap);
title(sprintf("Whole Video Mean | AUC Map"));
axis image;
colormap('hot');
colorbar;
drawnow;

end