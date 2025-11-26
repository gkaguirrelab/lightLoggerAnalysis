function [slopeMap, aucMap, spdByRegion, frq] = mapSPDs(video_hdf5_path, fps, window, step, options)
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
    options.chunk_size_seconds {mustBeNumeric} = 1;
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

% Calculate the step size in frames
framesPerChunk = floor(options.chunk_size_seconds * fps);

% Find the chunk starts
chunkStarts = 1:framesPerChunk/2:(end_frame-framesPerChunk);

% Find the number of chunks in the video
nChunks = length(chunkStarts);

% Find the row and column starts
rowStarts = 1:step:(nRows - window(1) + 1);
colStarts = 1:step:(nCols - window(2) + 1);
nRowPatches = length(rowStarts);
nColPatches = length(colStarts);

% Total number of patches (window positions across the image)
nPatches = nRowPatches * nColPatches;

% Allocate average Slope3D and AUC3D variables across all chunks
slopeMap = nan(nRows, nCols, nChunks);
aucMap = nan(nRows, nCols, nChunks);

% Allocate storage for the spds
spdByRegion = nan(nRowPatches, nColPatches, nChunks, framesPerChunk/2);

% Turn off a warning that occurs during robust linear fitting
warnState = warning();
warning('off','stats:statrobustfit:IterationLimit');

% Move over the chunks of the video
for ff = 1:nChunks

    % Counter for layer index (each patch corresponds to one layer)
    layer = 0;

    % Initialize 3D arrays for slope and AUC values per window layer
    slope3D = nan(nRows, nCols, nPatches);
    auc3D = nan(nRows, nCols, nPatches);

    % Inform the users
    tic;
    fprintf("Processing chunk: %d/%d...", ff, nChunks);

    % Initialize a frame chunk we wil populate
    frameChunk = h5read(video_hdf5_path, "/video", [1, 1, chunkStarts(ff)], [inf, inf, framesPerChunk]);
    frameChunk = permute(frameChunk, [3 2 1]);  % flip back to nFrames x nRows x nCols

    % Slide the analysis window across the image in row and column directions
    for rr = 1:nRowPatches
        for cc = 1:nColPatches

            % Get this row and column
            row = rowStarts(rr);
            col = colStarts(cc);

            % Increment patch counter
            layer = layer + 1;

            % Create a binary mask selecting the current region of interest
            regionMatrix = zeros(nRows, nCols);
            regionMatrix(row:row+window(1)-1, col:col+window(2)-1) = 1;
            try
                % Compute temporal SPD restricted to this region
                [spd, fLoc] = calcTemporalSPD(frameChunk, fps, 'lineResolution', false, 'regionMatrix', regionMatrix);
            catch
                % If computation failes, skip that patch
                continue;
            end

            % If spd or frq are nan, skip this patch
            if(any(isnan(spd)) || any(isnan(fLoc)))
                continue
            else
                frq = fLoc;
            end

            % Discard the zeroeth frequency
            frq = frq(2:end);
            spd = spd(2:end);

            % Save the raw spd
            spdByRegion(rr,cc,ff,:) = spd;

            % Fit a straight line in log-log space for the slope
            C = robustfit(log10(frq'), log10(spd) );

            % Calculate the auc
            auc = (mean(polyval(C,options.aucRange))/2)*diff(options.aucRange);

            % Assign slope and AUC values to current region layer
            slope3D(row:row+window(1)-1, col:col+window(2)-1, layer) = C(1);
            auc3D(row:row+window(1)-1, col:col+window(2)-1, layer) = auc;

        end % col

    end % row

    % Add this to the growing average for the whole video
    slopeMap(:,:,ff) = mean(slope3D, 3,'omitnan');
    aucMap(:,:,ff) =  mean(auc3D,3,'omitnan');

    % Finish the console report
    fprintf("%2.2f seconds\n", toc)

end % end chunks

% Restore the warning state
warning(warnState);

% Finally, finish the average slope map and AUC map calculations
slopeMap = mean(slopeMap,3,'omitmissing');
aucMap = mean(aucMap,3,'omitmissing');

% Get the mean spd by region
spdByRegion = squeeze(mean(spdByRegion,3,'omitmissing'));

figure;
imagesc(slopeMap,[-2.5 -2])
title(sprintf("Whole Video Mean | Slope Map"));
axis image;
colormap('hot');
colorbar;

figure;
imagesc(aucMap,[-3.75 -3.25])
title(sprintf("Whole Video Mean | AUC Map"));
axis image;
colormap('hot');
colorbar;
drawnow;

end