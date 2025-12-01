function [slopeMap, aucMap, spdByRegion, frq, medianImage] = mapSPDs(videoPath, options)
% Computes slope and Area Under the Curve (AUC) maps of temporal SPD across
% image regions and projects them onto a 1 m visual field surface, then
% plots both maps.
%
% Required Inputs:
%   video             - String. path to video in hdf5 format
%   fps               - Scalar. Sampling rate (Hz)
%   windowSpacePixels            - [height width] of square region
%   stepSpacePixels              - stepSpacePixels size for moving windowSpacePixels.
%   doPlot            - (boolean) Visualize the SPD maps or not
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
    videoPath {mustBeText}
    options.fps             (1,1) {mustBeNumeric}   = 120
    options.frameIdxToProcess {mustBeNumeric} = [1, inf];
    options.windowSpacePixels          (1,2) {mustBeNumeric}   = [24 24]
    options.stepSpacePixels            (1,1) {mustBeNumeric}   = 12
    options.windowTimeSecs {mustBeNumeric} = 10
    options.stepTimeSecs {mustBeNumeric} = 0.5
    options.aucFreqRangeHz {mustBeNumeric} = [log10(0.1),log10(60)]
    options.doPlot  (1,1) logical           = true
    options.nWorkers (1,1) {mustBeNumeric}   = 6
end

% Load in some info about the video to get us started
video_info = h5info(videoPath, "/video");
video_size = video_info.Dataspace.Size;
nRows = video_size(1);
nCols = video_size(2);
nFrames = video_size(3);

if(nFrames <= 0)
    error("Video has no frames");
end

% Place some options in variables
fps = options.fps;
windowSpacePixels = options.windowSpacePixels;
stepSpacePixels = options.stepSpacePixels;
windowTimeSecs = options.windowTimeSecs;
stepTimeSecs = options.stepTimeSecs;

startFrameIdx = options.frameIdxToProcess(1);
if(options.frameIdxToProcess(2) == inf)
    endFrameIdx = nFrames;
else
    endFrameIdx = options.frameIdxToProcess(2);
end

% Calculate the video chunk size in frames
framesPerChunk = floor(windowTimeSecs * fps);
framesPerStep = floor(stepTimeSecs * fps);

% Find the chunk starts
chunkStarts = startFrameIdx:framesPerStep:(endFrameIdx-framesPerChunk);

% Find the number of chunks in the video
nChunks = length(chunkStarts);

% Find the row and column starts
rowStarts = 1:stepSpacePixels:(nRows - windowSpacePixels(1) + 1);
colStarts = 1:stepSpacePixels:(nCols - windowSpacePixels(2) + 1);
nRowPatches = length(rowStarts);
nColPatches = length(colStarts);

% Total number of patches (windowSpacePixels positions across the image)
nPatches = nRowPatches * nColPatches;

% Allocate average Slope3D and AUC3D variables across all chunks
slopeMap = nan(nRows, nCols, nChunks);
aucMap = nan(nRows, nCols, nChunks);

% Allocate storage for the spds
spdByRegion = nan(nRowPatches, nColPatches, nChunks, floor(framesPerChunk/2));

% Allocate storage for the median image
medianImage = nan(nRows, nCols, nChunks);

% Allocate storage for the frq variable
frq = nan(1,floor(framesPerChunk/2)+1);

% Move over the chunks of the video
for ff = 1:nChunks

    % Turn off a warning that occurs during robust linear fitting
    warnState = warning();
    warning('off','stats:statrobustfit:IterationLimit');

    % Counter for layer index (each patch corresponds to one layer)
    layer = 0;

    % Initialize 3D arrays for slope and AUC values per windowSpacePixels layer
    slope3D = nan(nRows, nCols, nPatches);
    auc3D = nan(nRows, nCols, nPatches);

    % Inform the user
    startTime = datetime('now');
    fprintf("Processing chunk: %d/%d...", ff, nChunks);

    % Initialize a frame chunk we will populate
    frameChunk = h5read(videoPath, "/video", [1, 1, chunkStarts(ff)], [inf, inf, framesPerChunk]);
    frameChunk = permute(frameChunk, [3 2 1]);  % flip back to nFrames x nRows x nCols

    % Get the median image across time for this chunk
    medianImage(:,:,ff) = median(frameChunk,1,'omitmissing');

    % Slide the analysis windowSpacePixels across the image in row and column directions
    for rr = 1:nRowPatches
        for cc = 1:nColPatches

            % Get this row and column
            row = rowStarts(rr);
            col = colStarts(cc);

            % Increment patch counter
            layer = layer + 1;

            % Create a binary mask selecting the current region of interest
            regionMatrix = zeros(nRows, nCols);
            regionMatrix(row:row+windowSpacePixels(1)-1, col:col+windowSpacePixels(2)-1) = 1;
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

            % Save the raw spd
            spdByRegion(rr,cc,ff,:) = spd(2:end);

            % Fit a straight line in log-log space for the slope
            C = robustfit(log10(frq(2:end)'), log10(spd(2:end)) );

            % Calculate the auc
            auc = (mean(polyval(C,options.aucFreqRangeHz))/2)*diff(options.aucFreqRangeHz);

            % Assign slope and AUC values to current region layer
            slope3D(row:row+windowSpacePixels(1)-1, col:col+windowSpacePixels(2)-1, layer) = C(1);
            auc3D(row:row+windowSpacePixels(1)-1, col:col+windowSpacePixels(2)-1, layer) = auc;

        end % col
    end % row

    % Add this to the growing average for the whole video
    slopeMap(:,:,ff) = mean(slope3D, 3,'omitnan');
    aucMap(:,:,ff) =  mean(auc3D,3,'omitnan');

    % Finish the console report
    endTime = datetime('now');
    fprintf("%2.2f seconds\n", seconds(endTime-startTime))

    % Restore the warning state
    warning(warnState);

end % end chunks

% Drop the zeroeth frequency
frq = frq(2:end);

% Obtain the average slope and  AUC map, and the median image
slopeMap = mean(slopeMap,3,'omitmissing');
aucMap = mean(aucMap,3,'omitmissing');
medianImage = median(medianImage,3,'omitmissing');

% Get the mean spd by region
spdByRegion = squeeze(mean(spdByRegion,3,'omitmissing'));

% Show the results if requested
if options.doPlot

    figure('Position',[100,100,900,300]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile
    imagesc(slopeMap)
    title(sprintf("Slope Map"));
    axis image;
    colormap('hot');
    colorbar;

    nexttile
    imagesc(aucMap)
    title(sprintf("AUC Map"));
    axis image;
    colormap('hot');
    colorbar;
    drawnow;

    nexttile
    title('Median image');
    image(medianImage);
    axis image;
    drawnow;

end

end