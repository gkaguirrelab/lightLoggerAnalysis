function [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = mapSPDs(videoPath, options)
% Compute temporal SPD exponent and variance maps across image regions
%
% Syntax:
%   [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = mapSPDs(videoPath)
%   [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = mapSPDs(videoPath, options)
%
% Description:
%   This function computes temporal spectral power distribution (SPD)
%   summaries across local image regions of a video. The video is divided
%   into overlapping temporal chunks and overlapping spatial patches. For
%   each patch within each chunk, the function estimates the temporal SPD
%   and fits a line in the log-log domain to characterize the slope of the
%   spectrum. The fitted slope is stored as the exponent map, and the area
%   under the fitted spectrum is used to derive a variance map.
%
%   The function also stores the regional SPDs, the frequency support used
%   in the spectral analysis, the median image across time, and a vector
%   indicating dropped or invalid frames. After processing all chunks, the
%   function averages the exponent, variance, and SPD summaries across
%   time. The resulting maps can optionally be displayed.
%
% Inputs:
%   videoPath            - Char/string. Path to the video file to be
%                          analyzed.
%
% Optional key/value pairs:
%   fps                  - Scalar. Sampling rate of the video in Hz.
%   frameIdxToProcess    - Two-element numeric vector. Start and end frame
%                          indices to process. Use [1, inf] to process the
%                          full video.
%   windowSpacePixels    - Two-element numeric vector. Height and width of
%                          the spatial analysis window in pixels.
%   stepSpacePixels      - Scalar. Step size in pixels between adjacent
%                          spatial analysis windows.
%   windowTimeSecs       - Scalar. Duration in seconds of each temporal
%                          chunk used for SPD estimation.
%   stepTimeSecs         - Scalar. Step size in seconds between adjacent
%                          temporal chunks.
%   doPlot               - Logical. If true, display the exponent map,
%                          variance map, and median image at the end of
%                          processing.
%   nWorkers             - Scalar. Number of workers requested for
%                          parallel processing. Included for compatibility
%                          with parallelized versions of this function.
%   frameDropVector      - Numeric vector. Optional logical/numeric vector
%                          indicating frames that should be treated as
%                          dropped and replaced with NaNs before analysis.
%
% Outputs:
%   exponentMap          - Numeric matrix. Spatial map of the fitted log-log
%                          temporal SPD slope, averaged across chunks.
%   varianceMap          - Numeric matrix. Spatial map of the variance
%                          derived from the fitted temporal SPD, averaged
%                          across chunks.
%   spdByRegion          - Numeric array. Mean temporal SPD for each
%                          spatial patch, averaged across chunks.
%   frq                  - Numeric vector. Frequency support used in the
%                          SPD computation, excluding the zeroeth and
%                          Nyquist frequencies.
%   medianImage          - Numeric matrix. Median image across processed
%                          chunks.
%   frameDropVector      - Numeric vector. Updated frame drop vector
%                          indicating frames that were treated as invalid.
%
% Examples:
%{
    videoPath = "/path/to/task_video.avi";

    [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = ...
        mapSPDs(videoPath, ...
                "fps", 120, ...
                "windowSpacePixels", [24 24], ...
                "stepSpacePixels", 12, ...
                "windowTimeSecs", 1, ...
                "stepTimeSecs", 0.5, ...
                "doPlot", true);
%}
    arguments
        videoPath {mustBeText}
        options.fps             (1,1) {mustBeNumeric}   = 120
        options.frameIdxToProcess {mustBeNumeric} = [1, inf];
        options.windowSpacePixels          (1,2) {mustBeNumeric}   = [24 24]
        options.stepSpacePixels            (1,1) {mustBeNumeric}   = 12
        options.windowTimeSecs {mustBeNumeric} = 1
        options.stepTimeSecs {mustBeNumeric} = 0.5
        options.doPlot  (1,1) logical           = true
        options.nWorkers (1,1) {mustBeNumeric}   = 6
        options.frameDropVector {mustBeNumeric}  = [];
        options.color_mode {mustBeMember(options.color_mode, ["L+M+S", "L-M", "GRAY", "a"])} = "L+M+S";
        options.verbose = false; 
    end

    % Load in some info about the video to get us started
    video_reader = videoIOWrapper(videoPath, 'ioAction', 'read');
    nRows = video_reader.Height;
    nCols = video_reader.Width;
    nFrames = video_reader.NumFrames;

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

    % Allocate average exponent and AUC variables across all chunks
    exponentMap = nan(nRows, nCols, nChunks);
    varianceMap = nan(nRows, nCols, nChunks);

    % Allocate storage for the spds
    spdByRegion = nan(nRowPatches, nColPatches, nChunks, floor(framesPerChunk/2)-1);

    % Allocate storage for the median image
    medianImage = nan(nRows, nCols, nChunks);

    % Allocate storage for the frq variable
    frq = nan(1,floor(framesPerChunk/2)+1);

    % Define the fmincon options
    optSearch = optimoptions('fmincon');
    optSearch.Display = 'off';

    % Define or extract a frameDropVector
    if ~isempty(options.frameDropVector)
        if length(options.frameDropVector) >= nFrames
            frameDropVector = options.frameDropVector(1:nFrames);
        else
            frameDropVector = zeros(1,nFrames);
            frameDropVector(1:length(options.frameDropVector)) = options.frameDropVector;
        end
    else
        frameDropVector = nan(1,nFrames);
    end


    % Move over the chunks of the video
    for ff = 1:numel(chunkStarts)
        % Turn off a warning that occurs during robust linear fitting
        warnState = warning();
        warning('off','stats:statrobustfit:IterationLimit');

        % Counter for layer index (each patch corresponds to one layer)
        layer = 0;

        % Initialize 3D arrays for slope and AUC values per windowSpacePixels layer
        exponent3D = nan(nRows, nCols, nPatches);
        variance3D = nan(nRows, nCols, nPatches);

        % Inform the user
        startTime = datetime('now');
        fprintf("Processing chunk: %d/%d...", ff, nChunks);

        % Load in a desired amount of frames for this chunk 
        chunk_start_frame = chunkStarts(ff); 
        frameChunk = load_frame_chunk(video_reader, chunk_start_frame, framesPerChunk, options.color_mode, options.verbose); 
        readTime = datetime('now');

        % Nan any frames that we have already specified are dropped
        frameChunkDropVector = frameDropVector(1,chunkStarts(ff):chunkStarts(ff)+framesPerChunk-1);
        frameChunk(frameChunkDropVector==1,:,:) = nan;

        % Write frame drop back into the vector for storage
        for ii = 1:framesPerChunk
            thisFrame = squeeze(frameChunk(ii,:,:));
            frameDropVector(1,chunkStarts(ff)-1+ii) = ...
                all(isnan(thisFrame(:)));
        end

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
                    % If computation fails, skip that patch
                    continue;
                end

                % If spd or frq are nan, skip this patch
                if(any(isnan(spd)) || any(isnan(fLoc)))
                    continue
                else
                    frq = fLoc;
                end

                % Save the raw spd
                spdByRegion(rr,cc,ff,:) = spd(2:end-1);

                % Obtain the slope of the data in the log-log domain, which
                % amounts to the exponent of 1/f^2 model
                p = polyfit(log10(frq(2:end-1)),log10(spd(2:end-1)),1);

                % Derive the total variance of the fitted signal in units of
                % contrast (which is the same as the area under the curve)
                varVal = sum(10.^polyval(p,log10(frq(2:end-1)))) * diff(frq(1:2));

                % Assign exponent and AUC values to current region layer
                variance3D(row:row+windowSpacePixels(1)-1, col:col+windowSpacePixels(2)-1, layer) = varVal;
                exponent3D(row:row+windowSpacePixels(1)-1, col:col+windowSpacePixels(2)-1, layer) = p(1);

            end % col
        end % row

        % Add this to the growing average for the whole video
        exponentMap(:,:,ff) = mean(exponent3D, 3,'omitnan');
        varianceMap(:,:,ff) =  mean(variance3D,3,'omitnan');

        % Finish the console report
        endTime = datetime('now');
        fprintf("read: %2.2fs, process: %2.2fs\n", ...
            seconds(readTime-startTime),seconds(endTime-readTime));

        % Restore the warning state
        warning(warnState);

    end % end chunks

    % Drop the zeroeth and nyquist frequencies
    frq = frq(2:end-1);

    % Obtain the average slope and  AUC map, and the median image
    exponentMap = mean(exponentMap,3,'omitmissing');
    varianceMap = mean(varianceMap,3,'omitmissing');
    medianImage = median(medianImage,3,'omitmissing');

    % Get the mean spd by region
    spdByRegion = squeeze(mean(spdByRegion,3,'omitmissing'));

    % Show the results if requested
    if options.doPlot

        figure('Position',[100,100,900,300]);
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

        nexttile
        imagesc(exponentMap)
        title(sprintf("Exponent Map"));
        axis image;
        colormap('hot');
        colorbar;

        nexttile
        imagesc(varianceMap)
        title(sprintf("Variance Map"));
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

% Define location function to read in a number of frames
% from the video 
function frame_chunk = load_frame_chunk(video_reader, start_frame, num_frames_to_read, color_mode, verbose)
    if(start_frame > video_reader.NumFrames)
        error(sprintf("Start frame %d is out of bounds for video with NumFrames %d", start_frame, video_reader.NumFrames));
    end 
    
    % Determine the true number of frames to read (e.g. on the last chunk
    % there may not be an equal number to read) 
    num_frames_to_read = min([num_frames_to_read, video_reader.NumFrames - start_frame + 1]);

    frame_chunk = zeros(num_frames_to_read, ...
                        video_reader.Height, ...
                        video_reader.Width, ...
                        'double' ...
                        );
    % Read the target amount of frames 
    insertion_index = 1; 
    for ii = start_frame : start_frame + num_frames_to_read - 1
        read_frame = video_reader.readFrame('frameNum', ii, "zeros_as_nans", true, 'color', color_mode, 'verbose', verbose);
        [read_height, read_width] = size(read_frame);
        if(read_height ~= video_reader.Height || read_width ~= video_reader.Width)
            fprintf("Frame shape (%d, %d) does not match video shape (%d, %d)\n", read_height, read_width, video_reader.Height, video_reader.Width); 
            error("FRAME SHAPE HAS BECOME INHOMOGENOUS\n"); 
        end 

        % For SPD analysis, we want to ensure absolutely NO zeros can be in the image
        if (ndims(read_frame) == 2 || size(read_frame,3) == 1)
            % Grayscale: check any zero
            if any(read_frame(:) == 0)
                error("Zero value detected in grayscale frame at frame %d", ii);
            end

        elseif (ndims(read_frame) == 3)
            % Multi-channel: check any pixel that is entirely zero across channels
            zeroMask = true(size(read_frame,1), size(read_frame,2));
            for c = 1:size(read_frame,3)
                zeroMask = zeroMask & (read_frame(:,:,c) == 0);
            end

            if any(zeroMask(:))
                error("Zero pixel (all channels == 0) detected in frame %d", ii);
            end
        end

        frame_chunk(insertion_index, :, :) = read_frame; 
        insertion_index = insertion_index + 1; 
    end 



end 