function outputPaths = makeVirtualFoveationProcessFigures(varargin)
% Make annotated frames illustrating virtual foveation/projection videos.
%
% Examples:
%{
    makeVirtualFoveationProcessFigures();
    makeVirtualFoveationProcessFigures( ...
        "frameIdx", time2frame([3 48 333], 120), ...
        "lineStyle", "dashed");
%}
    % Parse name-value pairs or single structure input
    options = struct();
    if numel(varargin) == 1 && isstruct(varargin{1})
        options = varargin{1};
    elseif mod(numel(varargin), 2) == 0
        for ii = 1:2:numel(varargin)
            options.(varargin{ii}) = varargin{ii+1};
        end
    end

    % Retrieve dropbox base directory dynamically with fallback
    try
        dropBoxBaseDir = getpref('combiExperiments', 'dropboxBaseDir');
    catch
        dropBoxBaseDir = '';
    end
    
    if isempty(dropBoxBaseDir)
        dropBoxBaseDir = '/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya';
    end
    
    % Default parameters with updated line style option
    defaultOptions = struct( ...
        'justProjectionVideo', fullfile(dropBoxBaseDir, 'FLIC_admin', 'Presentations', 'LightLoggerDemoVideos', 'FLIC_2001_walkIndoor_task_justProjection.avi'), ...
        'virtuallyFoveatedVideo', fullfile(dropBoxBaseDir, 'FLIC_admin', 'Presentations', 'LightLoggerDemoVideos', 'FLIC_2001_walkIndoor_task_virtuallyFoveated.avi'), ...
        'frameIdx', [], ...
        'timestamp', [3 48 333], ...
        'fps', 120, ...
        'fovDegrees', 120, ...
        'eccentricitiesDeg', [5 40], ...
        'lineWidth', 2.0, ...
        'lineStyle', 'solid', ... % Options: 'solid', 'dashed'
        'outputDir', fullfile(dropBoxBaseDir, 'FLIC_admin', 'Presentations', 'Light Logger Figures', 'virtual foveation figures'), ...
        'writeAnnotatedVideo', false ...
    );
    
    % Populate any missing fields with defaults
    fn = fieldnames(defaultOptions);
    for i = 1:numel(fn)
        if ~isfield(options, fn{i})
            options.(fn{i}) = defaultOptions.(fn{i});
        end
    end
    
    if isempty(options.frameIdx)
        frameIdx = round(time2frame(options.timestamp, options.fps));
    else
        frameIdx = round(options.frameIdx(:)');
    end
    
    if ~exist(options.outputDir, "dir")
        mkdir(options.outputDir);
    end
    
    videoSpecs = struct( ...
        "label", {"justProjection", "virtuallyFoveated"}, ...
        "path", {string(options.justProjectionVideo), string(options.virtuallyFoveatedVideo)} ...
    );
    
    outputPaths = struct();
    for vv = 1:numel(videoSpecs)
        videoLabel = videoSpecs(vv).label;
        videoPath = videoSpecs(vv).path;
        assert(isfile(videoPath), "Could not find video: %s", videoPath);
        
        videoReader = iOpenVideoReader(videoPath, max(1000, numel(frameIdx)));
        nVideoFrames = iNumFrames(videoReader);
        frameIdxThisVideo = frameIdx(frameIdx >= 1 & frameIdx <= nVideoFrames);
        
        assert(~isempty(frameIdxThisVideo), ...
            "Requested frame(s) are outside %s frame range 1:%d.", videoPath, nVideoFrames);
        
        outputPaths.(videoLabel).png = strings(1, numel(frameIdxThisVideo));
        outputPaths.(videoLabel).pdf = strings(1, numel(frameIdxThisVideo));
        
        if options.writeAnnotatedVideo
            annotatedVideoPath = fullfile(options.outputDir, sprintf("%s_annotatedFrames.avi", videoLabel));
            videoWriter = VideoWriter(annotatedVideoPath);
            videoWriter.FrameRate = options.fps;
            open(videoWriter);
            outputPaths.(videoLabel).video = string(annotatedVideoPath);
        end
        
        for ff = 1:numel(frameIdxThisVideo)
            thisFrameIdx = frameIdxThisVideo(ff);
            frame = iReadVideoFrame(videoReader, thisFrameIdx);
            annotatedFrame = iAddEccentricityOverlay(frame, ...
                options.fovDegrees, options.eccentricitiesDeg, options.lineWidth, options.lineStyle);
            
            % Pull up the image in a figure as it is created
            figName = 'makeVirtualFoveationProcessFigures_Display';
            figHandle = findobj('Type', 'Figure', 'Name', figName);
            if isempty(figHandle)
                figHandle = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
            else
                figure(figHandle);
            end
            
            imshow(annotatedFrame);
            title(sprintf('%s - Frame %d', videoLabel, thisFrameIdx), 'Interpreter', 'none');
            drawnow;
            
            pngPath = fullfile(options.outputDir, ...
                sprintf("%s_frame%06d_overlay.png", videoLabel, thisFrameIdx));
            pdfPath = fullfile(options.outputDir, ...
                sprintf("%s_frame%06d_overlay.pdf", videoLabel, thisFrameIdx));
            
            imwrite(annotatedFrame, pngPath);
            iExportFramePdf(annotatedFrame, pdfPath);
            
            outputPaths.(videoLabel).png(ff) = string(pngPath);
            outputPaths.(videoLabel).pdf(ff) = string(pdfPath);
            
            if options.writeAnnotatedVideo
                writeVideo(videoWriter, annotatedFrame);
            end
        end
        iCloseVideoReader(videoReader);
        if options.writeAnnotatedVideo
            close(videoWriter);
        end
    end
end

function videoReader = iOpenVideoReader(videoPath, readAheadBufferSize)
    try
        videoReader.reader = videoIOWrapper(videoPath, "ioAction", 'read', ...
            "readAheadBufferSize", readAheadBufferSize);
        videoReader.type = "videoIOWrapper";
    catch err
        warning("makeVirtualFoveationProcessFigures:VideoIOWrapperFallback", ...
            "Falling back to MATLAB VideoReader for %s. videoIOWrapper error: %s", ...
            videoPath, err.message);
        try
            videoReader.reader = VideoReader(videoPath);
            videoReader.type = "VideoReader";
        catch videoReaderErr
            warning("makeVirtualFoveationProcessFigures:FfmpegFallback", ...
                "Falling back to ffmpeg frame extraction for %s. VideoReader error: %s", ...
                videoPath, videoReaderErr.message);
            videoReader.reader = string(videoPath);
            videoReader.type = "ffmpeg";
        end
    end
end

function nFrames = iNumFrames(videoReader)
    switch videoReader.type
        case "videoIOWrapper"
            nFrames = videoReader.reader.NumFrames;
        case "VideoReader"
            nFrames = floor(videoReader.reader.Duration * videoReader.reader.FrameRate);
        case "ffmpeg"
            nFrames = inf;
    end
end

function frame = iReadVideoFrame(videoReader, frameIdx)
    switch videoReader.type
        case "videoIOWrapper"
            frame = videoReader.reader.readFrame("frameNum", frameIdx, "color", "RGB");
        case "VideoReader"
            frame = read(videoReader.reader, frameIdx);
        case "ffmpeg"
            frame = iReadFrameWithFfmpeg(videoReader.reader, frameIdx);
    end
end

function iCloseVideoReader(videoReader)
    if videoReader.type == "videoIOWrapper"
        close(videoReader.reader);
    end
end

function frame = iReadFrameWithFfmpeg(videoPath, frameIdx)
    ffmpegPath = iFindExecutable("ffmpeg");
    tempPng = string(tempname) + ".png";
    cleanupObj = onCleanup(@() iDeleteIfExists(tempPng));
    zeroBasedFrameIdx = frameIdx - 1;
    cmd = sprintf("%s -v error -y -i %s -vf %s -vframes 1 %s", ...
        iShellQuote(ffmpegPath), ...
        iShellQuote(videoPath), ...
        iShellQuote(sprintf("select=eq(n\\,%d)", zeroBasedFrameIdx)), ...
        iShellQuote(tempPng));
    [status, msg] = system(cmd);
    if status ~= 0 || ~isfile(tempPng)
        error("Could not extract frame %d from %s with ffmpeg. %s", ...
            frameIdx, videoPath, msg);
    end
    frame = imread(tempPng);
end

function executablePath = iFindExecutable(executableName)
    homebrewPath = "/opt/homebrew/bin/" + executableName;
    if isfile(homebrewPath)
        executablePath = homebrewPath;
        return;
    end
    [status, pathText] = system("command -v " + executableName);
    if status ~= 0
        error("Could not find required executable: %s", executableName);
    end
    executablePath = strtrim(string(pathText));
end

function quoted = iShellQuote(value)
    value = char(string(value));
    quoted = "'" + string(strrep(value, "'", "'\''")) + "'";
end

function iDeleteIfExists(filePath)
    if isfile(filePath)
        delete(filePath);
    end
end

function annotatedFrame = iAddEccentricityOverlay(frame, fovDegrees, eccentricitiesDeg, lineWidth, lineStyle)
    annotatedFrame = frame;
    if ~isa(annotatedFrame, "uint8")
        annotatedFrame = uint8(annotatedFrame);
    end
    [nRows, nCols, ~] = size(annotatedFrame);
    centerX = (nCols + 1) / 2;
    centerY = (nRows + 1) / 2;
    degPerPix = fovDegrees / nRows;
    
    white = uint8([255 255 255]);
    grey = uint8([128 128 128]);
    
    % Compute the arm pixel length corresponding to 10 degrees
    armPix10 = 10 / degPerPix;
    
    % Draw only the 10-degree arms of the meridian lines at the 4 outer ends
    annotatedFrame = iDrawLine(annotatedFrame, centerX, 1, centerX, 1 + armPix10, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, centerX, nRows, centerX, nRows - armPix10, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, 1, centerY, 1 + armPix10, centerY, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, nCols, centerY, nCols - armPix10, centerY, white, lineWidth, lineStyle);
    
    % Dynamically draw rings using the eccentricitiesDeg array
    for ee = 1:numel(eccentricitiesDeg)
        radiusPix = eccentricitiesDeg(ee) / degPerPix;
        % All eccentricity rings are grey
        annotatedFrame = iDrawCircle(annotatedFrame, centerX, centerY, radiusPix, grey, lineWidth, lineStyle);
    end
end

function frame = iDrawCircle(frame, centerX, centerY, radiusPix, color, lineWidth, lineStyle)
    circumference = 2 * pi * radiusPix;
    nPoints = ceil(circumference);
    nPoints = max(nPoints, 36); % Ensure at least 36 points
    theta = linspace(0, 2*pi, nPoints);
    x = centerX + radiusPix * cos(theta);
    y = centerY + radiusPix * sin(theta);
    
    for pp = 1:numel(x)
        switch lower(lineStyle)
            case 'solid'
                frame = iSetDisc(frame, x(pp), y(pp), color, lineWidth);
            case 'dashed'
                % Half length and closer together (period 18, length 6)
                if mod(pp - 1, 18) < 6
                    frame = iSetDisc(frame, x(pp), y(pp), color, lineWidth);
                end
        end
    end
end

function frame = iDrawLine(frame, x1, y1, x2, y2, color, lineWidth, lineStyle)
    nPoints = ceil(max(abs(x2 - x1), abs(y2 - y1))) + 1;
    x = linspace(x1, x2, nPoints);
    y = linspace(y1, y2, nPoints);
    
    for pp = 1:nPoints
        switch lower(lineStyle)
            case 'solid'
                frame = iSetDisc(frame, x(pp), y(pp), color, lineWidth);
            case 'dashed'
                % Half length and closer together (period 18, length 6)
                if mod(pp - 1, 18) < 6
                    frame = iSetDisc(frame, x(pp), y(pp), color, lineWidth);
                end
        end
    end
end

function frame = iSetDisc(frame, x, y, color, lineWidth)
    [nRows, nCols, ~] = size(frame);
    radius = max(0, floor((lineWidth - 1) / 2));
    colCenter = round(x);
    rowCenter = round(y);
    rows = max(1, rowCenter-radius):min(nRows, rowCenter+radius);
    cols = max(1, colCenter-radius):min(nCols, colCenter+radius);
    for cc = 1:3
        channel = frame(:,:,cc);
        channel(rows, cols) = color(cc);
        frame(:,:,cc) = channel;
    end
end

function iExportFramePdf(frame, pdfPath)
    fig = figure("Visible", "off", "Color", "w");
    ax = axes("Parent", fig);
    image(ax, frame);
    axis(ax, "image");
    axis(ax, "off");
    set(ax, "Position", [0 0 1 1]);
    exportgraphics(fig, pdfPath, "ContentType", "image");
    close(fig);
end