function outputPaths = makeVirtualFoveationProcessFigures(varargin)
% Make virtual foveation process figures.
%
% Syntax:
%   outputPaths = makeVirtualFoveationProcessFigures(varargin)
%
% Description:
%   This function make virtual foveation process figures.
% Inputs:
%   varargin                 - Input used by the function.
%
% Outputs:
%   outputPaths              - Path-like value returned by the function.
%
% Examples:
%{
    outputPaths = makeVirtualFoveationProcessFigures(varargin)
%}

    makeVirtualFoveationProcessFigures();
    makeVirtualFoveationProcessFigures( ...
        "frameIdx", time2frame([1 31 508], 120), ...
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
        'justProjectionVideo', fullfile(dropBoxBaseDir, 'FLIC_admin', 'Presentations', 'LightLoggerDemoVideos', 'FLIC_2001_walkBiopond_task_justProjection.avi'), ...
        'virtuallyFoveatedVideo', fullfile(dropBoxBaseDir, 'FLIC_admin', 'Presentations', 'LightLoggerDemoVideos', 'FLIC_2001_walkBiopond_task_virtuallyFoveated.avi'), ...
        'frameIdx', [], ...
        'timestamp', [3 48 333], ...
        'fps', 120, ...
        'fovDegrees', 120, ...
        'eccentricitiesDeg', [5 40], ...
        'lineWidth', 3.0, ...
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
% Internal helper to i open video reader.
%
% Syntax:
%   videoReader = iOpenVideoReader(videoPath, readAheadBufferSize)
%
% Description:
%   This local helper function internal helper to i open video reader within its parent workflow.
% Inputs:
%   videoPath                - Path-like input used by the function.
%   readAheadBufferSize      - Input used by the function.
%
% Outputs:
%   videoReader              - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

        videoReader.reader = videoIOWrapper(videoPath, "ioAction", 'read', ...
            "readAheadBufferSize", readAheadBufferSize);
        videoReader.type = "videoIOWrapper";
end

function nFrames = iNumFrames(videoReader)
% Internal helper to i num frames.
%
% Syntax:
%   nFrames = iNumFrames(videoReader)
%
% Description:
%   This local helper function internal helper to i num frames within its parent workflow.
% Inputs:
%   videoReader              - Input used by the function.
%
% Outputs:
%   nFrames                  - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

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
% Internal helper to i read video frame.
%
% Syntax:
%   frame = iReadVideoFrame(videoReader, frameIdx)
%
% Description:
%   This local helper function internal helper to i read video frame within its parent workflow.
% Inputs:
%   videoReader              - Input used by the function.
%   frameIdx                 - Input used by the function.
%
% Outputs:
%   frame                    - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    switch videoReader.type
        case "videoIOWrapper"
            frame = videoReader.reader.readFrame("frameNum", frameIdx, "color", "RGB");
    end
end

function iCloseVideoReader(videoReader)
% Internal helper to i close video reader.
%
% Syntax:
%   iCloseVideoReader(videoReader)
%
% Description:
%   This local helper function internal helper to i close video reader within its parent workflow.
% Inputs:
%   videoReader              - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    if videoReader.type == "videoIOWrapper"
        close(videoReader.reader);
    end
end


function executablePath = iFindExecutable(executableName)
% Internal helper to i find executable.
%
% Syntax:
%   executablePath = iFindExecutable(executableName)
%
% Description:
%   This local helper function internal helper to i find executable within its parent workflow.
% Inputs:
%   executableName           - Input used by the function.
%
% Outputs:
%   executablePath           - Path-like value returned by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

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
% Internal helper to i shell quote.
%
% Syntax:
%   quoted = iShellQuote(value)
%
% Description:
%   This local helper function internal helper to i shell quote within its parent workflow.
% Inputs:
%   value                    - Input used by the function.
%
% Outputs:
%   quoted                   - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    value = char(string(value));
    quoted = "'" + string(strrep(value, "'", "'\''")) + "'";
end

function iDeleteIfExists(filePath)
% Internal helper to i delete if exists.
%
% Syntax:
%   iDeleteIfExists(filePath)
%
% Description:
%   This local helper function internal helper to i delete if exists within its parent workflow.
% Inputs:
%   filePath                 - Path-like input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    if isfile(filePath)
        delete(filePath);
    end
end

function annotatedFrame = iAddEccentricityOverlay(frame, fovDegrees, eccentricitiesDeg, lineWidth, lineStyle)
% Internal helper to i add eccentricity overlay.
%
% Syntax:
%   annotatedFrame = iAddEccentricityOverlay(frame, fovDegrees, eccentricitiesDeg, lineWidth, lineStyle)
%
% Description:
%   This local helper function internal helper to i add eccentricity overlay within its parent workflow.
% Inputs:
%   frame                    - Input used by the function.
%   fovDegrees               - Input used by the function.
%   eccentricitiesDeg        - Input used by the function.
%   lineWidth                - Input used by the function.
%   lineStyle                - Input used by the function.
%
% Outputs:
%   annotatedFrame           - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    annotatedFrame = frame;
    if ~isa(annotatedFrame, "uint8")
        annotatedFrame = uint8(annotatedFrame);
    end
    [nRows, nCols, ~] = size(annotatedFrame);
    centerX = (nCols + 1) / 2;
    centerY = (nRows + 1) / 2;
    degPerPix = fovDegrees / nRows;
    
    % Colors
    black = uint8([0 0 0]);
    grey  = uint8([128 128 128]);
    white = uint8([255 255 255]);
    
    % Outline thickness
    outlineWidth = lineWidth + 2; 
    
    % Compute the arm pixel length corresponding to 10 degrees
    armPix10 = 10 / degPerPix;
    
    %% 1. Draw Meridian Arms (White with Black Outlines)
    % We draw all 4 black outlines first, then all 4 white cores
    
    % --- Black Outlines (Thicker) ---
    annotatedFrame = iDrawLine(annotatedFrame, centerX, 1, centerX, 1 + armPix10, black, outlineWidth, lineStyle); % Top
    annotatedFrame = iDrawLine(annotatedFrame, centerX, nRows, centerX, nRows - armPix10, black, outlineWidth, lineStyle); % Bottom
    annotatedFrame = iDrawLine(annotatedFrame, 1, centerY, 1 + armPix10, centerY, black, outlineWidth, lineStyle); % Left
    annotatedFrame = iDrawLine(annotatedFrame, nCols, centerY, nCols - armPix10, centerY, black, outlineWidth, lineStyle); % Right
    
    % --- White Cores (Original Width) ---
    annotatedFrame = iDrawLine(annotatedFrame, centerX, 1, centerX, 1 + armPix10, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, centerX, nRows, centerX, nRows - armPix10, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, 1, centerY, 1 + armPix10, centerY, white, lineWidth, lineStyle);
    annotatedFrame = iDrawLine(annotatedFrame, nCols, centerY, nCols - armPix10, centerY, white, lineWidth, lineStyle);
    
    %% 2. Draw Eccentricity Rings (Gray with Black Outlines)
    for ee = 1:numel(eccentricitiesDeg)
        radiusPix = eccentricitiesDeg(ee) / degPerPix;
        
        % Black Outline
        annotatedFrame = iDrawCircle(annotatedFrame, centerX, centerY, radiusPix, black, outlineWidth, lineStyle);
        
        % Gray Core
        annotatedFrame = iDrawCircle(annotatedFrame, centerX, centerY, radiusPix, grey, lineWidth, lineStyle);
    end
end

function frame = iDrawCircle(frame, centerX, centerY, radiusPix, color, lineWidth, lineStyle)
% Internal helper to i draw circle.
%
% Syntax:
%   frame = iDrawCircle(frame, centerX, centerY, radiusPix, color, lineWidth, lineStyle)
%
% Description:
%   This local helper function internal helper to i draw circle within its parent workflow.
% Inputs:
%   frame                    - Input used by the function.
%   centerX                  - Input used by the function.
%   centerY                  - Input used by the function.
%   radiusPix                - Input used by the function.
%   color                    - Input used by the function.
%   lineWidth                - Input used by the function.
%   lineStyle                - Input used by the function.
%
% Outputs:
%   frame                    - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

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
% Internal helper to i draw line.
%
% Syntax:
%   frame = iDrawLine(frame, x1, y1, x2, y2, color, lineWidth, lineStyle)
%
% Description:
%   This local helper function internal helper to i draw line within its parent workflow.
% Inputs:
%   frame                    - Input used by the function.
%   x1                       - Input used by the function.
%   y1                       - Input used by the function.
%   x2                       - Input used by the function.
%   y2                       - Input used by the function.
%   color                    - Input used by the function.
%   lineWidth                - Input used by the function.
%   lineStyle                - Input used by the function.
%
% Outputs:
%   frame                    - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

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
% Internal helper to i set disc.
%
% Syntax:
%   frame = iSetDisc(frame, x, y, color, lineWidth)
%
% Description:
%   This local helper function internal helper to i set disc within its parent workflow.
% Inputs:
%   frame                    - Input used by the function.
%   x                        - Input used by the function.
%   y                        - Input used by the function.
%   color                    - Input used by the function.
%   lineWidth                - Input used by the function.
%
% Outputs:
%   frame                    - Output produced by the function.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

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
% Internal helper to i export frame pdf.
%
% Syntax:
%   iExportFramePdf(frame, pdfPath)
%
% Description:
%   This local helper function internal helper to i export frame pdf within its parent workflow.
% Inputs:
%   frame                    - Input used by the function.
%   pdfPath                  - Path-like input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See makeVirtualFoveationProcessFigures.m for usage context.
%}

    fig = figure("Visible", "off", "Color", "w");
    ax = axes("Parent", fig);
    image(ax, frame);
    axis(ax, "image");
    axis(ax, "off");
    set(ax, "Position", [0 0 1 1]);
    exportgraphics(fig, pdfPath, "ContentType", "image");
    close(fig);
end
