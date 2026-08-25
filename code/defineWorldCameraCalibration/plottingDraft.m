function [fig, plotAxes, steps] = plottingDraft(inputImage, varargin)
% plottingDraft  Draft visualization of the world-camera processing stages.
%
%   fig = plottingDraft(inputImage) draws the diagnostic layout used by the
%   Python visualize_pipeline routine for one input image. For now, each
%   processing-stage function is an identity function, so the same image is
%   shown in every Before and After panel.
%
%   [fig, plotAxes, steps] = plottingDraft(inputImage, Name, Value) also
%   returns the 6-by-4 array of main plot axes and the stage image structure.
%
%   Name-value options:
%       FigureTitle            Figure heading.
%       AGCSettings            Frame settings in the order cameraAgain,
%                              AGCDgain, cameraExposure, AGCAgain,
%                              AGCExposure. A 3-value legacy vector is also
%                              accepted.
%       StageTitles            One title per pipeline stage.
%       StageFunctions         One function handle per stage. Each function
%                              accepts the preceding image and returns the
%                              image after that stage. The defaults are
%                              identity functions and are the intended
%                              replacement points for the real pipeline.
%       ZoomBounds             Saturation inset bounds [top bottom left right]
%                              using one-based, inclusive image coordinates.
%       MicroZoomBounds        Pixel-value inset bounds in the same format.
%       FringeZoomBounds       Green-fringe inset bounds in the same format.
%       FringeMicroZoomBounds  Pixel-value fringe bounds in the same format.
%
%   The figure mirrors the Python diagnostic: Before image and statistics,
%   frame AGC settings, before/after R/G/B relative-intensity distributions,
%   and After image. Histograms contain 64 bins over [0, 1], normalize each
%   image by its own nonsaturated maximum, and plot log10(pixel count).

validateInputImage(inputImage);

stageTitles = { ...
    'Linearized Camera Response', ...
    'Fielding Function', ...
    'Color Correction', ...
    'Digital Gain', ...
    'Debayered Image', ...
    'Floor / Ceiling Propagation'};

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'inputImage');
addParameter(parser, 'FigureTitle', 'World Camera Processing Pipeline', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(parser, 'AGCSettings', [1, 1, 8333, 1, 8333], ...
    @(x) isnumeric(x) && isvector(x) && ~isempty(x));
addParameter(parser, 'StageTitles', stageTitles, @isValidStageTitles);
addParameter(parser, 'StageFunctions', defaultStageFunctions(), ...
    @isValidStageFunctions);
addParameter(parser, 'ZoomBounds', [], @isValidBounds);
addParameter(parser, 'MicroZoomBounds', [], @isValidBounds);
addParameter(parser, 'FringeZoomBounds', [], @isValidBounds);
addParameter(parser, 'FringeMicroZoomBounds', [], @isValidBounds);
parse(parser, inputImage, varargin{:});
options = parser.Results;

stageTitles = cellstr(options.StageTitles);
stageFunctions = options.StageFunctions;
if numel(stageTitles) ~= numel(stageFunctions)
    error('plottingDraft:StageCountMismatch', ...
        'StageTitles and StageFunctions must contain the same number of entries.');
end

steps = buildPipelineSteps(inputImage, stageTitles, stageFunctions);
nSteps = numel(steps);

histogramData = repmat(emptyHistogramData(), nSteps, 1);
for ii = 1:nSteps
    histogramData(ii) = getMiddleHistogramData( ...
        steps(ii).beforeImage, steps(ii).afterImage);
end
verifyStageHandoffs(histogramData, steps);

sharedYAxisMax = max([histogramData.yMax]);
[zoomBounds, microZoomBounds, fringeZoomBounds, ...
    fringeMicroZoomBounds, addedSaturationCount] = ...
    resolveDiagnosticBounds(steps, options);

fig = figure( ...
    'Color', 'w', ...
    'Name', char(options.FigureTitle), ...
    'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [40, 40, 1700, 1800]);

annotation(fig, 'textbox', [0.20, 0.965, 0.77, 0.025], ...
    'String', char(options.FigureTitle), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none', ...
    'EdgeColor', 'none');

axisPositions = makeAxisPositions(nSteps);
plotAxes = gobjects(nSteps, 4);
for rowIdx = 1:nSteps
    for colIdx = 1:4
        plotAxes(rowIdx, colIdx) = axes(fig, ...
            'Units', 'normalized', ...
            'Position', axisPositions{rowIdx, colIdx});
    end

    plotImagePanel(fig, plotAxes(rowIdx, 1), ...
        steps(rowIdx).beforeImage, 'Before', '', ...
        zoomBounds, microZoomBounds, fringeZoomBounds, ...
        fringeMicroZoomBounds);

    statsText = formatBeforeAfterStats( ...
        steps(rowIdx).beforeImage, steps(rowIdx).afterImage);
    text(plotAxes(rowIdx, 1), -0.04, 0.5, statsText, ...
        'Units', 'normalized', ...
        'FontSize', 7, ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'right', ...
        'Interpreter', 'tex', ...
        'Clipping', 'off');

    plotAGCSettings(options.AGCSettings, plotAxes(rowIdx, 2));
    plotMiddleHistogram(histogramData(rowIdx), plotAxes(rowIdx, 3), ...
        steps(rowIdx).title, sharedYAxisMax);

    titleExtra = '';
    if rowIdx == nSteps
        titleExtra = sprintf('\n%d new saturated pixels', addedSaturationCount);
    end
    plotImagePanel(fig, plotAxes(rowIdx, 4), ...
        steps(rowIdx).afterImage, 'After', titleExtra, ...
        zoomBounds, microZoomBounds, fringeZoomBounds, ...
        fringeMicroZoomBounds);
end

end


function stageFunctions = defaultStageFunctions()
% Replace these identity handles with the individual pipeline functions.
stageFunctions = { ...
    @(imageData) imageData, ... % Linearize camera response.
    @(imageData) imageData, ... % Apply the fielding correction.
    @(imageData) imageData, ... % Apply Bayer-site color weights.
    @(imageData) imageData, ... % Apply digital gain.
    @(imageData) imageData, ... % De-Bayer to RGB.
    @(imageData) imageData};    % Propagate floor / ceiling contributors.
end


function steps = buildPipelineSteps(inputImage, stageTitles, stageFunctions)
% Run stages sequentially. The defaults leave currentImage unchanged.
nSteps = numel(stageFunctions);
steps = repmat(struct( ...
    'title', '', ...
    'beforeImage', [], ...
    'afterImage', []), nSteps, 1);

currentImage = inputImage;
for ii = 1:nSteps
    beforeImage = currentImage;
    afterImage = stageFunctions{ii}(beforeImage);
    validateInputImage(afterImage);

    steps(ii).title = stageTitles{ii};
    steps(ii).beforeImage = beforeImage;
    steps(ii).afterImage = afterImage;
    currentImage = afterImage;
end
end


function positions = makeAxisPositions(nRows)
% Match Python's gridspec width ratios [1.55, 0.8, 0.8, 1.55].
leftMargin = 0.215;
rightMargin = 0.025;
topMargin = 0.055;
bottomMargin = 0.025;
columnGap = 0.018;
% Leave room for the rotated AGC labels and histogram axis labels between
% rows, as tight_layout does in the Python figure.
rowGap = 0.035;
columnRatios = [1.55, 0.8, 0.8, 1.55];

availableWidth = 1 - leftMargin - rightMargin - 3 * columnGap;
columnWidths = availableWidth * columnRatios / sum(columnRatios);
availableHeight = 1 - topMargin - bottomMargin - (nRows - 1) * rowGap;
rowHeight = availableHeight / nRows;

xPositions = zeros(1, 4);
xPositions(1) = leftMargin;
for colIdx = 2:4
    xPositions(colIdx) = xPositions(colIdx - 1) + ...
        columnWidths(colIdx - 1) + columnGap;
end

positions = cell(nRows, 4);
for rowIdx = 1:nRows
    yPosition = bottomMargin + (nRows - rowIdx) * (rowHeight + rowGap);
    for colIdx = 1:4
        positions{rowIdx, colIdx} = [ ...
            xPositions(colIdx), yPosition, columnWidths(colIdx), rowHeight];
    end
end
end


function plotImagePanel(fig, targetAxis, imageData, imageTitle, titleExtra, ...
        zoomBounds, microZoomBounds, fringeZoomBounds, fringeMicroZoomBounds)
displayImage = imageForDisplay(imageData);
image(targetAxis, displayImage);
axis(targetAxis, 'image');
axis(targetAxis, 'off');
title(targetAxis, sprintf('%s (%s)%s', imageTitle, class(imageData), titleExtra), ...
    'Interpreter', 'none');

addZoomInset(fig, imageData, targetAxis, zoomBounds, microZoomBounds, ...
    [0.01, 0.02, 0.48, 0.50], 'Saturation', [1, 0, 0]);
addZoomInset(fig, imageData, targetAxis, fringeZoomBounds, ...
    fringeMicroZoomBounds, [0.51, 0.02, 0.48, 0.50], ...
    'Fringe', [0, 0.4470, 0.7410]);
end


function addZoomInset(fig, imageData, targetAxis, zoomBounds, ...
        microZoomBounds, insetBounds, zoomTitle, edgeColor)
if isempty(zoomBounds)
    return
end

zoomBounds = clipBounds(zoomBounds, size(imageData));
if isempty(zoomBounds)
    return
end

mainPosition = targetAxis.Position;
insetPosition = [ ...
    mainPosition(1) + insetBounds(1) * mainPosition(3), ...
    mainPosition(2) + insetBounds(2) * mainPosition(4), ...
    insetBounds(3) * mainPosition(3), ...
    insetBounds(4) * mainPosition(4)];
zoomAxis = axes(fig, 'Units', 'normalized', 'Position', insetPosition);

rowStart = zoomBounds(1);
rowEnd = zoomBounds(2);
colStart = zoomBounds(3);
colEnd = zoomBounds(4);
displayImage = imageForDisplay(imageData);
image(zoomAxis, displayImage(rowStart:rowEnd, colStart:colEnd, :));
axis(zoomAxis, 'image');
set(zoomAxis, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box', 'on', ...
    'LineWidth', 1.5, ...
    'XColor', edgeColor, ...
    'YColor', edgeColor);
title(zoomAxis, zoomTitle, 'FontSize', 7, 'Color', edgeColor);

rectangle(targetAxis, ...
    'Position', [colStart - 0.5, rowStart - 0.5, ...
        colEnd - colStart + 1, rowEnd - rowStart + 1], ...
    'EdgeColor', edgeColor, ...
    'FaceColor', 'none', ...
    'LineWidth', 1.5, ...
    'Clipping', 'off');

if isempty(microZoomBounds)
    return
end

microZoomBounds = clipBounds(microZoomBounds, size(imageData));
if isempty(microZoomBounds)
    return
end

microRowStart = microZoomBounds(1);
microRowEnd = microZoomBounds(2);
microColStart = microZoomBounds(3);
microColEnd = microZoomBounds(4);
if microRowStart < rowStart || microRowEnd > rowEnd || ...
        microColStart < colStart || microColEnd > colEnd
    return
end

rectangle(zoomAxis, ...
    'Position', [microColStart - colStart + 0.5, ...
        microRowStart - rowStart + 0.5, ...
        microColEnd - microColStart + 1, ...
        microRowEnd - microRowStart + 1], ...
    'EdgeColor', edgeColor, ...
    'FaceColor', 'none', ...
    'LineWidth', 1.2);

microPosition = [ ...
    insetPosition(1) + 0.22 * insetPosition(3), ...
    insetPosition(2) + 0.03 * insetPosition(4), ...
    0.76 * insetPosition(3), ...
    0.76 * insetPosition(4)];
microAxis = axes(fig, 'Units', 'normalized', 'Position', microPosition);
microDisplay = bayerTintedMicroDisplay(displayImage, microZoomBounds, 0.35);
image(microAxis, microDisplay);
axis(microAxis, 'image');
set(microAxis, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box', 'on', ...
    'LineWidth', 1.2, ...
    'XColor', edgeColor, ...
    'YColor', edgeColor);

microHeight = microRowEnd - microRowStart + 1;
microWidth = microColEnd - microColStart + 1;
title(microAxis, sprintf('%dx%d', microHeight, microWidth), ...
    'FontSize', 9, 'Color', edgeColor);
hold(microAxis, 'on');
for xx = 0.5:1:(microWidth + 0.5)
    plot(microAxis, [xx, xx], [0.5, microHeight + 0.5], ...
        'Color', 'w', 'LineWidth', 0.4);
end
for yy = 0.5:1:(microHeight + 0.5)
    plot(microAxis, [0.5, microWidth + 0.5], [yy, yy], ...
        'Color', 'w', 'LineWidth', 0.4);
end
if mod(microHeight, 2) == 1 && mod(microWidth, 2) == 1
    rectangle(microAxis, ...
        'Position', [floor(microWidth / 2) + 0.5, ...
            floor(microHeight / 2) + 0.5, 1, 1], ...
        'EdgeColor', edgeColor, ...
        'FaceColor', 'none', ...
        'LineWidth', 1.6);
end
addMicroValueLabels(imageData, displayImage, microAxis, microZoomBounds);
end


function displayImage = bayerTintedMicroDisplay(displayImage, bounds, alpha)
rowStart = bounds(1);
rowEnd = bounds(2);
colStart = bounds(3);
colEnd = bounds(4);
displayImage = displayImage(rowStart:rowEnd, colStart:colEnd, :);

for localRow = 1:size(displayImage, 1)
    for localCol = 1:size(displayImage, 2)
        globalRow = rowStart + localRow - 1;
        globalCol = colStart + localCol - 1;
        siteRGB = getBayerSiteRGB(globalRow, globalCol);
        pixelRGB = reshape(displayImage(localRow, localCol, 1:3), 1, 3);
        displayImage(localRow, localCol, 1:3) = ...
            reshape((1 - alpha) * pixelRGB + alpha * siteRGB, 1, 1, 3);
    end
end
displayImage = min(max(displayImage, 0), 1);
end


function siteRGB = getBayerSiteRGB(rowIdx, colIdx)
% The raw camera pattern is BGGR: B at odd/odd and R at even/even.
if mod(rowIdx, 2) == 1 && mod(colIdx, 2) == 1
    siteRGB = [0, 0.25, 1];
elseif mod(rowIdx, 2) == 0 && mod(colIdx, 2) == 0
    siteRGB = [1, 0, 0];
else
    siteRGB = [0, 0.85, 0];
end
end


function addMicroValueLabels(imageData, displayImage, targetAxis, bounds)
rowStart = bounds(1);
rowEnd = bounds(2);
colStart = bounds(3);
colEnd = bounds(4);
microValues = imageData(rowStart:rowEnd, colStart:colEnd, :);
microDisplay = displayImage(rowStart:rowEnd, colStart:colEnd, :);

for rowIdx = 1:size(microValues, 1)
    for colIdx = 1:size(microValues, 2)
        if ndims(microValues) == 3 && size(microValues, 3) >= 3
            pixelValues = reshape(microValues(rowIdx, colIdx, 1:3), 1, 3);
            label = sprintf('R:%s\nG:%s\nB:%s', ...
                formatMicroPixelValue(pixelValues(1)), ...
                formatMicroPixelValue(pixelValues(2)), ...
                formatMicroPixelValue(pixelValues(3)));
        else
            label = formatMicroPixelValue(microValues(rowIdx, colIdx));
        end

        pixelBrightness = mean(microDisplay(rowIdx, colIdx, 1:3), 'all');
        if pixelBrightness > 0.55
            textColor = 'k';
            boxColor = 'w';
        else
            textColor = 'w';
            boxColor = 'k';
        end
        text(targetAxis, colIdx, rowIdx, label, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 7, ...
            'FontWeight', 'bold', ...
            'Color', textColor, ...
            'BackgroundColor', boxColor, ...
            'Margin', 0.2, ...
            'Interpreter', 'none');
    end
end
end


function histogramData = getMiddleHistogramData(beforeImage, afterImage)
rawBeforeChannels = getRGBChannelValues(beforeImage, true);
rawAfterChannels = getRGBChannelValues(afterImage, true);
beforeImageMax = getPositiveFiniteMax(rawBeforeChannels);
afterImageMax = getPositiveFiniteMax(rawAfterChannels);
beforeChannels = normalizeChannelsByImageMax(rawBeforeChannels, beforeImageMax);
afterChannels = normalizeChannelsByImageMax(rawAfterChannels, afterImageMax);

combinedValues = [beforeChannels.R; beforeChannels.G; beforeChannels.B; ...
    afterChannels.R; afterChannels.G; afterChannels.B];
if isempty(combinedValues)
    histogramData = emptyHistogramData();
    return
end

binEdges = linspace(0, 1, 65);
histogramData = emptyHistogramData();
histogramData.isEmpty = false;
histogramData.x = (binEdges(1:end-1) + binEdges(2:end)) / 2;

groupNames = {'Before', 'After'};
groupChannels = {beforeChannels, afterChannels};
channelNames = {'R', 'G', 'B'};
allLogCounts = [];
for groupIdx = 1:2
    for channelIdx = 1:3
        channelName = channelNames{channelIdx};
        values = groupChannels{groupIdx}.(channelName);
        values = min(max(values, 0), 1);
        counts = histcounts(values, binEdges);
        logCounts = nan(size(counts));
        positiveCounts = counts > 0;
        logCounts(positiveCounts) = log10(counts(positiveCounts));
        histogramData.logCounts.(groupNames{groupIdx}).(channelName) = ...
            logCounts;
        allLogCounts = [allLogCounts, logCounts(isfinite(logCounts))]; %#ok<AGROW>
    end
end

if ~isempty(allLogCounts)
    histogramData.yMax = max(allLogCounts);
end
end


function histogramData = emptyHistogramData()
emptyChannels = struct('R', [], 'G', [], 'B', []);
histogramData = struct( ...
    'isEmpty', true, ...
    'x', [], ...
    'logCounts', struct('Before', emptyChannels, 'After', emptyChannels), ...
    'yMax', 0);
end


function plotMiddleHistogram(histogramData, targetAxis, stepTitle, yAxisMax)
if histogramData.isEmpty
    title(targetAxis, stepTitle, 'FontSize', 11, 'FontWeight', 'bold');
    axis(targetAxis, 'off');
    return
end

channelNames = {'R', 'G', 'B'};
channelColors = {[1, 0, 0], [0, 0.5020, 0], [0, 0, 1]};
groupNames = {'Before', 'After'};
plotStyles = {'-x', '-o'};
hold(targetAxis, 'on');
for groupIdx = 1:2
    for channelIdx = 1:3
        channelName = channelNames{channelIdx};
        plot(targetAxis, histogramData.x, ...
            histogramData.logCounts.(groupNames{groupIdx}).(channelName), ...
            plotStyles{groupIdx}, ...
            'LineWidth', 1.2, ...
            'MarkerSize', 8, ...
            'Color', channelColors{channelIdx}, ...
            'DisplayName', sprintf('%s %s', ...
                groupNames{groupIdx}, channelName));
    end
end
title(targetAxis, stepTitle, 'FontSize', 11, 'FontWeight', 'bold');
xlabel(targetAxis, 'Intensity / image max, excluding saturation', 'FontSize', 8);
ylabel(targetAxis, 'log10(pixel count at intensity)', 'FontSize', 8);
xlim(targetAxis, [0, 1]);
if yAxisMax > 0
    ylim(targetAxis, [0, 1.05 * yAxisMax]);
else
    ylim(targetAxis, [0, 1]);
end
set(targetAxis, 'FontSize', 7);
legend(targetAxis, 'Location', 'eastoutside', 'FontSize', 7);
end


function plotAGCSettings(agcSettings, targetAxis)
values = double(agcSettings(:)');
modernNames = {'cameraAgain', 'AGCDgain', 'cameraExposure', ...
    'AGCAgain', 'AGCExposure'};
if numel(values) == numel(modernNames)
    columnNames = modernNames;
elseif numel(values) == 3
    columnNames = {'Again', 'Dgain', 'Exposure'};
else
    columnNames = cell(1, numel(values));
    for ii = 1:numel(values)
        if ii <= numel(modernNames)
            columnNames{ii} = modernNames{ii};
        else
            columnNames{ii} = sprintf('col%d', ii - 1);
        end
    end
end

isExposure = contains(lower(string(columnNames)), 'exposure');
barHeights = values;
positiveExposure = isExposure & values > 0;
barHeights(positiveExposure) = log10(values(positiveExposure));
barHeights(isExposure & ~(values > 0)) = 0;

bars = bar(targetAxis, 1:numel(values), barHeights, 'FaceColor', 'flat');
blue = [0.1216, 0.4667, 0.7059];
orange = [1.0000, 0.4980, 0.0549];
bars.CData = repmat(blue, numel(values), 1);
bars.CData(isExposure, :) = repmat(orange, nnz(isExposure), 1);

title(targetAxis, 'Frame AGC Settings', 'FontSize', 11, 'FontWeight', 'bold');
set(targetAxis, ...
    'XTick', 1:numel(values), ...
    'XTickLabel', columnNames, ...
    'FontSize', 7);
xtickangle(targetAxis, 45);
ylabel(targetAxis, 'value (log10 for exposure)', 'FontSize', 8);

for ii = 1:numel(values)
    text(targetAxis, ii, barHeights(ii) / 2, formatScalar(values(ii)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Rotation', 90, ...
        'FontSize', 7, ...
        'FontWeight', 'bold', ...
        'Color', 'w');
end
end


function channelValues = getRGBChannelValues(imageData, excludeSaturated)
if ismatrix(imageData)
    channelValues.R = imageData(2:2:end, 2:2:end);
    greenOne = imageData(1:2:end, 2:2:end);
    greenTwo = imageData(2:2:end, 1:2:end);
    channelValues.G = [greenOne(:); greenTwo(:)];
    channelValues.B = imageData(1:2:end, 1:2:end);
    channelValues.R = channelValues.R(:);
    channelValues.B = channelValues.B(:);
elseif ndims(imageData) == 3 && size(imageData, 3) >= 3
    redValues = imageData(:, :, 1);
    greenValues = imageData(:, :, 2);
    blueValues = imageData(:, :, 3);
    channelValues = struct( ...
        'R', redValues(:), ...
        'G', greenValues(:), ...
        'B', blueValues(:));
else
    error('plottingDraft:UnsupportedImageShape', ...
        'The image must be a 2-D Bayer image or an RGB image.');
end

channelNames = {'R', 'G', 'B'};
for ii = 1:3
    channelName = channelNames{ii};
    values = channelValues.(channelName);
    validValues = isfinite(values);
    if excludeSaturated
        dtypeMax = getDtypeMaxValue(imageData);
        if ~isempty(dtypeMax)
            validValues = validValues & values ~= dtypeMax;
        end
    end
    channelValues.(channelName) = double(values(validValues));
end
end


function imageMax = getPositiveFiniteMax(channelValues)
values = [channelValues.R; channelValues.G; channelValues.B];
values = values(isfinite(values) & values > 0);
if isempty(values)
    imageMax = 1;
else
    imageMax = max(values);
end
end


function channelValues = normalizeChannelsByImageMax(channelValues, imageMax)
channelNames = {'R', 'G', 'B'};
for ii = 1:3
    channelName = channelNames{ii};
    values = channelValues.(channelName);
    channelValues.(channelName) = values(isfinite(values)) / imageMax;
end
end


function displayImage = imageForDisplay(imageData)
finiteValues = flattenImageValues(imageData, true);
displayData = double(imageData);
if isempty(finiteValues)
    displayBase = zeros(size(imageData, 1), size(imageData, 2));
else
    finiteMin = min(finiteValues);
    finiteMax = max(finiteValues);
    if finiteMax == finiteMin
        finiteMax = finiteMin + 1;
    end

    displayData(isnan(displayData) | isinf(displayData) & displayData < 0) = finiteMin;
    displayData(isinf(displayData) & displayData > 0) = finiteMax;
    if ismatrix(imageData)
        displayBase = (displayData - finiteMin) / (finiteMax - finiteMin);
    else
        displayBase = displayData(:, :, 1:3) / finiteMax;
    end
    displayBase = min(max(displayBase, 0), 1);
end

if ismatrix(displayBase)
    displayImage = repmat(displayBase, 1, 1, 3);
else
    displayImage = displayBase;
end
end


function values = flattenImageValues(imageData, excludeSaturated)
values = imageData(:);
validValues = isfinite(values);
if excludeSaturated
    dtypeMax = getDtypeMaxValue(imageData);
    if ~isempty(dtypeMax)
        validValues = validValues & values ~= dtypeMax;
    end
end
values = double(values(validValues));
end


function statsText = formatBeforeAfterStats(beforeImage, afterImage)
beforeStats = getImageStats(beforeImage);
afterStats = getImageStats(afterImage);
statsText = sprintf([ ...
    '\\bfValues (before / after)\\rm\n' ...
    'observed min, excluding saturation: %s / %s\n' ...
    'observed max, excluding saturation: %s / %s\n' ...
    '\\bfPixel Counts (before / after)\\rm\n' ...
    'zero value: %d / %d\n' ...
    'at observed max: %d / %d\n' ...
    'Inf: %d / %d\n' ...
    'at dtype max: %d / %d\n' ...
    'NaN: %d / %d'], ...
    formatScalar(beforeStats.minExcludingSaturated), ...
    formatScalar(afterStats.minExcludingSaturated), ...
    formatScalar(beforeStats.maxExcludingSaturated), ...
    formatScalar(afterStats.maxExcludingSaturated), ...
    beforeStats.zeroCount, afterStats.zeroCount, ...
    beforeStats.observedMaxCount, afterStats.observedMaxCount, ...
    beforeStats.infCount, afterStats.infCount, ...
    beforeStats.dtypeMaxCount, afterStats.dtypeMaxCount, ...
    beforeStats.nanCount, afterStats.nanCount);
end


function stats = getImageStats(imageData)
nonsaturatedValues = flattenImageValues(imageData, true);
if isempty(nonsaturatedValues)
    observedMin = [];
    observedMax = [];
    observedMaxCount = 0;
else
    observedMin = min(nonsaturatedValues);
    observedMax = max(nonsaturatedValues);
    observedMaxCount = nnz(pixelMaskFromValueMask(imageData == observedMax));
end

dtypeMax = getDtypeMaxValue(imageData);
if isempty(dtypeMax)
    dtypeMaxCount = 0;
else
    dtypeMaxCount = nnz(pixelMaskFromValueMask(imageData == dtypeMax));
end

stats = struct( ...
    'minExcludingSaturated', observedMin, ...
    'maxExcludingSaturated', observedMax, ...
    'observedMaxCount', observedMaxCount, ...
    'zeroCount', nnz(pixelMaskFromValueMask(imageData == 0)), ...
    'infCount', nnz(pixelMaskFromValueMask(isinf(imageData))), ...
    'dtypeMaxCount', dtypeMaxCount, ...
    'nanCount', nnz(pixelMaskFromValueMask(isnan(imageData))));
end


function dtypeMax = getDtypeMaxValue(imageData)
if isinteger(imageData)
    dtypeMax = double(intmax(class(imageData)));
elseif isfloat(imageData)
    dtypeMax = realmax(class(imageData));
elseif islogical(imageData)
    dtypeMax = 1;
else
    dtypeMax = [];
end
end


function pixelMask = pixelMaskFromValueMask(valueMask)
if ndims(valueMask) == 3
    pixelMask = any(valueMask(:, :, 1:min(3, size(valueMask, 3))), 3);
else
    pixelMask = valueMask;
end
end


function [zoomBounds, microZoomBounds, fringeZoomBounds, ...
        fringeMicroZoomBounds, addedSaturationCount] = ...
        resolveDiagnosticBounds(steps, options)
finalBefore = steps(end).beforeImage;
finalAfter = steps(end).afterImage;

changedMask = getChangedPixelMask(finalBefore, finalAfter);
addedSaturationMask = getSaturatedPixelMask(finalAfter) & ...
    ~getSaturatedPixelMask(finalBefore);
addedSaturationCount = nnz(addedSaturationMask);
if any(addedSaturationMask, 'all')
    zoomMask = addedSaturationMask;
else
    zoomMask = changedMask;
end

zoomBounds = options.ZoomBounds;
if isempty(zoomBounds)
    zoomBounds = getZoomBounds(zoomMask, size(finalAfter), 32);
end

microZoomBounds = options.MicroZoomBounds;
if isempty(microZoomBounds)
    microCenter = getGreenPixelCenter(finalBefore, zoomBounds);
    if isempty(microCenter)
        microCenter = getGreenPixelCenter(finalAfter, zoomBounds);
    end
    microZoomBounds = getPixelCropBounds(microCenter, size(finalAfter), 3);
end

fringeCenter = getGreenFringePixelCenter(finalBefore, finalAfter);
fringeZoomBounds = options.FringeZoomBounds;
if isempty(fringeZoomBounds)
    fringeZoomBounds = getZoomBoundsAroundPixel( ...
        fringeCenter, size(finalAfter), 32);
end

fringeMicroZoomBounds = options.FringeMicroZoomBounds;
if isempty(fringeMicroZoomBounds)
    fringeMicroZoomBounds = getPixelCropBounds( ...
        fringeCenter, size(finalAfter), 3);
end
end


function changedMask = getChangedPixelMask(beforeImage, afterImage)
if ~isequal(size(beforeImage), size(afterImage))
    changedMask = false(size(afterImage, 1), size(afterImage, 2));
    return
end
matchingValues = beforeImage == afterImage | ...
    isnan(beforeImage) & isnan(afterImage);
changedMask = pixelMaskFromValueMask(~matchingValues);
end


function saturatedMask = getSaturatedPixelMask(imageData)
infMask = pixelMaskFromValueMask(isinf(imageData));
dtypeMax = getDtypeMaxValue(imageData);
if isempty(dtypeMax)
    dtypeMaxMask = false(size(imageData, 1), size(imageData, 2));
else
    dtypeMaxMask = pixelMaskFromValueMask(imageData == dtypeMax);
end
saturatedMask = infMask | dtypeMaxMask;
end


function bounds = getZoomBounds(pixelMask, imageSize, cropSize)
[rows, cols] = find(pixelMask);
if isempty(rows)
    bounds = [];
    return
end

height = imageSize(1);
width = imageSize(2);
cropHeight = min(cropSize, height);
cropWidth = min(cropSize, width);
medianRow = median(rows);
medianCol = median(cols);
[~, centerIdx] = min((rows - medianRow).^2 + (cols - medianCol).^2);
centerRow = rows(centerIdx);
centerCol = cols(centerIdx);
rowStart = max(1, min(centerRow - floor(cropHeight / 2), ...
    height - cropHeight + 1));
colStart = max(1, min(centerCol - floor(cropWidth / 2), ...
    width - cropWidth + 1));
bounds = [rowStart, rowStart + cropHeight - 1, ...
    colStart, colStart + cropWidth - 1];
end


function bounds = getPixelCropBounds(centerPixel, imageSize, cropSize)
if isempty(centerPixel)
    bounds = [];
    return
end

height = imageSize(1);
width = imageSize(2);
cropSize = min([cropSize, height, width]);
if cropSize > 1 && mod(cropSize, 2) == 0
    cropSize = cropSize - 1;
end
rowStart = max(1, min(centerPixel(1) - floor(cropSize / 2), ...
    height - cropSize + 1));
colStart = max(1, min(centerPixel(2) - floor(cropSize / 2), ...
    width - cropSize + 1));
bounds = [rowStart, rowStart + cropSize - 1, ...
    colStart, colStart + cropSize - 1];
end


function centerPixel = getGreenPixelCenter(imageData, zoomBounds)
if isempty(zoomBounds)
    centerPixel = [];
    return
end

zoomBounds = clipBounds(zoomBounds, size(imageData));
rowStart = zoomBounds(1);
rowEnd = zoomBounds(2);
colStart = zoomBounds(3);
colEnd = zoomBounds(4);
displayCrop = imageForDisplay(imageData);
displayCrop = displayCrop(rowStart:rowEnd, colStart:colEnd, 1:3);
if isempty(displayCrop)
    centerPixel = [];
    return
end

greenScore = displayCrop(:, :, 2) - ...
    max(displayCrop(:, :, 1), displayCrop(:, :, 3));
brightness = max(displayCrop, [], 3);
candidateMask = greenScore > 0.05 & brightness > 0.05;
if size(candidateMask, 1) >= 3 && size(candidateMask, 2) >= 3
    interiorMask = candidateMask;
    interiorMask([1, end], :) = false;
    interiorMask(:, [1, end]) = false;
    if any(interiorMask, 'all')
        candidateMask = interiorMask;
    end
end

if ~any(candidateMask, 'all')
    maximumGreenScore = max(greenScore, [], 'all');
    if maximumGreenScore <= 0
        centerPixel = [];
        return
    end
    candidateMask = greenScore == maximumGreenScore;
end

candidateScores = greenScore + 0.05 * brightness;
candidateScores(~candidateMask) = -Inf;
[~, linearIdx] = max(candidateScores, [], 'all', 'linear');
[cropRow, cropCol] = ind2sub(size(candidateScores), linearIdx);
centerPixel = [rowStart + cropRow - 1, colStart + cropCol - 1];
end


function centerPixel = getGreenFringePixelCenter(beforeImage, afterImage)
if ~isequal(size(beforeImage, 1), size(afterImage, 1)) || ...
        ~isequal(size(beforeImage, 2), size(afterImage, 2))
    centerPixel = [];
    return
end

afterDisplay = imageForDisplay(afterImage);
beforeDisplay = imageForDisplay(beforeImage);
afterGreenScore = afterDisplay(:, :, 2) - ...
    max(afterDisplay(:, :, 1), afterDisplay(:, :, 3));
beforeGreenScore = beforeDisplay(:, :, 2) - ...
    max(beforeDisplay(:, :, 1), beforeDisplay(:, :, 3));
brightness = max(afterDisplay, [], 3);
finiteAfterMask = ~getSaturatedPixelMask(afterImage);
candidateMask = finiteAfterMask & afterGreenScore > 0.05 & brightness > 0.05;

if ~any(candidateMask, 'all')
    candidateMask = finiteAfterMask & beforeGreenScore > 0.05;
end

if ~any(candidateMask, 'all')
    candidateScores = max(afterGreenScore, beforeGreenScore);
    candidateScores(~finiteAfterMask) = -Inf;
    if ~any(isfinite(candidateScores), 'all') || ...
            max(candidateScores, [], 'all') <= 0
        centerPixel = [];
        return
    end
else
    candidateScores = afterGreenScore + 0.25 * beforeGreenScore + ...
        0.05 * brightness;
    candidateScores(~candidateMask) = -Inf;
end

[~, linearIdx] = max(candidateScores, [], 'all', 'linear');
[centerRow, centerCol] = ind2sub(size(candidateScores), linearIdx);
centerPixel = [centerRow, centerCol];
end


function bounds = getZoomBoundsAroundPixel(centerPixel, imageSize, cropSize)
if isempty(centerPixel)
    bounds = [];
    return
end
pixelMask = false(imageSize(1), imageSize(2));
pixelMask(centerPixel(1), centerPixel(2)) = true;
bounds = getZoomBounds(pixelMask, imageSize, cropSize);
end


function bounds = clipBounds(bounds, imageSize)
if isempty(bounds)
    return
end
bounds = round(double(bounds(:)'));
bounds(1) = max(1, min(bounds(1), imageSize(1)));
bounds(2) = max(1, min(bounds(2), imageSize(1)));
bounds(3) = max(1, min(bounds(3), imageSize(2)));
bounds(4) = max(1, min(bounds(4), imageSize(2)));
if bounds(1) > bounds(2) || bounds(3) > bounds(4)
    bounds = [];
end
end


function verifyStageHandoffs(histogramData, steps)
for ii = 1:(numel(histogramData) - 1)
    if histogramData(ii).isEmpty || histogramData(ii + 1).isEmpty
        continue
    end
    channelNames = {'R', 'G', 'B'};
    for channelIdx = 1:3
        channelName = channelNames{channelIdx};
        previousAfter = histogramData(ii).logCounts.After.(channelName);
        nextBefore = histogramData(ii + 1).logCounts.Before.(channelName);
        assert(isequaln(previousAfter, nextBefore), ...
            ['Middle plot mismatch: After line for "%s" does not match ' ...
            'Before line for "%s" in %s.'], ...
            steps(ii).title, steps(ii + 1).title, channelName);
    end
end
end


function label = formatMicroPixelValue(value)
value = double(value);
if isnan(value)
    label = 'NaN';
elseif isinf(value) && value > 0
    label = 'Inf';
elseif isinf(value)
    label = '-Inf';
elseif abs(value - round(value)) <= eps(max(1, abs(value)))
    label = sprintf('%d', round(value));
else
    label = sprintf('%.3g', value);
end
end


function label = formatScalar(value)
if isempty(value)
    label = 'n/a';
elseif isfinite(value) && abs(value - round(value)) <= eps(max(1, abs(value)))
    label = sprintf('%d', round(value));
else
    label = sprintf('%.3g', value);
end
end


function validateInputImage(imageData)
if ~(isnumeric(imageData) || islogical(imageData)) || ~isreal(imageData) || ...
        isempty(imageData)
    error('plottingDraft:InvalidImage', ...
        'inputImage must be a nonempty, real numeric or logical image.');
end
if ~(ismatrix(imageData) || (ndims(imageData) == 3 && size(imageData, 3) >= 3))
    error('plottingDraft:InvalidImageShape', ...
        'inputImage must be a 2-D Bayer image or have at least three channels.');
end
end


function tf = isValidStageTitles(value)
tf = (isstring(value) && isvector(value) && ~isempty(value)) || ...
    (iscell(value) && ~isempty(value) && ...
    all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), value)));
end


function tf = isValidStageFunctions(value)
tf = iscell(value) && ~isempty(value) && ...
    all(cellfun(@(x) isa(x, 'function_handle'), value));
end


function tf = isValidBounds(value)
tf = isempty(value) || (isnumeric(value) && isreal(value) && ...
    numel(value) == 4 && all(isfinite(value)));
end
