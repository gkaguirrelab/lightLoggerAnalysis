function [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityData, options)
% Plot temporal SPD exponent, variance, and regional spectrum summaries
%
% Syntax:
%   [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityData)
%   [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityData, options)
%
% Description:
%   This function visualizes temporal SPD analysis results.
%
%   Supported input cases:
%
%   CASE 1:
%       activityData is a string/char path to a .mat file containing either
%       a variable named activityData or a single struct variable.
%
%   CASE 2:
%       activityData is a struct for a single activity where:
%           activityData.(activityName).exponentMap
%           activityData.(activityName).varianceMap
%           activityData.(activityName).spdByRegion
%           activityData.(activityName).frq
%           activityData.(activityName).medianImage
%           activityData.(activityName).frameDropVector
%
%   CASE 3:
%       activityData.(activityName).justProjection and/or
%       activityData.(activityName).virtuallyFoveated exist, and each one is
%       either:
%           (a) a struct with exponentMap / varianceMap / spdByRegion / frq
%       or  (b) a string/char path to a .mat file containing such a struct.
%
%       In this case:
%           - exponent maps are shown side-by-side in one figure
%           - variance maps are shown side-by-side in one figure
%           - SPD curves from both projection types are overlaid on one plot
%
% Optional key/value pairs:
%   fovDegrees         - Scalar FOV in degrees used to convert eccentricity
%                        ring radii into pixels.
%   exponent_clim      - false or 1x2 numeric for exponent map color limits.
%   variance_clim      - false or 1x2 numeric for variance map color limits.
%   spd_xlim           - false or 1x2 numeric for SPD x-limits.
%   spd_ylim           - false or 1x2 numeric for SPD y-limits.
%   output_dir         - false or string/char output directory for exports.
%   overwrite_existing - logical, whether to overwrite existing exports.
%   title              - title stem used in exported filenames.
%
% Outputs:
%   exponentMapHandle  - Figure handle for exponent figure.
%   varianceMapHandle  - Figure handle for variance figure.
%   spdByRegionHandle  - Figure handle for SPD-by-region figure.

    arguments
        activityData
        options.fovDegrees = 120
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.output_dir = false
        options.overwrite_existing = false
        options.title = ""
    end

    %% Case 1 support: if the top-level input is a path, load it
    if (ischar(activityData) || isstring(activityData))
        activityData = iLoadActivityDataStruct(activityData);
    end

    %% Determine the activity name
    activityNames = fieldnames(activityData);
    assert(~isempty(activityNames), 'activityData has no fields.');
    activityName = activityNames{1};

    % Main activity container
    activityStruct = activityData.(activityName);

    %% Determine whether this is case 2 or case 3
    isSingleProjectionCase = isfield(activityStruct, 'exponentMap');

    projectionNamesPreferredOrder = {'justProjection', 'virtuallyFoveated'};

    if (isSingleProjectionCase)
        projectionNames = {'singleProjection'};
    else
        projectionNames = {};
        for ii = 1:numel(projectionNamesPreferredOrder)
            if (isfield(activityStruct, projectionNamesPreferredOrder{ii}))
                projectionNames{end+1} = projectionNamesPreferredOrder{ii}; %#ok<AGROW>
            end
        end

        assert(~isempty(projectionNames), ...
            ['Could not detect supported input structure. Expected either direct data fields ' ...
             '(case 2) or projection-type fields such as justProjection / virtuallyFoveated (case 3).']);
    end

    %% Standardize all projection data into one struct
    % After this block, plotData.(projectionName) will always directly contain:
    %   exponentMap, varianceMap, spdByRegion, frq, ...
    plotData = struct();

    if (isSingleProjectionCase)
        plotData.singleProjection = iEnsureProjectionDataStruct(activityStruct);
    else
        for ii = 1:numel(projectionNames)
            projectionName = projectionNames{ii};
            plotData.(projectionName) = iEnsureProjectionDataStruct(activityStruct.(projectionName));
        end
    end

    %% Build FOV mask from the first available projection
    firstProjectionNames = fieldnames(plotData);
    firstProjectionName = firstProjectionNames{1};

    exampleExponentMap = plotData.(firstProjectionName).exponentMap;

    ellipseTransparentParams = [240, 240, 120000, .75, 0];
    p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
    myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
    [X, Y] = meshgrid(1:480, 1:480);
    mask = double(myEllipse(X,Y) < 1e-9);

    degPerPix = options.fovDegrees / size(exampleExponentMap, 1);
    idxStarts = 1:12:(480 - 24 + 1);

    %% CASE 2: single projection
    if (isSingleProjectionCase)

        exponentMap = plotData.singleProjection.exponentMap;
        varianceMap = plotData.singleProjection.varianceMap;
        spdByRegion = plotData.singleProjection.spdByRegion;
        frq = plotData.singleProjection.frq;

        exponentMap(mask == 0) = nan;
        varianceMap(mask == 0) = nan;

        % Exponent map
        exponentMapHandle = figure;
        imagesc(-exponentMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Exponent - %s', activityName))
        axis square
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar

        if (~islogical(options.output_dir))
            assert(~islogical(options.title));
            output_filepath = fullfile(options.output_dir, sprintf('%s_exponentMap.pdf', options.title));
            if (~exist(output_filepath, 'file') || options.overwrite_existing)
                exportgraphics(exponentMapHandle, output_filepath, 'ContentType', 'vector');
            end
            close(exponentMapHandle);
        end

        % Variance map
        varianceMapHandle = figure;
        imagesc(varianceMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Contrast Variance - %s', activityName))
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar

        if (~islogical(options.output_dir))
            assert(~islogical(options.title));
            output_filepath = fullfile(options.output_dir, sprintf('%s_varianceMap.pdf', options.title));
            if (~exist(output_filepath, 'file') || options.overwrite_existing)
                exportgraphics(varianceMapHandle, output_filepath, 'ContentType', 'vector');
            end
            close(varianceMapHandle);
        end

        % SPD by region
        spdByRegionHandle = figure;
        loglog(frq, squeeze(spdByRegion(20,20,:)), '-k', 'LineWidth', 1.5);
        hold on
        loglog(frq, squeeze(spdByRegion(31,20,:)), '-r', 'LineWidth', 1.5);
        plot([10^0 10^1.5], [10^-2 10^-5], ':k')
        legend({'center','periphery'}, 'Location', 'best');
        ylabel('Power [contrast^2/Hz]');
        xlabel('Frequency [log Hz]');
        title(sprintf('SPDs from the center and periphery - %s', activityName));

        if (~islogical(options.spd_xlim))
            xlim(double(options.spd_xlim));
        end
        if (~islogical(options.spd_ylim))
            ylim(double(options.spd_ylim));
        end

        if (~islogical(options.output_dir))
            assert(~islogical(options.title));
            output_filepath = fullfile(options.output_dir, sprintf('%s_spdByRegion.pdf', options.title));
            if (~exist(output_filepath, 'file') || options.overwrite_existing)
                exportgraphics(spdByRegionHandle, output_filepath, 'ContentType', 'vector');
            end
            close(spdByRegionHandle);
        end

        return
    end

    %% CASE 3: multiple projection types
    nProjections = numel(projectionNames);

    % Exponent maps side-by-side
    exponentMapHandle = figure;
    tiledlayout(1, nProjections, 'Padding', 'compact', 'TileSpacing', 'compact');

    for pp = 1:nProjections
        projectionName = projectionNames{pp};
        exponentMap = plotData.(projectionName).exponentMap;
        exponentMap(mask == 0) = nan;

        nexttile;
        imagesc(-exponentMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Exponent - %s', projectionName), 'Interpreter', 'none');
        axis square
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar
    end

    sgtitle(sprintf('Exponent Maps - %s', activityName), 'Interpreter', 'none');

    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_exponentMap.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(exponentMapHandle, output_filepath, 'ContentType', 'vector');
        end
        close(exponentMapHandle);
    end

    % Variance maps side-by-side
    varianceMapHandle = figure;
    tiledlayout(1, nProjections, 'Padding', 'compact', 'TileSpacing', 'compact');

    for pp = 1:nProjections
        projectionName = projectionNames{pp};
        varianceMap = plotData.(projectionName).varianceMap;
        varianceMap(mask == 0) = nan;

        nexttile;
        imagesc(varianceMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Contrast Variance - %s', projectionName), 'Interpreter', 'none');
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar
    end

    sgtitle(sprintf('Variance Maps - %s', activityName), 'Interpreter', 'none');

    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_varianceMap.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(varianceMapHandle, output_filepath, 'ContentType', 'vector');
        end
        close(varianceMapHandle);
    end

    % SPD-by-region on same graph
    spdByRegionHandle = figure;
    hold on

    legendEntries = {};

    for pp = 1:nProjections
        projectionName = projectionNames{pp};
        spdByRegion = plotData.(projectionName).spdByRegion;
        frq = plotData.(projectionName).frq;

        switch projectionName
            case 'justProjection'
                centerStyle = '-k';
                peripheryStyle = '-r';
            case 'virtuallyFoveated'
                centerStyle = '--k';
                peripheryStyle = '--r';
            otherwise
                centerStyle = '-';
                peripheryStyle = '--';
        end

        loglog(frq, squeeze(spdByRegion(20,20,:)), centerStyle, 'LineWidth', 1.5);
        legendEntries{end+1} = sprintf('%s center', projectionName); %#ok<AGROW>

        loglog(frq, squeeze(spdByRegion(31,20,:)), peripheryStyle, 'LineWidth', 1.5);
        legendEntries{end+1} = sprintf('%s periphery', projectionName); %#ok<AGROW>
    end

    plot([10^0 10^1.5], [10^-2 10^-5], ':k');
    legend([legendEntries, {'reference slope'}], 'Location', 'best', 'Interpreter', 'none');
    ylabel('Power [contrast^2/Hz]');
    xlabel('Frequency [log Hz]');
    title(sprintf('SPDs from the center and periphery - %s', activityName), 'Interpreter', 'none');

    if (~islogical(options.spd_xlim))
        xlim(double(options.spd_xlim));
    end
    if (~islogical(options.spd_ylim))
        ylim(double(options.spd_ylim));
    end

    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_spdByRegion.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(spdByRegionHandle, output_filepath, 'ContentType', 'vector');
        end
        close(spdByRegionHandle);
    end

end


function activityDataStruct = iLoadActivityDataStruct(matPath)
% Load a MAT file and return the activityData-style struct inside it.

    loadedData = load(matPath);

    if (isfield(loadedData, 'activityData'))
        activityDataStruct = loadedData.activityData;
        return
    end

    loadedFieldNames = fieldnames(loadedData);
    assert(~isempty(loadedFieldNames), 'Loaded MAT-file contains no variables.');

    % If there is only one variable, use it
    if (numel(loadedFieldNames) == 1)
        activityDataStruct = loadedData.(loadedFieldNames{1});
        return
    end

    % Otherwise prefer the first struct variable
    for ii = 1:numel(loadedFieldNames)
        candidate = loadedData.(loadedFieldNames{ii});
        if (isstruct(candidate))
            activityDataStruct = candidate;
            return
        end
    end

    error('Could not find an activityData-style struct in MAT-file: %s', matPath);
end


function projectionData = iEnsureProjectionDataStruct(inputData)
% Standardize one projection input so the output directly contains:
%   exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector
%
% inputData may be:
%   1. a struct already containing exponentMap etc.
%   2. a path to a MAT file containing either:
%       - a struct with exponentMap etc.
%       - an activityData-style struct with one activity field whose value
%         contains exponentMap etc.

    % If it is a path, load it
    if (ischar(inputData) || isstring(inputData))
        loadedStruct = iLoadActivityDataStruct(inputData);
    else
        loadedStruct = inputData;
    end

    assert(isstruct(loadedStruct), 'Projection input must resolve to a struct.');

    % Case: already the direct projection data struct
    if (isfield(loadedStruct, 'exponentMap'))
        projectionData = loadedStruct;
        return
    end

    % Case: an activityData-style struct with one activity field
    loadedFieldNames = fieldnames(loadedStruct);
    assert(~isempty(loadedFieldNames), 'Resolved struct has no fields.');

    if (numel(loadedFieldNames) == 1)
        candidate = loadedStruct.(loadedFieldNames{1});

        if (isstruct(candidate) && isfield(candidate, 'exponentMap'))
            projectionData = candidate;
            return
        end
    end

    error(['Could not resolve projection data struct. Expected either a struct with exponentMap ' ...
           'or a single-activity struct whose child contains exponentMap.']);
end