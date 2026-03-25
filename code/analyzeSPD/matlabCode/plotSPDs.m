function [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(virtuallyFoveatedActivityData, options)
% Plot temporal SPD exponent, variance, and regional spectrum summaries
%
% Syntax:
%   [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(virtuallyFoveatedActivityData)
%   [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(virtuallyFoveatedActivityData, options)
%
% Description:
%   This function visualizes the outputs of temporal SPD analysis for a
%   single activity stored in a virtuallyFoveatedActivityData structure.
%   It extracts the exponent map, variance map, regional SPD summaries,
%   frequency support, and other metadata for the activity, then produces
%   three figures.
%
%   When justProjectionActivityData is provided, the exponent and variance
%   plots are displayed as side-by-side comparisons:
%       left  = just projection
%       right = virtually foveated
%
%   In that same case, the SPD-by-region plot overlays all regional curves
%   from both datasets onto a single axes.
%
%   If options.spd_target_axes is provided, the SPD-by-region plot is drawn
%   directly into that axes instead of creating a new figure. This is used
%   by groupSPDs to avoid fragile copyobj-based figure copying.
%
% Inputs:
%   virtuallyFoveatedActivityData - Struct or filepath. Contains temporal
%       SPD analysis results for a single activity. The first field is
%       assumed to be the activity name and to contain:
%           .exponentMap
%           .varianceMap
%           .spdByRegion
%           .frq
%           .medianImage
%           .frameDropVector
%
% Optional key/value pairs:
%   fovDegrees                - Scalar. Field of view in degrees used to
%                               convert eccentricity ring radii from
%                               degrees into pixels for plotting.
%   exponent_clim             - false or 2-element numeric vector for CLim.
%   variance_clim             - false or 2-element numeric vector for CLim.
%   spd_xlim                  - false or 2-element numeric vector for xlim.
%   spd_ylim                  - false or 2-element numeric vector for ylim.
%   output_dir                - false or output directory path.
%   overwrite_existing        - Logical. Whether to overwrite exported files.
%   title                     - String used in exported filenames.
%   justProjectionActivityData - false, struct, or filepath containing the
%                               matching activity results for just-projection
%                               plotting comparisons.
%   spd_target_axes           - false or valid axes handle. If provided,
%                               draw the SPD plot directly into this axes.
%
% Outputs:
%   exponentMapHandle         - Figure handle for the exponent map figure.
%   varianceMapHandle         - Figure handle for the variance map figure.
%   spdByRegionHandle         - Figure handle for the SPD summary figure,
%                               or the parent figure of spd_target_axes if
%                               spd_target_axes is supplied.

    arguments
        virtuallyFoveatedActivityData {mustBeStructOrText}
        options.fovDegrees = 120
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.output_dir = false
        options.overwrite_existing = false
        options.title = ""
        options.justProjectionActivityData = false
        options.spd_target_axes = false
    end

    % Pull optional datasets / handles out of options
    justProjectionActivityData = options.justProjectionActivityData;
    spdTargetAxes = options.spd_target_axes;

    % If the virtually foveated input was passed in as a file path, load it
    if (isstring(virtuallyFoveatedActivityData) || ischar(virtuallyFoveatedActivityData))
        try
            virtuallyFoveatedActivityData = load(virtuallyFoveatedActivityData).activityData;
        catch
            error("Could not load filepath: %s\n", virtuallyFoveatedActivityData);
        end
    end

    % If the just-projection input was passed in as a file path, load it too
    if (isstring(justProjectionActivityData) || ischar(justProjectionActivityData))
        try
            justProjectionActivityData = load(justProjectionActivityData).activityData;
        catch
            error("Could not load %s\n", justProjectionActivityData);
        end
    end

    % Extract the activity name
    field_names = fieldnames(virtuallyFoveatedActivityData);
    activityName = field_names{1};

    fovDegrees = options.fovDegrees;

    disp(virtuallyFoveatedActivityData)


    % Virtually foveated quantities
    virtuallyFoveatedExponentMap = virtuallyFoveatedActivityData.(activityName).exponentMap;
    virtuallyFoveatedVarianceMap = virtuallyFoveatedActivityData.(activityName).varianceMap;
    virtuallyFoveatedSpdByRegion = virtuallyFoveatedActivityData.(activityName).spdByRegion;
    virtuallyFoveatedFrq = virtuallyFoveatedActivityData.(activityName).frq;
    medianImage = virtuallyFoveatedActivityData.(activityName).medianImage; %#ok<NASGU>
    frameDropVector = virtuallyFoveatedActivityData.(activityName).frameDropVector; %#ok<NASGU>

    % Build elliptical FOV mask
    ellipseTransparentParams = [240, 240, 120000, .75, 0];
    p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
    myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
    [X, Y] = meshgrid(1:480, 1:480);
    mask = double(myEllipse(X,Y) < 1e-9);

    virtuallyFoveatedExponentMap(mask == 0) = nan;
    virtuallyFoveatedVarianceMap(mask == 0) = nan;

    % If just-projection data exists, extract it and apply the same mask
    if (~islogical(justProjectionActivityData))
        justProjectionExponentMap = justProjectionActivityData.(activityName).exponentMap;
        justProjectionVarianceMap = justProjectionActivityData.(activityName).varianceMap;
        justProjectionSpdByRegion = justProjectionActivityData.(activityName).spdByRegion;
        justProjectionFrq = justProjectionActivityData.(activityName).frq;

        justProjectionExponentMap(mask == 0) = nan;
        justProjectionVarianceMap(mask == 0) = nan;
    end

    degPerPix = fovDegrees / size(virtuallyFoveatedExponentMap,1);
    idxStarts = 1:12:(480 - 24 + 1);

    % ---------------------------------------------------------------------
    % Exponent map figure
    % ---------------------------------------------------------------------
    exponentMapHandle = figure;

    if (~islogical(justProjectionActivityData))
        tiledlayout(1,2);

        % Left = just projection
        nexttile;
        imagesc(-justProjectionExponentMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Exponent - %s - justProjection', activityName), 'Interpreter', 'none')
        axis square
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar

        % Right = virtually foveated
        nexttile;
        imagesc(-virtuallyFoveatedExponentMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Exponent - %s - virtuallyFoveated', activityName), 'Interpreter', 'none')
        axis square
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar
    else
        imagesc(-virtuallyFoveatedExponentMap);
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
    end

    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_exponentMap.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(exponentMapHandle, output_filepath, 'ContentType', 'vector')
        end
        close(exponentMapHandle);
    end

    % ---------------------------------------------------------------------
    % Variance map figure
    % ---------------------------------------------------------------------
    varianceMapHandle = figure;

    if (~islogical(justProjectionActivityData))
        tiledlayout(1,2);

        % Left = just projection
        nexttile;
        imagesc(justProjectionVarianceMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Contrast Variance - %s - justProjection', activityName), 'Interpreter', 'none')
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar

        % Right = virtually foveated
        nexttile;
        imagesc(virtuallyFoveatedVarianceMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Contrast Variance - %s - virtuallyFoveated', activityName), 'Interpreter', 'none')
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar
    else
        imagesc(virtuallyFoveatedVarianceMap);
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
    end

    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_varianceMap.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(varianceMapHandle, output_filepath, 'ContentType', 'vector')
        end
        close(varianceMapHandle);
    end

    % ---------------------------------------------------------------------
    % SPD-by-region plot
    % ---------------------------------------------------------------------
    % Draw directly into caller-provided axes if given; otherwise make a figure.
    if (~islogical(spdTargetAxes))
        assert(isgraphics(spdTargetAxes, 'axes'), 'options.spd_target_axes must be a valid axes handle.');
        targetAxes = spdTargetAxes;
        cla(targetAxes);
        hold(targetAxes, 'on');
        spdByRegionHandle = ancestor(targetAxes, 'figure');
    else
        spdByRegionHandle = figure;
        targetAxes = axes('Parent', spdByRegionHandle);
        hold(targetAxes, 'on');
    end

    if (~islogical(justProjectionActivityData))
        virtuallyFoveatedCenterSpd = squeeze(virtuallyFoveatedSpdByRegion(20,20,:));
        virtuallyFoveatedPeripherySpd = squeeze(virtuallyFoveatedSpdByRegion(31,20,:));
        justProjectionCenterSpd = squeeze(justProjectionSpdByRegion(20,20,:));
        justProjectionPeripherySpd = squeeze(justProjectionSpdByRegion(31,20,:));

        virtuallyFoveatedFrqVector = virtuallyFoveatedFrq(:);
        justProjectionFrqVector = justProjectionFrq(:);

        virtuallyFoveatedCenterSpd = virtuallyFoveatedCenterSpd(:);
        virtuallyFoveatedPeripherySpd = virtuallyFoveatedPeripherySpd(:);
        justProjectionCenterSpd = justProjectionCenterSpd(:);
        justProjectionPeripherySpd = justProjectionPeripherySpd(:);

        virtuallyFoveatedCenterValidIdx = isfinite(virtuallyFoveatedFrqVector) & ...
                                          isfinite(virtuallyFoveatedCenterSpd) & ...
                                          (virtuallyFoveatedFrqVector > 0) & ...
                                          (virtuallyFoveatedCenterSpd > 0);

        virtuallyFoveatedPeripheryValidIdx = isfinite(virtuallyFoveatedFrqVector) & ...
                                             isfinite(virtuallyFoveatedPeripherySpd) & ...
                                             (virtuallyFoveatedFrqVector > 0) & ...
                                             (virtuallyFoveatedPeripherySpd > 0);

        justProjectionCenterValidIdx = isfinite(justProjectionFrqVector) & ...
                                       isfinite(justProjectionCenterSpd) & ...
                                       (justProjectionFrqVector > 0) & ...
                                       (justProjectionCenterSpd > 0);

        justProjectionPeripheryValidIdx = isfinite(justProjectionFrqVector) & ...
                                          isfinite(justProjectionPeripherySpd) & ...
                                          (justProjectionFrqVector > 0) & ...
                                          (justProjectionPeripherySpd > 0);

        if (any(justProjectionCenterValidIdx))
            loglog(targetAxes, ...
                   justProjectionFrqVector(justProjectionCenterValidIdx), ...
                   justProjectionCenterSpd(justProjectionCenterValidIdx), ...
                   '-k', 'LineWidth', 1.5);
        end

        if (any(justProjectionPeripheryValidIdx))
            loglog(targetAxes, ...
                   justProjectionFrqVector(justProjectionPeripheryValidIdx), ...
                   justProjectionPeripherySpd(justProjectionPeripheryValidIdx), ...
                   '-r', 'LineWidth', 1.5);
        end

        if (any(virtuallyFoveatedCenterValidIdx))
            loglog(targetAxes, ...
                   virtuallyFoveatedFrqVector(virtuallyFoveatedCenterValidIdx), ...
                   virtuallyFoveatedCenterSpd(virtuallyFoveatedCenterValidIdx), ...
                   '--k', 'LineWidth', 1.5);
        end

        if (any(virtuallyFoveatedPeripheryValidIdx))
            loglog(targetAxes, ...
                   virtuallyFoveatedFrqVector(virtuallyFoveatedPeripheryValidIdx), ...
                   virtuallyFoveatedPeripherySpd(virtuallyFoveatedPeripheryValidIdx), ...
                   '--r', 'LineWidth', 1.5);
        end

        loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':k');

        legend(targetAxes, ...
            {'justProjection center', 'justProjection periphery', ...
             'virtuallyFoveated center', 'virtuallyFoveated periphery', ...
             'reference slope'}, ...
            'Location', 'best', 'Interpreter', 'none');

    else
        centerSPD = squeeze(virtuallyFoveatedSpdByRegion(20,20,:));
        peripherySPD = squeeze(virtuallyFoveatedSpdByRegion(31,20,:));

        frqVector = virtuallyFoveatedFrq(:);
        centerSPD = centerSPD(:);
        peripherySPD = peripherySPD(:);

        validCenter = isfinite(frqVector) & isfinite(centerSPD) & (frqVector > 0) & (centerSPD > 0);
        validPeriphery = isfinite(frqVector) & isfinite(peripherySPD) & (frqVector > 0) & (peripherySPD > 0);

        if (any(validCenter))
            loglog(targetAxes, frqVector(validCenter), centerSPD(validCenter), '-k', 'LineWidth', 1.5);
        end
        if (any(validPeriphery))
            loglog(targetAxes, frqVector(validPeriphery), peripherySPD(validPeriphery), '-r', 'LineWidth', 1.5);
        end

        loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':k');

        legend(targetAxes, {'center','periphery','reference slope'}, 'Location', 'best');
    end

    set(targetAxes, 'XScale', 'log', 'YScale', 'log');
    ylabel(targetAxes, 'Power [contrast^2/Hz]');
    xlabel(targetAxes, 'Frequency [Hz]');
    title(targetAxes, sprintf('SPDs from the center and periphery - %s', activityName), 'Interpreter', 'none');

    if (~islogical(options.spd_xlim))
        spdXLim = double(options.spd_xlim);
        spdXLim(spdXLim <= 0) = 1e-12;
        xlim(targetAxes, sort(spdXLim));
    end
    if (~islogical(options.spd_ylim))
        spdYLim = double(options.spd_ylim);
        spdYLim(spdYLim <= 0) = 1e-12;
        ylim(targetAxes, sort(spdYLim));
    end

    % Only export/close the SPD figure if we created it locally
    if (islogical(spdTargetAxes) && ~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_spdByRegion.pdf', options.title));
        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(spdByRegionHandle, output_filepath, 'ContentType', 'vector');
        end
        close(spdByRegionHandle);
    end
end

function mustBeStructOrText(x)
    if ~(isstruct(x) || ischar(x) || isstring(x))
        error('Input must be a struct, char, or string.');
    end
end