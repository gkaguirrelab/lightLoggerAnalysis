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
%   fovDegrees               - Scalar. Field of view in degrees used to
%                              convert eccentricity ring radii from
%                              degrees into pixels for plotting.
%   exponent_clim            - false or 2-element numeric vector for CLim.
%   variance_clim            - false or 2-element numeric vector for CLim.
%   spd_xlim                 - false or 2-element numeric vector for xlim.
%   spd_ylim                 - false or 2-element numeric vector for ylim.
%   output_dir               - false or output directory path.
%   overwrite_existing       - Logical. Whether to overwrite exported files.
%   title                    - String used in exported filenames.
%   justProjectionActivityData - false, struct, or filepath containing the
%                              matching activity results for just-projection
%                              plotting comparisons.
%
% Outputs:
%   exponentMapHandle        - Figure handle for the exponent map figure.
%   varianceMapHandle        - Figure handle for the variance map figure.
%   spdByRegionHandle        - Figure handle for the SPD summary figure.
%
% Examples:
%{
    load('/path/to/FLIC_2001_walkIndoorFoveate_SPDResults.mat', 'virtuallyFoveatedActivityData');

    [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = ...
        plotSPDs(virtuallyFoveatedActivityData, "fovDegrees", 120);
%}

    arguments
        virtuallyFoveatedActivityData;  % May be either the activity-data struct or a path to a .mat file containing activityData
        options.fovDegrees = 120;
        options.exponent_clim = false;
        options.variance_clim = false;
        options.spd_xlim = false;
        options.spd_ylim = false;
        options.output_dir = false;
        options.overwrite_existing = false;
        options.title = "";
        options.justProjectionActivityData = false;
    end

    % Pull the optional just-projection dataset out of the options struct so
    % we can treat it the same way as the virtually foveated data below.
    justProjectionActivityData = options.justProjectionActivityData;

    % If the virtually foveated input was passed in as a file path, load the
    % activityData variable from disk now.
    if (isstring(virtuallyFoveatedActivityData) || ischar(virtuallyFoveatedActivityData))
        virtuallyFoveatedActivityData = load(virtuallyFoveatedActivityData).activityData;
    end

    % If the just-projection input was passed in as a file path, load that
    % activityData variable from disk as well.
    if (isstring(justProjectionActivityData) || ischar(justProjectionActivityData))
        justProjectionActivityData = load(justProjectionActivityData).activityData;
    end

    % The activity-data struct is expected to contain exactly one top-level
    % field corresponding to the activity name. Pull that field name out so
    % the remaining code can access the data consistently.
    field_names = fieldnames(virtuallyFoveatedActivityData);
    activityName = field_names{1};

    % Cache the field-of-view setting locally for plotting conversions.
    fovDegrees = options.fovDegrees;

    % Extract all of the virtually foveated quantities using descriptive
    % variable names so it is always obvious which dataset they belong to.
    virtuallyFoveatedExponentMap = virtuallyFoveatedActivityData.(activityName).exponentMap;
    virtuallyFoveatedVarianceMap = virtuallyFoveatedActivityData.(activityName).varianceMap;
    virtuallyFoveatedSpdByRegion = virtuallyFoveatedActivityData.(activityName).spdByRegion;
    virtuallyFoveatedFrq = virtuallyFoveatedActivityData.(activityName).frq;
    medianImage = virtuallyFoveatedActivityData.(activityName).medianImage;
    frameDropVector = virtuallyFoveatedActivityData.(activityName).frameDropVector;

    % Build an elliptical mask corresponding to the valid field of view.
    % Values outside this ellipse are set to NaN so they do not appear in
    % the exponent/variance images.
    ellipseTransparentParams = [240, 240, 120000, .75, 0];
    p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
    myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
    [X, Y] = meshgrid(1:480, 1:480);
    mask = double(myEllipse(X,Y) < 1e-9);

    % Apply the mask to the virtually foveated maps.
    virtuallyFoveatedExponentMap(mask == 0) = nan;
    virtuallyFoveatedVarianceMap(mask == 0) = nan;

    % If just-projection data was provided, extract the analogous variables
    % and apply the same field-of-view mask so the comparison is visually
    % consistent across the two conditions.
    if (~islogical(justProjectionActivityData))
        justProjectionExponentMap = justProjectionActivityData.(activityName).exponentMap;
        justProjectionVarianceMap = justProjectionActivityData.(activityName).varianceMap;
        justProjectionSpdByRegion = justProjectionActivityData.(activityName).spdByRegion;
        justProjectionFrq = justProjectionActivityData.(activityName).frq;

        justProjectionExponentMap(mask == 0) = nan;
        justProjectionVarianceMap(mask == 0) = nan;
    end

    % Convert degrees of visual angle into pixels so that the concentric
    % eccentricity circles are drawn at the intended radii on the map.
    degPerPix = fovDegrees / size(virtuallyFoveatedExponentMap,1);

    % Precompute the top-left pixel indices of the 39x39 regional grid used
    % for overlaying the patch boundaries on the map images.
    idxStarts = 1:12:(480 - 24 + 1);

    % ---------------------------------------------------------------------
    % Exponent map figure
    % ---------------------------------------------------------------------
    % If just-projection data exists, show a 1x2 comparison:
    %   left  = just projection
    %   right = virtually foveated
    %
    % Otherwise preserve the original single-map behavior.
    exponentMapHandle = figure;

    if (~islogical(justProjectionActivityData))
        tiledlayout(1,2);

        % ---------------- Left tile: just projection exponent map ---------
        nexttile;
        imagesc(-justProjectionExponentMap);
        hold on

        % Mark the image center and overlay eccentricity rings plus the
        % underlying 39x39 regional grid.
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end

        title(sprintf('Exponent - Just Projection - %s', activityName));
        axis square

        % Only apply a custom color limit if the user supplied one.
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar

        % ------------- Right tile: virtually foveated exponent map --------
        nexttile;
        imagesc(-virtuallyFoveatedExponentMap);
        hold on
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end

        title(sprintf('Exponent - Virtually Foveated - %s', activityName));
        axis square

        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar

    else
        % Original single-dataset exponent map behavior.
        imagesc(-virtuallyFoveatedExponentMap);
        hold on
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Exponent - %s', activityName));
        axis square
        if (~islogical(options.exponent_clim))
            clim(double(options.exponent_clim));
        end
        colorbar
    end

    % Export the exponent figure if an output directory was supplied.
    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_exponentMap.pdf', options.title));
        if (~exist(output_filepath) || options.overwrite_existing)
            exportgraphics(exponentMapHandle, output_filepath, 'ContentType', 'vector')
        end
        close(exponentMapHandle);
    end

    % ---------------------------------------------------------------------
    % Variance map figure
    % ---------------------------------------------------------------------
    % Same logic as above: two-panel comparison if just-projection data is
    % available, otherwise preserve the original single-map behavior.
    varianceMapHandle = figure;

    if (~islogical(justProjectionActivityData))
        tiledlayout(1,2);

        % ---------------- Left tile: just projection variance map ---------
        nexttile;
        imagesc(justProjectionVarianceMap, [0.015 0.035]);
        hold on
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end

        title(sprintf('Contrast Variance - Just Projection - %s', activityName));
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar

        % ----------- Right tile: virtually foveated variance map ----------
        nexttile;
        imagesc(virtuallyFoveatedVarianceMap, [0.015 0.035]);
        hold on
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end

        title(sprintf('Contrast Variance - Virtually Foveated - %s', activityName));
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar

    else
        % Original single-dataset variance map behavior.
        imagesc(virtuallyFoveatedVarianceMap, [0.015 0.035]);
        hold on
        plot(240, 240, '+k');
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii / degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Contrast Variance - %s', activityName));
        axis square
        if (~islogical(options.variance_clim))
            clim(double(options.variance_clim));
        end
        colorbar
    end

    % Export the variance figure if requested.
    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_varianceMap.pdf', options.title));
        if (~exist(output_filepath) || options.overwrite_existing)
            exportgraphics(varianceMapHandle, output_filepath, 'ContentType', 'vector')
        end
        close(varianceMapHandle);
    end

        % ---------------------------------------------------------------------
    % SPD-by-region figure
    % ---------------------------------------------------------------------
    % Keep the original plotting style:
    %   - loglog plot
    %   - same axis labels
    %   - same title format
    %   - same reference line
    %   - same black/red center/periphery styling
    %
    % The only change is:
    %   - if justProjectionActivityData exists, also plot its center and
    %     periphery curves on the same axes.
    spdByRegionHandle = figure;

    if (~islogical(justProjectionActivityData))
        % Extract center/periphery SPD curves
        virtuallyFoveatedCenterSpd = squeeze(virtuallyFoveatedSpdByRegion(20,20,:));
        virtuallyFoveatedPeripherySpd = squeeze(virtuallyFoveatedSpdByRegion(31,20,:));
        justProjectionCenterSpd = squeeze(justProjectionSpdByRegion(20,20,:));
        justProjectionPeripherySpd = squeeze(justProjectionSpdByRegion(31,20,:));

        % Force column vectors (prevents logical indexing bugs)
        virtuallyFoveatedFrqVector = virtuallyFoveatedFrq(:);
        justProjectionFrqVector = justProjectionFrq(:);

        virtuallyFoveatedCenterSpd = virtuallyFoveatedCenterSpd(:);
        virtuallyFoveatedPeripherySpd = virtuallyFoveatedPeripherySpd(:);
        justProjectionCenterSpd = justProjectionCenterSpd(:);
        justProjectionPeripherySpd = justProjectionPeripherySpd(:);

        % Build valid masks
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

        % -----------------------------------------------------------------
        % SAFETY CHECKS
        % -----------------------------------------------------------------
        virtuallyFoveatedHasValid = any(virtuallyFoveatedCenterValidIdx) || ...
                                    any(virtuallyFoveatedPeripheryValidIdx);

        justProjectionHasValid = any(justProjectionCenterValidIdx) || ...
                                any(justProjectionPeripheryValidIdx);

        if (~virtuallyFoveatedHasValid)
            error('plotSPDs:NoValidVirtuallyFoveatedData', ...
                    'No valid SPD data for virtuallyFoveated (all values invalid or <= 0).');
        end

        if (~justProjectionHasValid)
            error('plotSPDs:NoValidJustProjectionData', ...
                    'No valid SPD data for justProjection (all values invalid or <= 0).');
        end

        % -----------------------------------------------------------------
        % Plot 
        % -----------------------------------------------------------------
        if (any(virtuallyFoveatedCenterValidIdx))
            loglog(virtuallyFoveatedFrqVector(virtuallyFoveatedCenterValidIdx), ...
                    virtuallyFoveatedCenterSpd(virtuallyFoveatedCenterValidIdx), '-k');
            hold on
        else
            hold on
        end

        if (any(virtuallyFoveatedPeripheryValidIdx))
            loglog(virtuallyFoveatedFrqVector(virtuallyFoveatedPeripheryValidIdx), ...
                    virtuallyFoveatedPeripherySpd(virtuallyFoveatedPeripheryValidIdx), '-r');
        end

        if (any(justProjectionCenterValidIdx))
            loglog(justProjectionFrqVector(justProjectionCenterValidIdx), ...
                    justProjectionCenterSpd(justProjectionCenterValidIdx), '--k');
        end

        if (any(justProjectionPeripheryValidIdx))
            loglog(justProjectionFrqVector(justProjectionPeripheryValidIdx), ...
                    justProjectionPeripherySpd(justProjectionPeripheryValidIdx), '--r');
        end

        plot([10^0 10^1.5], [10^-2 10^-5], ':k')

        ylabel('Power [contrast^2/Hz]');
        xlabel('Frequency [log Hz]');
        title(sprintf('SPDs from the center and periphery - %s', activityName));

        legend({'virtuallyFoveated center', ...
                'virtuallyFoveated periphery', ...
                'justProjection center', ...
                'justProjection periphery'});

    else
        % Preserve the exact original single-dataset behavior.
        virtuallyFoveatedCenterSpd = squeeze(virtuallyFoveatedSpdByRegion(20,20,:));
        virtuallyFoveatedPeripherySpd = squeeze(virtuallyFoveatedSpdByRegion(31,20,:));

        virtuallyFoveatedFrqVector = virtuallyFoveatedFrq(:);
        virtuallyFoveatedCenterSpd = virtuallyFoveatedCenterSpd(:);
        virtuallyFoveatedPeripherySpd = virtuallyFoveatedPeripherySpd(:);

        virtuallyFoveatedCenterValidIdx = isfinite(virtuallyFoveatedFrqVector) & ...
                                          isfinite(virtuallyFoveatedCenterSpd) & ...
                                          (virtuallyFoveatedFrqVector > 0) & ...
                                          (virtuallyFoveatedCenterSpd > 0);

        virtuallyFoveatedPeripheryValidIdx = isfinite(virtuallyFoveatedFrqVector) & ...
                                             isfinite(virtuallyFoveatedPeripherySpd) & ...
                                             (virtuallyFoveatedFrqVector > 0) & ...
                                             (virtuallyFoveatedPeripherySpd > 0);

        if (any(virtuallyFoveatedCenterValidIdx))
            loglog(virtuallyFoveatedFrqVector(virtuallyFoveatedCenterValidIdx), ...
                   virtuallyFoveatedCenterSpd(virtuallyFoveatedCenterValidIdx), '-k');
            hold on
        else
            hold on
        end

        if (any(virtuallyFoveatedPeripheryValidIdx))
            loglog(virtuallyFoveatedFrqVector(virtuallyFoveatedPeripheryValidIdx), ...
                   virtuallyFoveatedPeripherySpd(virtuallyFoveatedPeripheryValidIdx), '-r');
        end

        plot([10^0 10^1.5], [10^-2 10^-5], ':k')
        legend({'center','periphery'});
        ylabel('Power [contrast^2/Hz]');
        xlabel('Frequency [log Hz]');
        title(sprintf('SPDs from the center and periphery - %s', activityName));
    end

    if (~islogical(options.spd_xlim))
        xlim(double(options.spd_xlim));
    end
    if (~islogical(options.spd_ylim))
        ylim(double(options.spd_ylim));
    end

    % Export the SPD figure if requested.
    if (~islogical(options.output_dir))
        assert(~islogical(options.title));
        output_filepath = fullfile(options.output_dir, sprintf('%s_spdByRegion.pdf', options.title));
        if (~exist(output_filepath) || options.overwrite_existing)
            exportgraphics(spdByRegionHandle, output_filepath, 'ContentType', 'vector')
        end
        close(spdByRegionHandle);
    end

end