function [exponentMapHandle, varianceMapHandle, spdByRegionHandle, regionAverages] = plotSPDs(virtuallyFoveatedActivityData, options)
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
        options.num_participants (1,1) double {mustBePositive} = 1
        options.return_region_averages_only (1,1) logical = false
        options.acrossParticipantRegionSTD = false
    end

    % Pull optional datasets / handles out of options
    justProjectionActivityData = options.justProjectionActivityData;
    spdTargetAxes = options.spd_target_axes;

    regionAverages = struct();

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

    % Virtually foveated quantities
    virtuallyFoveatedExponentMap = virtuallyFoveatedActivityData.(activityName).exponentMap;
    virtuallyFoveatedVarianceMap = virtuallyFoveatedActivityData.(activityName).varianceMap;
    virtuallyFoveatedSpdByRegion = virtuallyFoveatedActivityData.(activityName).spdByRegion;
    virtuallyFoveatedFrq = virtuallyFoveatedActivityData.(activityName).frq;
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

   % Build center/periphery masks for SPD averaging
    center_spd_eccentricity = 5; % Degrees
    centerRadiusPix = center_spd_eccentricity / degPerPix;

    % Do this for each subject, take the std across these averages, 
    % then use this as the standard deviation in the SEM 

    spdHeight = size(virtuallyFoveatedSpdByRegion, 1);
    spdWidth = size(virtuallyFoveatedSpdByRegion, 2);

    % Region-center locations expressed in image-pixel coordinates
    spdXCenters = ((1:spdWidth) - 0.5) * (size(virtuallyFoveatedExponentMap, 2) / spdWidth);
    spdYCenters = ((1:spdHeight) - 0.5) * (size(virtuallyFoveatedExponentMap, 1) / spdHeight);

    [spdXGrid, spdYGrid] = meshgrid(spdXCenters, spdYCenters);

    imageCenterX = size(virtuallyFoveatedExponentMap, 2) / 2;
    imageCenterY = size(virtuallyFoveatedExponentMap, 1) / 2;

    spdDistanceFromCenter = sqrt((spdXGrid - imageCenterX).^2 + (spdYGrid - imageCenterY).^2);

    centerMask = (spdDistanceFromCenter <= centerRadiusPix);
    peripheryMask = (spdDistanceFromCenter > centerRadiusPix);

    % ---------------------------------------------------------------------
    % Compute region averages once up front
    % ---------------------------------------------------------------------
    virtuallyFoveatedCenterSpd = ...
        iComputeRegionMean(virtuallyFoveatedSpdByRegion, centerMask);

    virtuallyFoveatedPeripherySpd = ...
        iComputeRegionMean(virtuallyFoveatedSpdByRegion, peripheryMask);

    regionAverages.virtuallyFoveated.center = virtuallyFoveatedCenterSpd;
    regionAverages.virtuallyFoveated.periphery = virtuallyFoveatedPeripherySpd;

    if (~islogical(justProjectionActivityData))
        justProjectionCenterSpd = ...
            iComputeRegionMean(justProjectionSpdByRegion, centerMask);

        justProjectionPeripherySpd = ...
            iComputeRegionMean(justProjectionSpdByRegion, peripheryMask);

        regionAverages.justProjection.center = justProjectionCenterSpd;
        regionAverages.justProjection.periphery = justProjectionPeripherySpd;
    end

    % ---------------------------------------------------------------------
    % OPTIONAL: return region averages only (skip all plotting)
    % ---------------------------------------------------------------------
    if (options.return_region_averages_only)
        exponentMapHandle = [];
        varianceMapHandle = [];
        spdByRegionHandle = [];
        return;
    end


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
            climVals = double(options.exponent_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colormap(pink);
            colorbar
        end 

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
            climVals = double(options.exponent_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colormap(pink);
            colorbar
        end 

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
            climVals = double(options.exponent_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colormap(pink);
            colorbar
        end 
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
            climVals = double(options.variance_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colormap(pink);
            colorbar
        end 

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
            climVals = double(options.variance_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colorbar
        end 
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
            climVals = double(options.variance_clim);
            clim(climVals);
            colormap(pink);
            cb = colorbar;
            cb.Ticks = linspace(climVals(1), climVals(2), 6);
            cb.TickLabels = arrayfun(@(x) sprintf('%.3g', x), cb.Ticks, 'UniformOutput', false);
        else 
            colormap(pink);
            colorbar
        end 
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

    numParticipants = options.num_participants;
    plotSem = (~islogical(options.acrossParticipantRegionSTD)) && (numParticipants > 1);

    % Color families:
    %   justProjection     = blue family
    %   virtuallyFoveated  = red family
    %
    % Within each family:
    %   center     = dark, saturated
    %   periphery  = much lighter / desaturated

    % justProjection (blue family)
    jpCenterColor     = [0.00, 0.20, 0.70];   % deep blue
    jpPeripheryColor  = [0.55, 0.75, 0.95];   % very light blue

    % virtuallyFoveated (red family)
    vfCenterColor     = [0.75, 0.05, 0.05];   % deep red
    vfPeripheryColor  = [0.95, 0.65, 0.65];   % very light red

    % Reference slope
    referenceColor = [0.25, 0.25, 0.25];

        % DOUBLE DATASET BRANCH
    if (~islogical(justProjectionActivityData))

        justProjectionCenterSpd = ...
            iComputeRegionMean(justProjectionSpdByRegion, centerMask);

        justProjectionPeripherySpd = ...
            iComputeRegionMean(justProjectionSpdByRegion, peripheryMask);

        virtuallyFoveatedCenterSpd = ...
            iComputeRegionMean(virtuallyFoveatedSpdByRegion, centerMask);

        virtuallyFoveatedPeripherySpd = ...
            iComputeRegionMean(virtuallyFoveatedSpdByRegion, peripheryMask);

        regionAverages.justProjection.center = justProjectionCenterSpd;
        regionAverages.justProjection.periphery = justProjectionPeripherySpd;

        regionAverages.virtuallyFoveated.center = virtuallyFoveatedCenterSpd;
        regionAverages.virtuallyFoveated.periphery = virtuallyFoveatedPeripherySpd;

        justProjectionFrqVector = justProjectionFrq(:);
        virtuallyFoveatedFrqVector = virtuallyFoveatedFrq(:);

        if (plotSem)
            [jpCenterPatchHandle, jpCenterLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, justProjectionFrqVector, justProjectionCenterSpd, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.justProjection.center, numParticipants), ...
                jpCenterColor, '-', 1.8);

            [jpPeripheryPatchHandle, jpPeripheryLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, justProjectionFrqVector, justProjectionPeripherySpd, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.justProjection.periphery, numParticipants), ...
                jpPeripheryColor, '--', 1.8);

            [vfCenterPatchHandle, vfCenterLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, virtuallyFoveatedFrqVector, virtuallyFoveatedCenterSpd, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.virtuallyFoveated.center, numParticipants), ...
                vfCenterColor, '-', 1.8);

            [vfPeripheryPatchHandle, vfPeripheryLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, virtuallyFoveatedFrqVector, virtuallyFoveatedPeripherySpd, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.virtuallyFoveated.periphery, numParticipants), ...
                vfPeripheryColor, '--', 1.8);

            referenceHandle = loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':', 'Color', referenceColor, 'LineWidth', 1.5);

            legendHandles = [ ...
                jpCenterLineHandle, jpPeripheryLineHandle, ...
                vfCenterLineHandle, vfPeripheryLineHandle, ...
                referenceHandle ...
            ];

            legendEntries = { ...
                'Non-foveated center', 'Non-foveated periphery', ...
                'Foveated center', 'Foveated periphery', ...
                'Reference (1/f^2)' ...
            };
        else
            jpCenterLineHandle = loglog(targetAxes, justProjectionFrqVector, justProjectionCenterSpd, ...
                '-', 'Color', jpCenterColor, 'LineWidth', 1.8);

            jpPeripheryLineHandle = loglog(targetAxes, justProjectionFrqVector, justProjectionPeripherySpd, ...
                '--', 'Color', jpPeripheryColor, 'LineWidth', 1.8);

            vfCenterLineHandle = loglog(targetAxes, virtuallyFoveatedFrqVector, virtuallyFoveatedCenterSpd, ...
                '-', 'Color', vfCenterColor, 'LineWidth', 1.8);

            vfPeripheryLineHandle = loglog(targetAxes, virtuallyFoveatedFrqVector, virtuallyFoveatedPeripherySpd, ...
                '--', 'Color', vfPeripheryColor, 'LineWidth', 1.8);

            referenceHandle = loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':', 'Color', referenceColor, 'LineWidth', 1.5);

            legendHandles = [ ...
                jpCenterLineHandle, jpPeripheryLineHandle, ...
                vfCenterLineHandle, vfPeripheryLineHandle, ...
                referenceHandle ...
            ];

            legendEntries = { ...
                'Non-foveated center', 'Non-foveated periphery', ...
                'Foveated center', 'Foveated periphery', ...
                'Reference (1/f^2)' ...
            };
        end

        legend(targetAxes, legendHandles, legendEntries, ...
            'Location', 'best', 'Interpreter', 'none');

    % SINGLE DATASET BRANCH
    else
        centerSPD = iComputeRegionMean(virtuallyFoveatedSpdByRegion, centerMask);
        peripherySPD = iComputeRegionMean(virtuallyFoveatedSpdByRegion, peripheryMask);

        regionAverages.virtuallyFoveated.center = centerSPD;
        regionAverages.virtuallyFoveated.periphery = peripherySPD;

        frqVector = virtuallyFoveatedFrq(:);

        if (plotSem)
            [centerPatchHandle, centerLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, frqVector, centerSPD, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.virtuallyFoveated.center, numParticipants), ...
                vfCenterColor, '-', 1.8);

            [peripheryPatchHandle, peripheryLineHandle] = iPlotSpdWithSemPatch( ...
                targetAxes, frqVector, peripherySPD, ...
                iComputeRegionSemFromAcrossParticipantStd( ...
                    options.acrossParticipantRegionSTD.virtuallyFoveated.periphery, numParticipants), ...
                vfPeripheryColor, '--', 1.8);

           referenceHandle = loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':', 'Color', referenceColor, 'LineWidth', 1.5);

            legendHandles = [ ...
                centerPatchHandle, centerLineHandle, ...
                peripheryPatchHandle, peripheryLineHandle, ...
                referenceHandle ...
            ];

            legendEntries = { ...
                'center SEM', 'center', ...
                'periphery SEM', 'periphery', ...
                'reference slope' ...
            };
        else
            centerLineHandle = loglog(targetAxes, frqVector, centerSPD, ...
                '-', 'Color', vfCenterColor, 'LineWidth', 1.8);

            peripheryLineHandle = loglog(targetAxes, frqVector, peripherySPD, ...
                '--', 'Color', vfPeripheryColor, 'LineWidth', 1.8);

            referenceHandle = loglog(targetAxes, [10^0; 10^1.5], [10^-2; 10^-5], ':', 'Color', referenceColor, 'LineWidth', 1.5);

            legendHandles = [centerLineHandle, peripheryLineHandle, referenceHandle];

            legendEntries = {'center', 'periphery', 'reference slope'};
        end

        legend(targetAxes, legendHandles, legendEntries, ...
            'Location', 'best', 'Interpreter', 'none');
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

    % Get current limits
    xL = xlim(targetAxes);
    yL = ylim(targetAxes);

    %% =========================
    % X AXIS (fixed powers of 2)
    %% =========================
    desiredTicks = [1, 2, 4, 8, 16, 32, 64];
    xticksVals = desiredTicks(desiredTicks >= xL(1) & desiredTicks <= xL(2));
    set(targetAxes, 'XTick', xticksVals);
    xticklabels(targetAxes, arrayfun(@(x) sprintf('%d', x), xticksVals, 'UniformOutput', false));

    %% =========================
    % Y AXIS (1e-5, 10e-1)
    %% =========================
    yticksVals = [1e-5 1e-4 1e-3 1e-2 1e-1];
    set(targetAxes, 'YTick', yticksVals);
    yticklabels(targetAxes, { ...
        '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}' ...
    });

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

function regionMeanSpd = iComputeRegionMean(spdByRegion, regionMask)
% Compute mean SPD at each frequency for a spatial region.

    numFrequencies = size(spdByRegion, 3);
    regionMeanSpd = nan(numFrequencies, 1);

    for ff = 1:numFrequencies
        currentSlice = spdByRegion(:,:,ff);
        regionValues = currentSlice(regionMask);
        regionValues = regionValues(isfinite(regionValues));

        if (~isempty(regionValues))
            regionMeanSpd(ff) = mean(regionValues, 'omitmissing');
        end
    end
end

function regionSemSpd = iComputeRegionSemFromAcrossParticipantStd(acrossParticipantRegionSTD, numParticipants)
% Convert across-participant region STD into SEM.
%
% Expected input:
%   acrossParticipantRegionSTD is a column vector (or vector-like) with one
%   STD value per frequency.

    regionSemSpd = acrossParticipantRegionSTD(:) ./ sqrt(numParticipants);
end

function [patchHandle, lineHandle] = iPlotSpdWithSemPatch(targetAxes, frqVector, meanSpd, semSpd, lineColor, lineStyle, lineWidth)
% Plot mean SPD with a shaded SEM background band.

    frqVector = frqVector(:);
    meanSpd = meanSpd(:);
    semSpd = semSpd(:);

    validIdx = isfinite(frqVector) & isfinite(meanSpd) & isfinite(semSpd) & ...
               (frqVector > 0) & (meanSpd > 0);

    if (~any(validIdx))
        patchHandle = gobjects(1);
        lineHandle = gobjects(1);
        return;
    end

    x = frqVector(validIdx);
    y = meanSpd(validIdx);
    sem = semSpd(validIdx);

    lowerBand = y - sem;
    upperBand = y + sem;

    % Keep patch valid on log axes
    lowerBand(lowerBand <= 0) = eps;
    upperBand(upperBand <= 0) = eps;

    patchHandle = patch(targetAxes, ...
        [x; flipud(x)], ...
        [lowerBand; flipud(upperBand)], ...
        lineColor, ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'on');

    lineHandle = loglog(targetAxes, x, y, ...
        'LineStyle', lineStyle, ...
        'Color', lineColor, ...
        'LineWidth', lineWidth);
end

function mustBePositive(x)
    if (any(x <= 0))
        error('Value must be positive.');
    end
end