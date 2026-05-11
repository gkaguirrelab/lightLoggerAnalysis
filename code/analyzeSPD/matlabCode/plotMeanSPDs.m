function plotMeanSPDs(spds, options)
% Plot averaged SPD summaries using the combineSPDs visual style
%
% Syntax:
%   plotMeanSPDs(spds)
%   plotMeanSPDs(spds, options)
%
% Description:
%   This function plots mean SPD summaries that were previously generated
%   by `meanSPDs.m`. For now, this function only implements the
%   "acrossSubject" plotting case, where each activity has a single
%   averaged SPD result per color mode. The plotting style intentionally
%   mirrors the line styling, legend formatting, reference line, and
%   formula annotation used in `combineSPDs.m`, but uses the already
%   averaged SPD curves rather than recomputing them from per-subject
%   files.
%
% Inputs:
%   spds                  - Struct or filepath. If a filepath is passed,
%                           it should point to a `.mat` file containing a
%                           top-level `spds` variable produced by the
%                           Python `plot_mean_spds` helper.
%
% Optional key/value pairs:
%   output_dir            - String. Target output directory. If empty,
%                           figures are created but not exported. The
%                           directory is created automatically if needed.
%
% Notes:
%   The current Python caller may provide either of these top-level
%   structures:
%
%       Intended:
%           spds.(color_mode).(activity_name) = "/path/to/meanSPDs.mat"
%
%       Flat fallback:
%           spds.(activity_name) = "/path/to/meanSPDs.mat"
%
%   This MATLAB function accepts both forms and normalizes them into a
%   common internal representation before plotting.

    arguments 
        spds
        options.output_dir = ""
    end 

    % If the input is a filepath, load the top-level `spds` variable from
    % disk so we can work with the MATLAB struct directly.
    if (isstring(spds) || ischar(spds))
        spds = load(spds).spds;
    end

    % Normalize the incoming Python-generated structure into a clean
    % `color_mode -> activity_name -> avg_spd_path` mapping. This lets the
    % plotting code below remain simple even if the Python side provides a
    % slightly flatter dictionary than intended.
    loaded_spds = iNormalizeAcrossSubjectInput(spds);

    % Gather the color modes we will plot. Each color mode will contribute
    % its own set of projection-type curves to the final figure.
    color_modes = fieldnames(loaded_spds);
    assert(~isempty(color_modes));

    % Let's get the activity names that were passed in from the 
    % first colormode 
    activity_names = fieldnames(loaded_spds.(color_modes{1})); 
    
    % Iterate over the activities and generate a combine-style SPD figure
    % for each one.
    for aa = 1:numel(activity_names)
        activity_name = activity_names{aa};

        % Build a lightweight activity struct of the form expected by the
        % local plotting helper:
        %
        %   activity_struct.(color_mode).(projection_type) = avgSPDStruct
        %
        activity_struct = struct();
        for cc = 1:numel(color_modes)
            % Retrieve the colormode 
            color_mode = color_modes{cc};

            % Load in the SPD
            avgSPDStruct = loaded_spds.(color_mode).(activity_name);
            activity_struct.(color_mode) = avgSPDStruct.mean;
        end

        % Plot this activity using the same visual conventions as
        % combineSPDs.
        figure_handle = iPlotMeanActivity(activity_struct, activity_name);

        % Also plot exponent and variance maps for each color mode using
        % the same side-by-side map logic used elsewhere in the SPD
        % plotting code. The helper returns the full handle struct so the
        % export loop below can simply iterate over it.
        map_handles = iPlotMeanMaps(activity_struct, activity_name);

        % Output the figure if requested. `options.output_dir` is treated
        % as a directory path, so create it if it does not already exist.
        if (~(options.output_dir == ""))
            if (~isfolder(options.output_dir))
                mkdir(options.output_dir);
            end

            output_filepath = fullfile(options.output_dir, sprintf('%s_spdByRegion.pdf', activity_name));
            exportgraphics(figure_handle, output_filepath, 'ContentType', 'vector');
            close(figure_handle);

            for cc = 1:numel(color_modes)
                color_mode = color_modes{cc};
                exponent_output_filepath = fullfile(options.output_dir, sprintf('%s_%s_exponentMap.pdf', activity_name, color_mode));
                variance_output_filepath = fullfile(options.output_dir, sprintf('%s_%s_varianceMap.pdf', activity_name, color_mode));

                exportgraphics(map_handles.(color_mode).exponentMap, exponent_output_filepath, 'ContentType', 'vector');
                exportgraphics(map_handles.(color_mode).varianceMap, variance_output_filepath, 'ContentType', 'vector');
                close(map_handles.(color_mode).exponentMap);
                close(map_handles.(color_mode).varianceMap);
            end
        end
    end

end


function normalized_spds = iNormalizeAcrossSubjectInput(spds)
% Normalize the Python-generated acrossSubject input into:
%
%   normalized_spds.(color_mode).(activity_name) = avg_spd_path

    normalized_spds = struct();
    color_modes = fieldnames(spds);

    % First, detect whether a top-level field already looks like a color
    % mode whose value is itself a struct of activities.
    recognized_color_modes = {"a", "c_lm", "c_s", "L-M", "L+M+S"};
    for ff = 1:numel(color_modes)
        color_mode = color_modes{ff};
        activities = spds.(color_mode);

        % Intended form: top-level color mode field whose children are
        % activities mapped to meanSPDs.mat filepaths OR spds themselves.
        activity_names = fieldnames(activities); 
        for aa = 1:numel(activity_names)
            activity_name = activity_names{aa};
            activity_value = activities.(activity_name);
            if (ischar(activity_value) || isstring(activity_value))
                normalized_spds.(color_mode).(activity_name) = load(activity_value).avgSPDStruct;
            end
        end

    end

    
end


function figure_handle = iPlotMeanActivity(activity_struct, activity_name)
% Plot one averaged activity using the same combineSPDs styling.

    axis_font_size = 14;
    label_font_size = 18;
    title_font_size = 20;
    legend_font_size = 15;

    fig = figure('Color', 'w', 'Position', [100 100 1100 650]);
    ax = axes('Parent', fig, 'Position', [0.08 0.14 0.60 0.76]);
    hold(ax, 'on');
    grid(ax, 'off');
    box(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    ax.FontSize = axis_font_size;

    legend_handles = gobjects(0);
    legend_labels = {};

    color_modes = fieldnames(activity_struct);
    has_multiple_color_modes = numel(color_modes) > 1;
    for cc = 1:numel(color_modes)
        color_mode = color_modes{cc};
        has_added_compact_color_mode_legend_item = false;
        projection_struct = activity_struct.(color_mode);
        projection_types = fieldnames(projection_struct);

        for pp = 1:numel(projection_types)
            projection_type = projection_types{pp};
            avg_spd_struct = projection_struct.(projection_type);

            [frq, region_averages] = iLoadMeanSpdRegionAverages(avg_spd_struct);
            if (isempty(frq) || isempty(region_averages))
                continue;
            end

            line_style = iGetProjectionLineStyle(projection_type);
            base_color = iGetColorModeColor(color_mode);
            center_color = base_color;
            periphery_color = iLightenColor(base_color, 0.55);

            [center_freq, center_spd] = iCleanSpdForPlot(frq, region_averages.center);
            if (~isempty(center_freq))
                center_handle = loglog(ax, center_freq, center_spd, 'Color', center_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                if (~has_multiple_color_modes || cc == 1)
                    legend_handles(end+1) = center_handle; %#ok<AGROW>
                    legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'center'); %#ok<AGROW>
                elseif (~has_added_compact_color_mode_legend_item)
                    legend_handles(end+1) = center_handle; %#ok<AGROW>
                    legend_labels{end+1} = iFormatColorModeLabel(color_mode); %#ok<AGROW>
                    has_added_compact_color_mode_legend_item = true;
                end
            end

            [periphery_freq, periphery_spd] = iCleanSpdForPlot(frq, region_averages.periphery);
            if (~isempty(periphery_freq))
                periphery_handle = loglog(ax, periphery_freq, periphery_spd, 'Color', periphery_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                if (~has_multiple_color_modes || cc == 1)
                    legend_handles(end+1) = periphery_handle; %#ok<AGROW>
                    legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'periphery'); %#ok<AGROW>
                end
            end
        end
    end

    xlabel(ax, 'Frequency (Hz)', 'FontSize', label_font_size);
    ylabel(ax, 'Spectral power density (contrast^2/Hz)', 'FontSize', label_font_size);
    title(ax, activity_name, 'Interpreter', 'none', 'FontSize', title_font_size);
    axis(ax, 'square');
    xlim(ax, [1 64]);
    ylim(ax, [1e-10 1e-1]);
    xticks(ax, [1 2 4 8 16 32 64]);
    xticklabels(ax, {'1', '2', '4', '8', '16', '32', '64'});

    reference_handle = iPlotReferenceLine(ax);
    if (isgraphics(reference_handle))
        legend_handles(end+1) = reference_handle;
        legend_labels{end+1} = '1/f^2';
    end

    if (~isempty(legend_handles))
        legend_handle = legend(ax, legend_handles, legend_labels, ...
            'Interpreter', 'tex', 'Location', 'northeastoutside');
        legend_handle.FontSize = legend_font_size;
    end

    annotation(fig, 'textbox', [0.72 0.08 0.24 0.40], ...
        'String', iGetFormulaAnnotationText(), ...
        'Interpreter', 'tex', ...
        'FitBoxToText', 'off', ...
        'EdgeColor', [0.7 0.7 0.7], ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    figure_handle = fig;
end


function [frq, region_averages] = iLoadMeanSpdRegionAverages(avg_spd_struct)
% Extract the averaged frequency support and compute center/periphery
% averages from the already-averaged `spdByRegion` map.

    frq = avg_spd_struct.frq(:);
    spd_by_region = avg_spd_struct.spdByRegion;
    exponent_map = avg_spd_struct.exponentMap;

    if (isempty(frq) || isempty(spd_by_region) || isempty(exponent_map))
        region_averages = [];
        return;
    end

    fov_degrees = 120;
    deg_per_pix = fov_degrees / size(exponent_map, 1);
    center_spd_eccentricity = 5;
    center_radius_pix = center_spd_eccentricity / deg_per_pix;

    spd_height = size(spd_by_region, 1);
    spd_width = size(spd_by_region, 2);

    spd_x_centers = ((1:spd_width) - 0.5) * (size(exponent_map, 2) / spd_width);
    spd_y_centers = ((1:spd_height) - 0.5) * (size(exponent_map, 1) / spd_height);
    [spd_x_grid, spd_y_grid] = meshgrid(spd_x_centers, spd_y_centers);

    image_center_x = size(exponent_map, 2) / 2;
    image_center_y = size(exponent_map, 1) / 2;
    spd_distance_from_center = sqrt((spd_x_grid - image_center_x).^2 + (spd_y_grid - image_center_y).^2);

    center_mask = (spd_distance_from_center <= center_radius_pix);
    periphery_mask = (spd_distance_from_center > center_radius_pix);

    region_averages.center = iComputeRegionMean(spd_by_region, center_mask);
    region_averages.periphery = iComputeRegionMean(spd_by_region, periphery_mask);
end


function map_handles = iPlotMeanMaps(activity_struct, activity_name)
% Plot exponent and variance maps for each color mode using the averaged
% justProjection / virtuallyFoveated mean payloads. Return the figure
% handles in a struct so the caller can export and close them outside.

    color_modes = fieldnames(activity_struct);
    map_handles = struct();

    for cc = 1:numel(color_modes)
        color_mode = color_modes{cc};
        color_mode_struct = activity_struct.(color_mode);

        just_projection_exponent_map = color_mode_struct.justProjection.exponentMap;
        just_projection_variance_map = color_mode_struct.justProjection.varianceMap;
        virtually_foveated_exponent_map = color_mode_struct.virtuallyFoveated.exponentMap;
        virtually_foveated_variance_map = color_mode_struct.virtuallyFoveated.varianceMap;

        ellipseTransparentParams = [240, 240, 120000, .75, 0];
        p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
        myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
        [X, Y] = meshgrid(1:480, 1:480);
        mask = double(myEllipse(X,Y) < 1e-9);

        just_projection_exponent_map(mask == 0) = nan;
        just_projection_variance_map(mask == 0) = nan;
        virtually_foveated_exponent_map(mask == 0) = nan;
        virtually_foveated_variance_map(mask == 0) = nan;

        fov_degrees = 120;
        deg_per_pix = fov_degrees / size(virtually_foveated_exponent_map, 1);

        map_handles.(color_mode).exponentMap = figure('Color', 'w');
        tiledlayout(1,2);

        nexttile;
        imagesc(-just_projection_exponent_map);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/deg_per_pix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        title(sprintf('Exponent - %s - %s - justProjection', activity_name, color_mode), 'Interpreter', 'none')
        axis square
        set(gca, 'XTick', [], 'YTick', []);
        colormap(hot);
        cb = colorbar;
        cb.Ticks = iGetQuarterStepTicks(clim);
        cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

        nexttile;
        imagesc(-virtually_foveated_exponent_map);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/deg_per_pix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        title(sprintf('Exponent - %s - %s - virtuallyFoveated', activity_name, color_mode), 'Interpreter', 'none')
        axis square
        set(gca, 'XTick', [], 'YTick', []);
        colormap(hot);
        cb = colorbar;
        cb.Ticks = iGetQuarterStepTicks(clim);
        cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

        map_handles.(color_mode).varianceMap = figure('Color', 'w');
        tiledlayout(1,2);

        nexttile;
        imagesc(just_projection_variance_map);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/deg_per_pix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        title(sprintf('Contrast Variance - %s - %s - justProjection', activity_name, color_mode), 'Interpreter', 'none')
        axis square
        set(gca, 'XTick', [], 'YTick', []);
        colormap(hot);
        cb = colorbar;
        cb.Ticks = iGetQuarterStepTicks(clim);
        cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);

        nexttile;
        imagesc(virtually_foveated_variance_map);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/deg_per_pix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        title(sprintf('Contrast Variance - %s - %s - virtuallyFoveated', activity_name, color_mode), 'Interpreter', 'none')
        axis square
        set(gca, 'XTick', [], 'YTick', []);
        colormap(hot);
        cb = colorbar;
        cb.Ticks = iGetQuarterStepTicks(clim);
        cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);
    end
end


function line_style = iGetProjectionLineStyle(projection_type)
    switch string(projection_type)
        case "virtuallyFoveated"
            line_style = '-';
        case "justProjection"
            line_style = '--';
        otherwise
            line_style = '-';
    end
end


function color = iGetColorModeColor(color_mode)
    switch string(color_mode)
        case {"L+M", "L+M+S", "a"}
            color = [0 0 0];
        case {"L-M", "c_lm"}
            color = [0.82 0.14 0.14];
        case {"S", "c_s"}
            color = [0.12 0.35 0.85];
        otherwise
            color = [0.25 0.25 0.25];
    end
end


function lighter_color = iLightenColor(color, amount)
    lighter_color = color + (1 - color) .* amount;
end


function [frq_clean, spd_clean] = iCleanSpdForPlot(frq, spd)
    frq = frq(:);
    spd = spd(:);

    n = min(numel(frq), numel(spd));
    frq = frq(1:n);
    spd = spd(1:n);

    valid = isfinite(frq) & isfinite(spd) & frq > 0 & spd > 0;
    frq_clean = frq(valid);
    spd_clean = spd(valid);
end


function region_mean_spd = iComputeRegionMean(spd_by_region, region_mask)
    num_frequencies = size(spd_by_region, 3);
    region_mean_spd = nan(num_frequencies, 1);

    for ff = 1:num_frequencies
        current_slice = spd_by_region(:, :, ff);
        region_values = current_slice(region_mask);
        region_values = region_values(isfinite(region_values));

        if (~isempty(region_values))
            region_mean_spd(ff) = mean(region_values, 'omitmissing');
        end
    end
end


function reference_handle = iPlotReferenceLine(ax)
    drawnow;
    axis_limits = axis(ax);
    reference_freq = axis_limits(1:2).';
    reference_spd = 10^-2 .* (reference_freq .^ -2);
    reference_handle = loglog(ax, reference_freq, reference_spd, ':', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1.5);
end


function legend_label = iFormatLegendLabel(color_mode, projection_type, region_name)
    color_mode_label = iFormatColorModeLabel(color_mode);
    projection_label = iFormatProjectionLabel(projection_type);
    legend_label = sprintf('%s %s %s', color_mode_label, projection_label, region_name);
end


function color_mode_label = iFormatColorModeLabel(color_mode)
    switch string(color_mode)
        case "a"
            color_mode_label = '\ita\rm';
        case "c_lm"
            color_mode_label = '\itc_{LM}\rm';
        case "c_s"
            color_mode_label = '\itc_{S}\rm';
        case "L+M"
            color_mode_label = '\itL\rm+\itM\rm';
        case "L+M+S"
            color_mode_label = '\itL\rm+\itM\rm+\itS\rm';
        case "L-M"
            color_mode_label = '\itL\rm-\itM\rm';
        case "S"
            color_mode_label = '\itS\rm';
        otherwise
            color_mode_label = color_mode;
    end
end


function projection_label = iFormatProjectionLabel(projection_type)
    switch string(projection_type)
        case "virtuallyFoveated"
            projection_label = 'Foveated';
        case "justProjection"
            projection_label = 'Non-foveated';
        otherwise
            projection_label = projection_type;
    end
end


function formula_text = iGetFormulaAnnotationText()
    formula_text = strjoin({
        '\bfFormula definitions\rm'
        ''
        'Cone preprocessing:'
        '\iti\rm_{hat} = log(\iti\rm) - mean(log(\iti\rm))'
        ''
        'Achromatic:'
        '\ita\rm = (\itl\rm_{hat} + \itm\rm_{hat}) / \surd2'
        ''
        'Red-green opponent:'
        '\itc\rm_{LM} = (\itl\rm_{hat} - \itm\rm_{hat}) / \surd2'
        ''
        'Blue-yellow opponent:'
        '\itc\rm_{S} = (2\its\rm_{hat} - (\itl\rm_{hat} + \itm\rm_{hat})) / \surd6'
        ''
        'Source: van Hateren et al.'
        }, newline);
end


function ticks = iGetQuarterStepTicks(climVals)
% Build colorbar ticks at 0.25 increments across the provided limits.

    tickStart = ceil(climVals(1) / 0.25) * 0.25;
    tickEnd = floor(climVals(2) / 0.25) * 0.25;

    if (~isfinite(tickStart) || ~isfinite(tickEnd) || tickStart > tickEnd)
        ticks = climVals(:).';
        return;
    end

    ticks = tickStart:0.25:tickEnd;
    if (isempty(ticks))
        ticks = climVals(:).';
    end
end
