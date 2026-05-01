function combineSPDs(spds, output_path, options)
    arguments
        spds
        output_path
        options.overwrite_existing = false
        options.verbose = true
    end

    % If given a path to the SPDs instead of the struct itself,
    % read it in.
    if (isstring(spds) || ischar(spds))
        spds = load(spds).spds;
    end

    %{
    Input dictionary should be of the form
    {subject: {activity:
                  {color_mode:
                              {projection_type: path_to_mat_spd_results}
                   }
              }
    }
    %}

    % If the output path does not exist, make it
    if (~isfolder(output_path))
        mkdir(output_path);
    end

    % Iterate over subjects and activities.
    subjects = fieldnames(spds);
    for ss = 1:numel(subjects)
        subject_id = subjects{ss};
        subject_struct = spds.(subject_id);

        activity_names = fieldnames(subject_struct);
        for aa = 1:numel(activity_names)
            activity_name = activity_names{aa};
            activity_struct = subject_struct.(activity_name);

            % Generate the output filepath
            output_filepath = fullfile(output_path, sprintf('%s_%s_combinedSPDs.pdf', subject_id, activity_name));

            % If we do not want to overwrite existing data, just skip
            if (isfile(output_filepath) && ~options.overwrite_existing)
                continue;
            end

            % Print out verbose information if we want it
            if (options.verbose)
                fprintf('Combining SPDs for %s | %s\n', subject_id, activity_name);
            end

            % Plot and output the combined SPDs
            figure_handle = localPlotSPD(activity_struct, subject_id, activity_name);
            exportgraphics(figure_handle, output_filepath, 'ContentType', 'vector');
            close(figure_handle);
        end
    end
end

% Local Function to Handle the SPD plotting
function figure_handle = localPlotSPD(activity_struct, subject_id, activity_name)
    % Initialize the figure that we will draw too
    fig = figure('Color', 'w', 'Position', [100 100 1100 650]);
    ax = axes('Parent', fig, 'Position', [0.08 0.14 0.60 0.76]);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'XScale', 'log', 'YScale', 'log');

    % Initialize containers for the legends of the figures
    legend_handles = gobjects(0);
    legend_labels = {};

    % Iterate over the color modes
    color_modes = fieldnames(activity_struct);
    for cc = 1:numel(color_modes)
        % Retireve the name of this color mode
        color_mode = color_modes{cc};
        
        % Retrieve the struct of projection types for this
        % activity 
        projection_struct = activity_struct.(color_mode);
        projection_types = fieldnames(projection_struct);

        % Iterate over the projection types 
        for pp = 1:numel(projection_types)

            % Gather the name of the projection type 
            projection_type = projection_types{pp};

            % Get the path to the SPD
            spd_path = projection_struct.(projection_type);

            if (~isfile(spd_path))
                continue;
            end

            [frq, region_averages] = iLoadSpdRegionAverages(spd_path, projection_type);
            if (isempty(frq) || isempty(region_averages))
                continue;
            end

            % Gather the line-based information for this graph
            line_style = iGetProjectionLineStyle(projection_type);
            base_color = iGetColorModeColor(color_mode);
            center_color = base_color;
            periphery_color = iLightenColor(base_color, 0.55);

            % Plot the center SPD with the saturated version
            % of the color mode color
            [center_freq, center_spd] = iCleanSpdForPlot(frq, region_averages.center);
            if (~isempty(center_freq))
                center_handle = loglog(ax, center_freq, center_spd, 'Color', center_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                legend_handles(end+1) = center_handle; 
                legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'center'); 
            end

            % Plot the periphery SPD with a lighter version
            % of the same color mode color
            [periphery_freq, periphery_spd] = iCleanSpdForPlot(frq, region_averages.periphery);
            if (~isempty(periphery_freq))
                periphery_handle = loglog(ax, periphery_freq, periphery_spd, 'Color', periphery_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                legend_handles(end+1) = periphery_handle; 
                legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'periphery');
            end
        end
    end

    % Label the Plot
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Spectral power density (contrast^2/Hz)');
    title(ax, sprintf('%s | %s', subject_id, activity_name), 'Interpreter', 'none');
    axis(ax, 'square');

    reference_handle = iPlotReferenceLine(ax);
    if (isgraphics(reference_handle))
        legend_handles(end+1) = reference_handle;
        legend_labels{end+1} = 'Reference (10^{-2} f^{-2})';
    end

    if (~isempty(legend_handles))
        legend(ax, legend_handles, legend_labels, ...
            'Interpreter', 'tex', 'Location', 'northeastoutside');
    end

    % Add a separate textbox with the formula definitions so the
    % line legend remains compact and easy to scan. Place this
    % beneath the legend to preserve more horizontal plot space.
    annotation(fig, 'textbox', [0.72 0.08 0.24 0.40], ...
        'String', iGetFormulaAnnotationText(), ...
        'Interpreter', 'tex', ...
        'FitBoxToText', 'off', ...
        'EdgeColor', [0.7 0.7 0.7], ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    % Return the figure handle
    figure_handle = fig;
end


function [frq, region_averages] = iLoadSpdRegionAverages(spd_path, projection_type)
    frq = [];
    region_averages = [];

    % Load in the SPD activity data from the targeted filepath
    activity_data = load(spd_path).activityData;
    activity_names = fieldnames(activity_data);
    if (isempty(activity_names))
        return;
    end

    % Retrieve the activity-level quantities needed to compute
    % the center and periphery regional averages
    activity_name = activity_names{1};
    frq = activity_data.(activity_name).frq(:);
    spd_by_region = activity_data.(activity_name).spdByRegion;
    exponent_map = activity_data.(activity_name).exponentMap;

    % Build the same center/periphery masks used elsewhere
    % in the SPD plotting code
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

    % Compute the mean SPD across the center and periphery
    % regions separately
    region_averages.center = iComputeRegionMean(spd_by_region, center_mask);
    region_averages.periphery = iComputeRegionMean(spd_by_region, periphery_mask);
end


function line_style = iGetProjectionLineStyle(projection_type)
    % Virtually foveated data should be solid and just
    % projection data should be dashed
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
    % Match the requested plotting colors for each
    % post-receptoral / color mode family
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
    % Create a lighter companion color for the periphery
    % traces while preserving the original hue
    lighter_color = color + (1 - color) .* amount;
end


function [frq_clean, spd_clean] = iCleanSpdForPlot(frq, spd)
    % Flatten, align, and remove invalid entries before
    % plotting on log axes
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
    % Compute the mean SPD across all patches belonging
    % to the provided spatial mask
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
    % Plot the 1/f^2 reference line across the visible
    % frequency range of the current axes
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
