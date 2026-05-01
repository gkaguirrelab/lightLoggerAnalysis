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
    fig = figure('Color', 'w', 'Position', [100 100 900 650]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'off');

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


            [center_freq, center_spd] = iCleanSpdForPlot(frq, region_averages.center);
            if (~isempty(center_freq))
                center_handle = loglog(ax, center_freq, center_spd, 'Color', center_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                legend_handles(end+1) = center_handle; 
                legend_labels{end+1} = sprintf('%s %s center', color_mode, projection_type); %
            end

            [periphery_freq, periphery_spd] = iCleanSpdForPlot(frq, region_averages.periphery);
            if (~isempty(periphery_freq))
                periphery_handle = loglog(ax, periphery_freq, periphery_spd, 'Color', periphery_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                legend_handles(end+1) = periphery_handle; 
                legend_labels{end+1} = sprintf('%s %s periphery', color_mode, projection_type);
            end
        end
    end

    % Label the Plot
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Spectral power density (contrast^2/Hz)');
    title(ax, sprintf('%s | %s', subject_id, activity_name), 'Interpreter', 'none');

    if (~isempty(legend_handles))
        legend(ax, legend_handles, legend_labels, ...
            'Interpreter', 'none', 'Location', 'best');
    end

    % Return the figure handle
    figure_handle = fig;
end


function [frq, region_averages] = iLoadSpdRegionAverages(spd_path, projection_type)
    frq = [];
    region_averages = [];

    activity_data = load(spd_path).activityData;
    activity_names = fieldnames(activity_data);
    if (isempty(activity_names))
        return;
    end

    activity_name = activity_names{1};
    frq = activity_data.(activity_name).frq(:);
    spd_by_region = activity_data.(activity_name).spdByRegion;
    exponent_map = activity_data.(activity_name).exponentMap;

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
