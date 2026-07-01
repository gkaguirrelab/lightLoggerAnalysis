function combineSPDs(spds, output_dir, options)
% Generate combined SPD plots and best-fit exports for each subject and activity
%
% Syntax:
%   combineSPDs(spds, output_dir)
%   combineSPDs(spds, output_dir, options)
%
% Description:
%   Iterates over a nested structure of SPD file paths organized by
%   subject, activity, color mode, and projection type. For each
%   subject-activity pair, the function loads the stored SPD results,
%   generates a combined SPD-by-region figure with center and periphery
%   curves across all color modes and projection types, and exports
%   exponent and variance spatial maps as PDF files. Best-fit line
%   coefficients (slope and intercept in log-log space) are saved as
%   .mat files for downstream analysis.
%
% Inputs:
%   spds                  - Struct or char/string filepath. Nested
%                           structure of SPD file paths organized as:
%                             spds.(subject).(activity).(color_mode).(projection_type)
%                           If a filepath is provided, the struct is
%                           loaded from the .mat file at that path.
%   output_dir            - Char/string. Directory where figures and
%                           best-fit .mat files are saved. Created
%                           automatically if it does not exist.
%
% Optional key/value pairs:
%   overwrite_existing    - Logical. If true, regenerate outputs even
%                           when they already exist. Defaults to false.
%   verbose               - Logical. If true, print progress information
%                           during processing. Defaults to true.
%
% Outputs:
%   none
%
% Examples:
%{
    spds = load("spd_paths.mat").spds;
    combineSPDs(spds, "/path/to/output", "overwrite_existing", true);
%}

    arguments
        spds
        output_dir
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
    if (~isfolder(output_dir))
        mkdir(output_dir);
    end

    % Iterate over subjects and activities.
    subjects = fieldnames(spds);
    num_subjects = numel(subjects);
    for ss = 1:num_subjects 
        subject_id = subjects{ss};
        subject_struct = spds.(subject_id);
        subject_output_dir = fullfile(output_dir, subject_id);
        if(~isfolder(subject_output_dir))
            mkdir(subject_output_dir); 
        end     

        activity_names = fieldnames(subject_struct);
        num_activities = numel(activity_names); 
        for aa = 1:num_activities
            activity_name = activity_names{aa};
            activity_struct = subject_struct.(activity_name);
            activity_output_dir = fullfile(subject_output_dir, activity_name);
            if(~isfolder(activity_output_dir))
                mkdir(activity_output_dir); 
            end 

            % Generate the output filepath for the SPD graph
            output_filepath = fullfile(activity_output_dir, sprintf('%s_%s_spdByRegion.pdf', subject_id, activity_name));

            % If we do not want to overwrite existing data, just skip
            if (isfile(output_filepath) && ~options.overwrite_existing)
                continue;
            end

            % Print out verbose information if we want it
            if (options.verbose)
                fprintf('Generating SPD plots for subject: %d / %d | activity: %d / %d\n', ss, num_subjects, aa, num_activities);
            end

            % Plot the SPDs and calculate their best fit lines as well
            [figure_handle, activity_best_fit_lines] = localPlotSPD(activity_struct, subject_id, activity_name);
            
            % Save the figure 
            exportgraphics(figure_handle, output_filepath, 'ContentType', 'vector');
            close(figure_handle);

            % Also plot exponent and variance maps for each color mode
            % using the same map-plotting style as plotSPDs_copy. The
            % helper returns figure handles only; export them here.
            map_handles = localPlotMaps(activity_struct, activity_name);
            activity_color_modes = fieldnames(map_handles);
            for cc = 1:numel(activity_color_modes)
                color_mode = activity_color_modes{cc};

                exponent_output_filepath = fullfile(activity_output_dir, sprintf('%s_%s_%s_exponentMap.pdf', subject_id, activity_name, color_mode));
                variance_output_filepath = fullfile(activity_output_dir, sprintf('%s_%s_%s_varianceMap.pdf', subject_id, activity_name, color_mode));

                exportgraphics(map_handles.(color_mode).exponentMap, exponent_output_filepath, 'ContentType', 'vector');
                exportgraphics(map_handles.(color_mode).varianceMap, variance_output_filepath, 'ContentType', 'vector');
                close(map_handles.(color_mode).exponentMap);
                close(map_handles.(color_mode).varianceMap);
            end

            % Output the best fit line information for each color mode and projection type of this activity 
            % as .mat files 
            activity_color_modes = fieldnames(activity_best_fit_lines); 
            for cc = 1:numel(activity_color_modes)
                color_mode = activity_color_modes{cc}; 
                color_mode_struct = activity_best_fit_lines.(color_mode); 

                projection_types = fieldnames(color_mode_struct); 
                for pp = 1:numel(projection_types)
                    projection_type = projection_types{pp};
                    projection_struct = color_mode_struct.(projection_type); 
                    bestFit = projection_struct; 

                    % Output the .mat file here. Include both the color
                    % mode and projection type in the filename so each
                    % saved fit remains distinct.
                    best_fit_output_path = fullfile(activity_output_dir, sprintf("%s_%s_bestFit.mat", color_mode, projection_type)); 
                    save(best_fit_output_path, 'bestFit'); 
                end 
            end     
        
        end
    end
end

% Local Function to Handle the SPD plotting
function [figure_handle, best_fit_lines] = localPlotSPD(activity_struct, subject_id, activity_name)
% Internal helper to local plot spd.
%
% Syntax:
%   figure_handle, best_fit_lines = localPlotSPD(activity_struct, subject_id, activity_name)
%
% Description:
%   This local helper function internal helper to local plot spd within its parent workflow.
% Inputs:
%   activity_struct          - Input used by the function.
%   subject_id               - Input used by the function.
%   activity_name            - Input used by the function.
%
% Outputs:
%   figure_handle            - Output produced by the function.
%   best_fit_lines           - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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

    % Initialize containers for the legends of the figures
    legend_handles = gobjects(0);
    legend_labels = {};
    best_fit_lines = struct();

    % Iterate over the color modes
    color_modes = fieldnames(activity_struct);
    has_multiple_color_modes = numel(color_modes) > 1;
    for cc = 1:numel(color_modes)
        % Retireve the name of this color mode
        color_mode = color_modes{cc};
        has_added_compact_color_mode_legend_item = false;
        
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
                if (~has_multiple_color_modes || cc == 1)
                    legend_handles(end+1) = center_handle; 
                    legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'center'); 
                elseif (~has_added_compact_color_mode_legend_item)
                    legend_handles(end+1) = center_handle;
                    legend_labels{end+1} = iFormatColorModeLabel(color_mode);
                    has_added_compact_color_mode_legend_item = true;
                end
                best_fit_lines.(color_mode).(projection_type).center = iComputeBestFitLine(center_freq, center_spd);
            end

            % Plot the periphery SPD with a lighter version
            % of the same color mode color
            [periphery_freq, periphery_spd] = iCleanSpdForPlot(frq, region_averages.periphery);
            if (~isempty(periphery_freq))
                periphery_handle = loglog(ax, periphery_freq, periphery_spd, 'Color', periphery_color, 'LineStyle', line_style, 'LineWidth', 2.0);
                if (~has_multiple_color_modes || cc == 1)
                    legend_handles(end+1) = periphery_handle; 
                    legend_labels{end+1} = iFormatLegendLabel(color_mode, projection_type, 'periphery');
                end
                best_fit_lines.(color_mode).(projection_type).periphery = iComputeBestFitLine(periphery_freq, periphery_spd);
            end
        end
    end

    % Label the Plot
    xlabel(ax, 'Frequency (Hz)', 'FontSize', label_font_size);
    ylabel(ax, 'Spectral power density (contrast^2/Hz)', 'FontSize', label_font_size);
    title(ax, sprintf('%s | %s', subject_id, activity_name), 'Interpreter', 'none', 'FontSize', title_font_size);
    axis(ax, 'square');
    xlim(ax, [1 64]);
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


function best_fit_line = iComputeBestFitLine(frq, spd)
% Internal helper to i compute best fit line.
%
% Syntax:
%   best_fit_line = iComputeBestFitLine(frq, spd)
%
% Description:
%   This local helper function internal helper to i compute best fit line within its parent workflow.
% Inputs:
%   frq                      - Input used by the function.
%   spd                      - Input used by the function.
%
% Outputs:
%   best_fit_line            - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    x = log10(frq(:));
    y = log10(spd(:));
    coefficients = polyfit(x, y, 1);
    y_fit = polyval(coefficients, x);

    best_fit_line = struct();
    best_fit_line.slope = coefficients(1);
    best_fit_line.intercept = coefficients(2);
    best_fit_line.log10_frequency = x;
    best_fit_line.log10_spd_fit = y_fit;
    best_fit_line.frequency = frq(:);
    best_fit_line.spd_fit = 10 .^ y_fit;
end


function [frq, region_averages] = iLoadSpdRegionAverages(spd_path, projection_type)
% Internal helper to i load spd region averages.
%
% Syntax:
%   frq, region_averages = iLoadSpdRegionAverages(spd_path, projection_type)
%
% Description:
%   This local helper function internal helper to i load spd region averages within its parent workflow.
% Inputs:
%   spd_path                 - Path-like input used by the function.
%   projection_type          - Input used by the function.
%
% Outputs:
%   frq                      - Output produced by the function.
%   region_averages          - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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


function map_handles = localPlotMaps(activity_struct, activity_name)
% Internal helper to local plot maps.
%
% Syntax:
%   map_handles = localPlotMaps(activity_struct, activity_name)
%
% Description:
%   This local helper function internal helper to local plot maps within its parent workflow.
% Inputs:
%   activity_struct          - Input used by the function.
%   activity_name            - Input used by the function.
%
% Outputs:
%   map_handles              - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    color_modes = fieldnames(activity_struct);
    map_handles = struct();

    for cc = 1:numel(color_modes)
        color_mode = color_modes{cc};
        color_mode_struct = activity_struct.(color_mode);

        just_projection_data = iLoadProjectionPayload(color_mode_struct.justProjection);
        virtually_foveated_data = iLoadProjectionPayload(color_mode_struct.virtuallyFoveated);

        just_projection_exponent_map = just_projection_data.exponentMap;
        just_projection_variance_map = just_projection_data.varianceMap;
        virtually_foveated_exponent_map = virtually_foveated_data.exponentMap;
        virtually_foveated_variance_map = virtually_foveated_data.varianceMap;

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


function loaded_projection = iLoadProjectionPayload(projection_entry)
% Internal helper to i load projection payload.
%
% Syntax:
%   loaded_projection = iLoadProjectionPayload(projection_entry)
%
% Description:
%   This local helper function internal helper to i load projection payload within its parent workflow.
% Inputs:
%   projection_entry         - Input used by the function.
%
% Outputs:
%   loaded_projection        - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    if (isstring(projection_entry) || ischar(projection_entry))
        loaded_mat = load(projection_entry);
        activity_names = fieldnames(loaded_mat.activityData);
        loaded_projection = loaded_mat.activityData.(activity_names{1});
    else
        loaded_projection = projection_entry;
    end
end


function line_style = iGetProjectionLineStyle(projection_type)
% Internal helper to i get projection line style.
%
% Syntax:
%   line_style = iGetProjectionLineStyle(projection_type)
%
% Description:
%   This local helper function internal helper to i get projection line style within its parent workflow.
% Inputs:
%   projection_type          - Input used by the function.
%
% Outputs:
%   line_style               - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i get color mode color.
%
% Syntax:
%   color = iGetColorModeColor(color_mode)
%
% Description:
%   This local helper function internal helper to i get color mode color within its parent workflow.
% Inputs:
%   color_mode               - Input used by the function.
%
% Outputs:
%   color                    - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i lighten color.
%
% Syntax:
%   lighter_color = iLightenColor(color, amount)
%
% Description:
%   This local helper function internal helper to i lighten color within its parent workflow.
% Inputs:
%   color                    - Input used by the function.
%   amount                   - Input used by the function.
%
% Outputs:
%   lighter_color            - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    lighter_color = color + (1 - color) .* amount;
end


function [frq_clean, spd_clean] = iCleanSpdForPlot(frq, spd)
% Internal helper to i clean spd for plot.
%
% Syntax:
%   frq_clean, spd_clean = iCleanSpdForPlot(frq, spd)
%
% Description:
%   This local helper function internal helper to i clean spd for plot within its parent workflow.
% Inputs:
%   frq                      - Input used by the function.
%   spd                      - Input used by the function.
%
% Outputs:
%   frq_clean                - Output produced by the function.
%   spd_clean                - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i compute region mean.
%
% Syntax:
%   region_mean_spd = iComputeRegionMean(spd_by_region, region_mask)
%
% Description:
%   This local helper function internal helper to i compute region mean within its parent workflow.
% Inputs:
%   spd_by_region            - Input used by the function.
%   region_mask              - Input used by the function.
%
% Outputs:
%   region_mean_spd          - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i plot reference line.
%
% Syntax:
%   reference_handle = iPlotReferenceLine(ax)
%
% Description:
%   This local helper function internal helper to i plot reference line within its parent workflow.
% Inputs:
%   ax                       - Input used by the function.
%
% Outputs:
%   reference_handle         - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    drawnow;
    axis_limits = axis(ax);
    reference_freq = axis_limits(1:2).';
    reference_spd = 10^-2 .* (reference_freq .^ -2);
    reference_handle = loglog(ax, reference_freq, reference_spd, ':', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1.5);
end


function legend_label = iFormatLegendLabel(color_mode, projection_type, region_name)
% Internal helper to i format legend label.
%
% Syntax:
%   legend_label = iFormatLegendLabel(color_mode, projection_type, region_name)
%
% Description:
%   This local helper function internal helper to i format legend label within its parent workflow.
% Inputs:
%   color_mode               - Input used by the function.
%   projection_type          - Input used by the function.
%   region_name              - Input used by the function.
%
% Outputs:
%   legend_label             - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

    color_mode_label = iFormatColorModeLabel(color_mode);
    projection_label = iFormatProjectionLabel(projection_type);
    legend_label = sprintf('%s %s %s', color_mode_label, projection_label, region_name);
end


function color_mode_label = iFormatColorModeLabel(color_mode)
% Internal helper to i format color mode label.
%
% Syntax:
%   color_mode_label = iFormatColorModeLabel(color_mode)
%
% Description:
%   This local helper function internal helper to i format color mode label within its parent workflow.
% Inputs:
%   color_mode               - Input used by the function.
%
% Outputs:
%   color_mode_label         - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i format projection label.
%
% Syntax:
%   projection_label = iFormatProjectionLabel(projection_type)
%
% Description:
%   This local helper function internal helper to i format projection label within its parent workflow.
% Inputs:
%   projection_type          - Input used by the function.
%
% Outputs:
%   projection_label         - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i get formula annotation text.
%
% Syntax:
%   formula_text = iGetFormulaAnnotationText()
%
% Description:
%   This local helper function internal helper to i get formula annotation text within its parent workflow.
% Inputs:
%   None.
%
% Outputs:
%   formula_text             - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
% Internal helper to i get quarter step ticks.
%
% Syntax:
%   ticks = iGetQuarterStepTicks(climVals)
%
% Description:
%   This local helper function internal helper to i get quarter step ticks within its parent workflow.
% Inputs:
%   climVals                 - Input used by the function.
%
% Outputs:
%   ticks                    - Output produced by the function.
%
% Examples:
%{
    % See combineSPDs.m for usage context.
%}

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
