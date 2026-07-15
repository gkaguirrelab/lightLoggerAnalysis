function figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
% Plot world-camera frame-mean linearity data by AGC target and NDF level.
%
% Syntax:
%   figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
%
% Description:
%   This intentionally simple analysis makes one figure per contrast AGC
%   target. Each NDF level is shown as a pair of plots: frame-mean response
%   with error bars on the left, and the fixed AGC settings on the right.
%   The plotted curves are also saved at the end of the function in
%   camera_linearity_curves_by_NDF.mat.

    arguments
        calibration_metadata
        measurements
        options.verbose logical = true
    end

    configure_default_plot_appearance();

    % Pull the matrix dimensions and metadata that define how the measurement
    % cell array is indexed: NDF, contrast target, settings level, replicate.
    [n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements] = size(measurements);
    settings_levels = double(calibration_metadata.background_scalars(:)');
    ndf_values = double(calibration_metadata.NDFs(:));
    contrast_targets = double(calibration_metadata.contrast_agc_targets(:));
    [sorted_ndfs, ndf_order] = sort(ndf_values);

    % Use one shared response range for every generated plot so the panels are
    % visually comparable across NDFs and contrast targets.
    [y_min, y_max] = frame_mean_response_limits(measurements);
    y_min = min(0, y_min - 25);
    y_max = y_max + 25;
    [error_y_min, error_y_max] = errorbar_response_limits(measurements, n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements);
    y_min = min(y_min, error_y_min - 25);
    y_max = max(y_max, error_y_max + 25);
    if(~isfinite(y_max))
        y_max = 300;
    end
    if(~isfinite(y_min))
        y_min = 0;
    end

    figureHandles = cell(3 * n_contrastTargets, 1);
    figure_idx = 0;
    curves_by_NDF = initialize_curve_export(settings_levels, sorted_ndfs, ndf_order, contrast_targets);

    for cc = 1:n_contrastTargets
        % Convert the raw replicate measurements into the plotted mean curve
        % and +/- 1 SD error bars for the current contrast target.
        frame_mean_by_measurement = extract_frame_mean_by_measurement(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        rgb_mean_by_measurement = extract_rgb_mean_by_measurement(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        frame_mean_by_NDF = mean(frame_mean_by_measurement, 3, 'omitnan');
        rgb_mean_by_NDF = mean(rgb_mean_by_measurement, 4, 'omitnan');
        frame_std_by_NDF = std(frame_mean_by_measurement, 0, 3, 'omitnan');
        agc_settings_by_NDF = extract_agc_settings_by_NDF(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        b_gamma_fit = fit_first_available_b_channel_gamma(settings_levels, rgb_mean_by_NDF, sorted_ndfs, ndf_order);
        curves_by_NDF.contrast_curves{cc} = build_curve_export(contrast_targets(cc), settings_levels, ...
            sorted_ndfs, ndf_order, frame_mean_by_NDF, frame_std_by_NDF, ...
            frame_mean_by_measurement, agc_settings_by_NDF);

        % Each NDF gets a three-column group: two columns for the response
        % curve and one narrower column for the AGC settings bars.
        n_groups_per_row = min(3, n_NDFs);
        n_group_columns = 3;
        n_tile_columns = n_groups_per_row * n_group_columns;
        n_tile_rows = ceil(n_NDFs / n_groups_per_row);

        figHandle = figure('Name', sprintf('Camera_Frame_Mean_Raw_AGC_Target_%0.3g', contrast_targets(cc)));
        set(figHandle, 'Color', 'w');
        set(figHandle, 'Units', 'pixels', 'Position', scaled_figure_position(n_tile_columns, n_tile_rows));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = figHandle;

        tiles = tiledlayout(figHandle, n_tile_rows, n_tile_columns, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tiles, sprintf('World Camera Frame Mean Raw | AGC Target %.3g', contrast_targets(cc)), ...
            'FontWeight', 'Bold', 'FontSize', 18);

        for oo = 1:numel(ndf_order)
            % Plot in sorted NDF order, but keep the original matrix index in
            % each title so the plotted data can be traced back to the input.
            nn = ndf_order(oo);
            tile_row = floor((oo - 1) / n_groups_per_row);
            tile_group = mod(oo - 1, n_groups_per_row);
            first_tile = tile_row * n_tile_columns + tile_group * n_group_columns + 1;

            responseAxes = nexttile(tiles, first_tile, [1, 2]);
            rgb_response = reshape(rgb_mean_by_NDF(nn, :, :), [n_settingsLevels, 3]);
            plot_raw_response(responseAxes, settings_levels, frame_mean_by_NDF(nn, :), ...
                frame_std_by_NDF(nn, :), rgb_response, sorted_ndfs(oo), nn, [y_min, y_max]);

            settingsAxes = nexttile(tiles, first_tile + 2);
            plot_agc_settings(settingsAxes, agc_settings_by_NDF(nn, :), sorted_ndfs(oo), nn);
        end

        if(options.verbose)
            fprintf('World Camera Frame Mean | AGC target %.3g | plotted %d NDF levels with y max %.3g.\n', ...
                contrast_targets(cc), n_NDFs, y_max);
        end

        drawnow;

        bGammaFig = figure('Name', sprintf('B_Channel_Raw_vs_Shared_Gamma_Fit_AGC_Target_%0.3g', contrast_targets(cc)));
        set(bGammaFig, 'Color', 'w');
        set(bGammaFig, 'Units', 'pixels', 'Position', scaled_figure_position(n_tile_columns + n_groups_per_row, n_tile_rows));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = bGammaFig;

        bGammaTiles = tiledlayout(bGammaFig, n_tile_rows, n_tile_columns, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(bGammaTiles, sprintf('B Channel Raw vs Shared Gamma Fit | AGC Target %.3g', contrast_targets(cc)), ...
            'FontWeight', 'Bold', 'FontSize', 18);

        for oo = 1:numel(ndf_order)
            nn = ndf_order(oo);
            tile_row = floor((oo - 1) / n_groups_per_row);
            tile_group = mod(oo - 1, n_groups_per_row);
            first_tile = tile_row * n_tile_columns + tile_group * n_group_columns + 1;

            bFitAxes = nexttile(bGammaTiles, first_tile, [1, 3]);
            b_response = reshape(rgb_mean_by_NDF(nn, :, 3), [1, n_settingsLevels]);
            plot_b_channel_raw_vs_shared_gamma_fit(bFitAxes, settings_levels, b_response, ...
                b_gamma_fit, sorted_ndfs(oo), nn, oo == 1);
        end

        drawnow;

        % Add one separate figure where the NDF response curves for this AGC
        % target are overlaid on the same axes.
        superimposedFig = figure('Name', sprintf('Curves_Superimposed_AGC_Target_%0.3g', contrast_targets(cc)));
        set(superimposedFig, 'Color', 'w');
        set(superimposedFig, 'Units', 'pixels', 'Position', scaled_figure_position(3, 1));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = superimposedFig;
        superimposedTiles = tiledlayout(superimposedFig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        superimposedAxes = nexttile(superimposedTiles);
        plot_curves_superimposed(superimposedAxes, settings_levels, frame_mean_by_NDF, ...
            sorted_ndfs, ndf_order, contrast_targets(cc), [y_min, y_max]);
        drawnow;
    end

    % Save the exact curve data used by the plots, grouped first by contrast
    % target and then by sorted NDF.
    curves_mat_filename = 'camera_linearity_curves_by_NDF.mat';
    save(curves_mat_filename, 'curves_by_NDF');

    if(options.verbose)
        fprintf('Saved camera linearity curves by NDF to %s.\n', curves_mat_filename);
    end
end

function curves_by_NDF = initialize_curve_export(settings_levels, sorted_ndfs, ndf_order, contrast_targets)
% Create the top-level export structure shared by all contrast targets.

    curves_by_NDF = struct();
    curves_by_NDF.created = datetime("now");
    curves_by_NDF.settings_levels = settings_levels;
    curves_by_NDF.sorted_NDFs = sorted_ndfs;
    curves_by_NDF.matrix_indices = ndf_order;
    curves_by_NDF.contrast_targets = contrast_targets;
    curves_by_NDF.agc_setting_names = world_agc_metadata_cols();
    curves_by_NDF.contrast_curves = cell(numel(contrast_targets), 1);
end

function contrast_curve = build_curve_export(contrast_target, settings_levels, sorted_ndfs, ndf_order, ...
    frame_mean_by_NDF, frame_std_by_NDF, frame_mean_by_measurement, agc_settings_by_NDF)
% Package the plotted curves for one contrast target into sorted NDF records.

    agc_setting_names = world_agc_metadata_cols();
    empty_ndf_curve = struct( ...
        'NDF', nan, ...
        'matrix_index', nan, ...
        'settings_levels', [], ...
        'frame_mean', [], ...
        'frame_std', [], ...
        'frame_mean_by_measurement', [], ...
        'agc_settings', [], ...
        'agc_setting_names', {agc_setting_names});
    ndf_curves = repmat(empty_ndf_curve, numel(ndf_order), 1);

    for oo = 1:numel(ndf_order)
        nn = ndf_order(oo);

        ndf_curves(oo).NDF = sorted_ndfs(oo);
        ndf_curves(oo).matrix_index = nn;
        ndf_curves(oo).settings_levels = settings_levels;
        ndf_curves(oo).frame_mean = frame_mean_by_NDF(nn, :);
        ndf_curves(oo).frame_std = frame_std_by_NDF(nn, :);
        ndf_curves(oo).frame_mean_by_measurement = reshape(frame_mean_by_measurement(nn, :, :), ...
            numel(settings_levels), size(frame_mean_by_measurement, 3));
        ndf_curves(oo).agc_settings = agc_settings_by_NDF(nn, :);
        ndf_curves(oo).agc_setting_names = agc_setting_names;
    end

    contrast_curve = struct();
    contrast_curve.contrast_target = contrast_target;
    contrast_curve.settings_levels = settings_levels;
    contrast_curve.sorted_NDFs = sorted_ndfs;
    contrast_curve.matrix_indices = ndf_order;
    contrast_curve.sorted_frame_mean_by_NDF = frame_mean_by_NDF(ndf_order, :);
    contrast_curve.sorted_frame_std_by_NDF = frame_std_by_NDF(ndf_order, :);
    contrast_curve.NDF_curves = ndf_curves;
end

function position = scaled_figure_position(n_tile_columns, n_tile_rows)
% Scale the figure canvas from fixed tile dimensions to keep plot aspect stable.

    tile_width_px = 170;
    tile_height_px = 300;
    width = n_tile_columns * tile_width_px;
    height = n_tile_rows * tile_height_px;

    screen_size = get(groot, 'ScreenSize');
    left = max(20, screen_size(1) + (screen_size(3) - width) / 2);
    bottom = max(50, screen_size(2) + (screen_size(4) - height) / 2);
    position = [left, bottom, width, height];
end

function agc_settings_by_NDF = extract_agc_settings_by_NDF(measurements, contrast_idx, n_NDFs, n_settingsLevels, n_measurements)
% Read one world AGC metadata row per NDF level.

    agc_settings_by_NDF = nan(n_NDFs, numel(world_agc_metadata_cols()));

    for nn = 1:n_NDFs
        % The AGC settings should be fixed within an NDF recording, so use the
        % first measurement that has a complete settings struct.
        for ss = 1:n_settingsLevels
            for mm = 1:n_measurements
                measurement = measurements{nn, contrast_idx, ss, mm};

                if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'settings'))
                    continue;
                end

                agc_settings_by_NDF(nn, :) = world_settings_triplet(measurement.W.settings);
                break;
            end

            if(all(isfinite(agc_settings_by_NDF(nn, :))))
                break;
            end
        end
    end
end

function settings_values = world_settings_triplet(settings)
% Convert one world-camera settings struct to the expected AGC metadata values.

    metadata_cols = world_agc_metadata_cols();
    settings_values = nan(1, numel(metadata_cols));

    for cc = 1:numel(metadata_cols)
        settings_values(cc) = setting_mean(settings, metadata_cols{cc});
    end
end

function metadata_cols = world_agc_metadata_cols()
% Canonical world AGC metadata fields.

    metadata_cols = {'cameraAgain', 'AGCDgain', 'cameraExposure', 'AGCAgain', 'AGCExposure'};
end

function value = setting_mean(settings, setting_name)
% Return the finite mean of one settings field, or NaN if unavailable.

    value = nan;

    if(~isfield(settings, setting_name))
        return;
    end

    setting_values = double(settings.(setting_name)(:));
    setting_values = setting_values(isfinite(setting_values));

    if(~isempty(setting_values))
        value = mean(setting_values, 'omitnan');
    end
end

function [y_min, y_max] = errorbar_response_limits(measurements, n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements)
% Find the plotted mean +/- SD range across every contrast target and NDF.

    y_min = inf;
    y_max = -inf;

    for cc = 1:n_contrastTargets
        frame_mean_by_measurement = extract_frame_mean_by_measurement(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        frame_mean_by_NDF = mean(frame_mean_by_measurement, 3, 'omitnan');
        frame_std_by_NDF = std(frame_mean_by_measurement, 0, 3, 'omitnan');

        y_min = min(y_min, min(frame_mean_by_NDF(:) - frame_std_by_NDF(:), [], 'omitnan'));
        y_max = max(y_max, max(frame_mean_by_NDF(:) + frame_std_by_NDF(:), [], 'omitnan'));
    end
end

function [y_min, y_max] = frame_mean_response_limits(measurements)
% Find the response range across every stored world-camera channel.

    y_min = inf;
    y_max = -inf;

    for ii = 1:numel(measurements)
        response_values = response_values_from_measurement(measurements{ii});

        if(any(isfinite(response_values)))
            y_min = min(y_min, min(response_values, [], 'omitnan'));
            y_max = max(y_max, max(response_values, [], 'omitnan'));
        end
    end
end

function configure_default_plot_appearance()
% Set plot defaults used by this analysis.

    set(groot, 'DefaultAxesFontSize', 14);
    set(groot, 'DefaultTextFontSize', 14);
    set(groot, 'DefaultAxesColor', 'w');
    set(groot, 'DefaultFigureColor', 'w');
end

function frame_mean_by_measurement = extract_frame_mean_by_measurement(measurements, contrast_idx, n_NDFs, n_settingsLevels, n_measurements)
% Average frame-mean camera values across frames, keeping repeated measurements separate.

    frame_mean_by_measurement = nan(n_NDFs, n_settingsLevels, n_measurements);

    for nn = 1:n_NDFs
        for ss = 1:n_settingsLevels
            for mm = 1:n_measurements
                measurement = measurements{nn, contrast_idx, ss, mm};
                frame_mean_by_measurement(nn, ss, mm) = frame_mean_from_measurement(measurement);
            end
        end
    end
end

function rgb_mean_by_measurement = extract_rgb_mean_by_measurement(measurements, contrast_idx, n_NDFs, n_settingsLevels, n_measurements)
% Average R/G/B camera values across frames, keeping repeated measurements separate.

    rgb_mean_by_measurement = nan(n_NDFs, n_settingsLevels, 3, n_measurements);

    for nn = 1:n_NDFs
        for ss = 1:n_settingsLevels
            for mm = 1:n_measurements
                measurement = measurements{nn, contrast_idx, ss, mm};
                rgb_mean_by_measurement(nn, ss, :, mm) = reshape(rgb_mean_from_measurement(measurement), [1, 1, 3, 1]);
            end
        end
    end
end

function frame_mean = frame_mean_from_measurement(measurement)
% Read the frame-mean column from one world-camera measurement.

    frame_mean = nan;

    if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'v'))
        return;
    end

    world_values = double(measurement.W.v);

    if(isempty(world_values))
        return;
    end

    validateattributes(world_values, {'numeric'}, {'2d'}, mfilename, 'measurement.W.v');

    if(size(world_values, 2) >= 4)
        frame_values = world_values(:, 4);
    else
        frame_values = mean(world_values, 2, 'omitnan');
    end

    frame_mean = mean(frame_values, 'omitnan');
end

function rgb_mean = rgb_mean_from_measurement(measurement)
% Read the R/G/B mean columns from one world-camera measurement.

    rgb_mean = nan(1, 3);

    if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'v'))
        return;
    end

    world_values = double(measurement.W.v);

    if(isempty(world_values) || size(world_values, 2) < 3)
        return;
    end

    validateattributes(world_values, {'numeric'}, {'2d'}, mfilename, 'measurement.W.v');
    rgb_mean = mean(world_values(:, 1:3), 1, 'omitnan');
end

function response_values = response_values_from_measurement(measurement)
% Read available R/G/B/frame mean values from one world-camera measurement.

    response_values = nan(1, 4);

    if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'v'))
        return;
    end

    world_values = double(measurement.W.v);

    if(isempty(world_values))
        return;
    end

    validateattributes(world_values, {'numeric'}, {'2d'}, mfilename, 'measurement.W.v');
    n_response_columns = min(size(world_values, 2), numel(response_values));
    response_values(1:n_response_columns) = mean(world_values(:, 1:n_response_columns), 1, 'omitnan');
end

function gamma_fit = fit_first_available_b_channel_gamma(settings_levels, rgb_mean_by_NDF, sorted_ndfs, ndf_order)
% Fit a shared gamma exponent from the first valid sorted-NDF B curve.

    [gamma_black_level, gamma_white_level] = gamma_response_anchors();
    gamma_fit = struct( ...
        'gamma', nan, ...
        'black_level', gamma_black_level, ...
        'white_level', gamma_white_level, ...
        'source_NDF', nan, ...
        'source_matrix_index', nan, ...
        'source_channel', 'B');

    for oo = 1:numel(ndf_order)
        nn = ndf_order(oo);
        b_response = reshape(rgb_mean_by_NDF(nn, :, 3), [1, size(rgb_mean_by_NDF, 2)]);
        [normalized_settings, normalized_response] = normalize_gamma_source_curve(settings_levels, b_response, ...
            gamma_black_level, gamma_white_level);

        if(numel(normalized_settings) < 2)
            continue;
        end

        gamma_fit.gamma = fit_gamma_exponent(normalized_settings, normalized_response);
        gamma_fit.source_NDF = sorted_ndfs(oo);
        gamma_fit.source_matrix_index = nn;
        return;
    end
end

function [black_level, white_level] = gamma_response_anchors()
% Fixed B-channel response anchors used for gamma fitting.

    black_level = 16;
    white_level = 254;
end

function [normalized_settings, normalized_response] = normalize_gamma_source_curve(settings_levels, response, black_level, white_level)
% Normalize one response curve using fixed black/white camera-count anchors.

    settings_levels = settings_levels(:);
    response = response(:);
    valid = isfinite(settings_levels) & isfinite(response) & settings_levels >= 0 & settings_levels <= 1;
    settings_levels = settings_levels(valid);
    response = response(valid);

    normalized_settings = [];
    normalized_response = [];

    if(numel(settings_levels) < 2 || white_level <= black_level)
        return;
    end

    [settings_levels, sort_order] = sort(settings_levels);
    response = response(sort_order);

    normalized_settings = settings_levels;
    normalized_response = (response - black_level) ./ (white_level - black_level);
    normalized_response = min(max(normalized_response, 0), 1);
end

function gamma_value = fit_gamma_exponent(normalized_settings, normalized_response)
% Fit normalized_response = normalized_settings^gamma using positive gamma.

    valid = isfinite(normalized_settings) & isfinite(normalized_response) & ...
        normalized_settings >= 0 & normalized_settings <= 1;
    normalized_settings = normalized_settings(valid);
    normalized_response = normalized_response(valid);

    if(numel(normalized_settings) < 2)
        gamma_value = 1;
        return;
    end

    fit_objective = @(log_gamma) sum((normalized_settings .^ exp(log_gamma) - normalized_response) .^ 2, 'omitnan');
    gamma_value = exp(fminsearch(fit_objective, 0, optimset('Display', 'off')));
end

function plot_b_channel_raw_vs_shared_gamma_fit(ax, settings_levels, b_response, gamma_fit, ndf_value, source_idx, show_legend)
% Plot one NDF B-channel raw curve against the shared source gamma fit.

    settings_levels = settings_levels(:);
    b_response = b_response(:);
    valid = isfinite(settings_levels) & isfinite(b_response);

    hold(ax, 'on');

    if(any(valid))
        plot(ax, settings_levels(valid), b_response(valid), '-o', ...
            'Color', [0, 0, 0], ...
            'MarkerFaceColor', [0, 0, 0], ...
            'MarkerEdgeColor', [0, 0, 0], ...
            'MarkerSize', 5, ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Raw B channel');
    end

    fitted_response = evaluate_shared_b_gamma_fit(settings_levels, b_response, gamma_fit);
    fit_valid = isfinite(settings_levels) & isfinite(fitted_response);

    if(any(fit_valid))
        plot(ax, settings_levels(fit_valid), fitted_response(fit_valid), '-', ...
            'Color', [0.85, 0.10, 0.10], ...
            'LineWidth', 2.0, ...
            'DisplayName', sprintf('Fit from %s channel | source NDF %.3g idx %d | gamma %.3g | %.0f-%.0f counts', ...
                gamma_fit.source_channel, gamma_fit.source_NDF, gamma_fit.source_matrix_index, gamma_fit.gamma, ...
                gamma_fit.black_level, gamma_fit.white_level));
    end

    title(ax, sprintf('NDF %.3g | matrix idx %d', ndf_value, source_idx), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings Level');
    ylabel(ax, 'B Mean');
    format_b_gamma_axes(ax, b_response, fitted_response);
    grid(ax, 'on');
    if(show_legend)
        legend(ax, 'Location', 'eastoutside', 'FontSize', 8, 'Box', 'off');
    end
    hold(ax, 'off');
end

function fitted_response = evaluate_shared_b_gamma_fit(settings_levels, ~, gamma_fit)
% Apply the fixed-anchor shared gamma curve in camera-count space.

    fitted_response = nan(size(settings_levels));

    if(~isfinite(gamma_fit.gamma) || ~isfinite(gamma_fit.black_level) || ...
            ~isfinite(gamma_fit.white_level) || gamma_fit.white_level <= gamma_fit.black_level)
        return;
    end

    normalized_settings = min(max(settings_levels, 0), 1);
    fitted_response = gamma_fit.black_level + (gamma_fit.white_level - gamma_fit.black_level) .* ...
        (normalized_settings .^ gamma_fit.gamma);
end

function format_b_gamma_axes(ax, b_response, fitted_response)
% Give B raw-vs-fit axes a compact response range.

    plotted_values = [b_response(:); fitted_response(:)];
    plotted_values = plotted_values(isfinite(plotted_values));

    if(isempty(plotted_values))
        return;
    end

    y_min = min(plotted_values);
    y_max = max(plotted_values);

    if(y_max <= y_min)
        y_padding = max(abs(y_max), 1) * 0.05;
    else
        y_padding = 0.08 * (y_max - y_min);
    end

    ylim(ax, [y_min - y_padding, y_max + y_padding]);
end

function plot_raw_response(ax, settings_levels, mean_response, error_response, rgb_response, ndf_value, source_idx, y_limits)
% Plot one NDF row and mark the raw minimum and maximum points.

    hold(ax, 'on');
    rgb_colors = [
        0.85, 0.10, 0.10
        0.10, 0.60, 0.20
        0.10, 0.25, 0.85
        ];
    rgb_labels = {'R mean', 'G mean', 'B mean'};

    for cc = 1:size(rgb_response, 2)
        plot(ax, settings_levels, rgb_response(:, cc), '-o', ...
            'Color', rgb_colors(cc, :), ...
            'MarkerFaceColor', rgb_colors(cc, :), ...
            'MarkerEdgeColor', rgb_colors(cc, :), ...
            'MarkerSize', 4, ...
            'LineWidth', 1.3, ...
            'DisplayName', rgb_labels{cc});
    end

    errorbar(ax, settings_levels, mean_response, error_response, '-o', ...
        'Color', [0, 0, 0], ...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 5, ...
        'LineWidth', 2.0, ...
        'CapSize', 6, ...
        'DisplayName', 'Mean +/- 1 SD');

    title(ax, sprintf('NDF %.3g | matrix idx %d', ndf_value, source_idx), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings Level');
    ylabel(ax, 'Frame Mean');
    ylim(ax, y_limits);
    label_min_max_lines(ax, mean_response);
    label_frame_max_rgb_values(ax, settings_levels, mean_response, rgb_response);
    grid(ax, 'on');
    hold(ax, 'off');
end

function plot_curves_superimposed(ax, settings_levels, frame_mean_by_NDF, sorted_ndfs, ndf_order, contrast_target, y_limits)
% Plot all NDF frame-mean curves on one shared axes.

    hold(ax, 'on');
    curve_colors = lines(numel(ndf_order));

    for oo = 1:numel(ndf_order)
        nn = ndf_order(oo);
        plot(ax, settings_levels, frame_mean_by_NDF(nn, :), '-o', ...
            'Color', curve_colors(oo, :), ...
            'MarkerFaceColor', curve_colors(oo, :), ...
            'MarkerEdgeColor', curve_colors(oo, :), ...
            'MarkerSize', 5, ...
            'LineWidth', 1.8, ...
            'DisplayName', sprintf('NDF %.3g | idx %d', sorted_ndfs(oo), nn));
    end

    label_min_max_lines(ax, frame_mean_by_NDF(:));
    title(ax, sprintf('Curves Superimposed | AGC Target %.3g', contrast_target), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings Level');
    ylabel(ax, 'Frame Mean');
    ylim(ax, y_limits);
    grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside', 'FontSize', 8, 'Box', 'off');
    hold(ax, 'off');
end

function plot_agc_settings(ax, agc_settings, ndf_value, source_idx)
% Plot the AGC settings for one NDF level with hard-coded maxima.

    % Exposure values are log-scaled so they can be compared on the same axes
    % as analog and digital gain values.
    setting_labels = {'cameraAgain', 'AGCDgain', 'log10 cameraExposure', 'AGCAgain', 'log10 AGCExposure'};
    agc_settings = [agc_settings(1), agc_settings(2), log10_positive(agc_settings(3)), ...
                    agc_settings(4), log10_positive(agc_settings(5))];
    setting_maxima = [10.333, 10, log10(8333), 10.333, log10(8333)];
    max_line_values = [10.333, 10, log10(8333)];
    max_line_labels = {'Again max', 'Dgain max', 'Exposure max'};
    setting_colors = lines(numel(setting_labels));
    max_line_colors = [
        0.00, 0.35, 0.75
        0.65, 0.20, 0.75
        0.95, 0.55, 0.10
        ];

    hold(ax, 'on');

    % Draw the measured settings as bars, then overlay shared max-reference
    % lines for gain, digital gain, and exposure.
    bar(ax, 1:numel(agc_settings), agc_settings, ...
        'FaceColor', 'flat', ...
        'CData', setting_colors, ...
        'DisplayName', 'AGC setting');

    for ss = 1:numel(max_line_values)
        yline(ax, max_line_values(ss), '--', ...
            'Color', max_line_colors(ss, :), ...
            'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end

    title(ax, sprintf('NDF %.3g AGC | matrix idx %d', ndf_value, source_idx), 'FontWeight', 'Bold');
    ylabel(ax, 'Setting Value');
    xticks(ax, 1:numel(setting_labels));
    xticklabels(ax, setting_labels);
    xtickangle(ax, 25);
    ylim(ax, [0, max([agc_settings(:); setting_maxima(:)], [], 'omitnan') * 1.10]);
    xlim(ax, [0.5, numel(setting_labels) + 0.5]);
    label_agc_bar_values(ax, agc_settings);
    label_max_lines_outside_right(ax, max_line_labels, max_line_values, max_line_colors);
    ax.FontSize = 11;
    grid(ax, 'on');
    hold(ax, 'off');
end

function label_agc_bar_values(ax, agc_settings)
% Put each plotted AGC value inside its bar.

    y_limits = ylim(ax);
    y_span = diff(y_limits);

    for ss = 1:numel(agc_settings)
        value = agc_settings(ss);

        if(~isfinite(value))
            continue;
        end

        if(value > 0.12 * y_span)
            y_position = value / 2;
        else
            y_position = value + 0.03 * y_span;
        end

        text(ax, ss, y_position, sprintf('%.3g', value), ...
            'Color', [0, 0, 0], ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Rotation', 90);
    end
end

function value = log10_positive(value)
% Log-transform positive values while leaving invalid values as NaN.

    if(~isfinite(value) || value <= 0)
        value = nan;
        return;
    end

    value = log10(value);
end

function label_max_lines_outside_right(ax, setting_labels, setting_maxima, setting_colors)
% Put compact max labels just outside the right side of the AGC axes.

    x_limits = xlim(ax);
    x_label = x_limits(2) + 0.05 * diff(x_limits);

    for ss = 1:numel(setting_maxima)
        if(~isfinite(setting_maxima(ss)))
            continue;
        end

        text(ax, x_label, setting_maxima(ss), setting_labels{ss}, ...
            'Color', setting_colors(ss, :), ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'Clipping', 'off');
    end
end

function label_min_max_lines(ax, raw_response)
% Draw horizontal min/max reference lines and label their values.

    raw_response = raw_response(:);

    valid = isfinite(raw_response);

    if(~any(valid))
        return;
    end

    y = raw_response(valid);
    min_value = min(y);
    max_value = max(y);
    x_limits = xlim(ax);
    x_label = x_limits(1) + 0.02 * diff(x_limits);

    yline(ax, min_value, '--', ...
        'Color', [0.85, 0.10, 0.10], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    yline(ax, max_value, '--', ...
        'Color', [0.15, 0.60, 0.25], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    text(ax, x_label, min_value, sprintf('min %.1f', min_value), ...
        'Color', [0.65, 0.15, 0.15], ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
    text(ax, x_label, max_value, sprintf('max %.1f', max_value), ...
        'Color', [0.05, 0.35, 0.15], ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
end

function label_frame_max_rgb_values(ax, settings_levels, mean_response, rgb_response)
% Label the RGB channel values at the frame-mean maximum.

    settings_levels = settings_levels(:);
    mean_response = mean_response(:);

    if(size(rgb_response, 2) ~= 3)
        return;
    end

    valid = isfinite(settings_levels) & isfinite(mean_response);

    if(~any(valid))
        return;
    end

    valid_indices = find(valid);
    [max_value, local_max_idx] = max(mean_response(valid_indices));
    max_idx = valid_indices(local_max_idx);
    rgb_values = rgb_response(max_idx, :);

    if(~all(isfinite(rgb_values)))
        return;
    end

    y_limits = ylim(ax);
    y_span = diff(y_limits);
    y_offset = 0.04 * y_span;
    if(max_value + 0.18 * y_span < y_limits(2))
        y_position = max_value + y_offset;
        vertical_alignment = 'bottom';
    else
        y_position = max_value - y_offset;
        vertical_alignment = 'top';
    end

    plot(ax, settings_levels(max_idx), max_value, 'p', ...
        'Color', [0.05, 0.35, 0.15], ...
        'MarkerFaceColor', [0.05, 0.35, 0.15], ...
        'MarkerSize', 8, ...
        'HandleVisibility', 'off');
    text(ax, settings_levels(max_idx), y_position, ...
        sprintf('max %.1f\nR %.1f  G %.1f  B %.1f', max_value, rgb_values(1), rgb_values(2), rgb_values(3)), ...
        'Color', [0.05, 0.25, 0.10], ...
        'VerticalAlignment', vertical_alignment, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2, ...
        'Clipping', 'on');
end
