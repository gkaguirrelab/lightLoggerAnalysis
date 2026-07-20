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
%   Superimposed NDF curves and an ND 0 / ND 0.4 B-channel gamma-shift comparison
%   are also shown. The plotted curves are saved at the end of the function in
%   camera_linearity_curves_by_NDF.mat, and simple ND 0 / ND 0.4 mean
%   measurements are saved in camera_linearity_ND0_ND0p4_rgb_means.mat.

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

    figureHandles = cell(5 * n_contrastTargets, 1);
    figure_idx = 0;
    curves_by_NDF = initialize_curve_export(settings_levels, sorted_ndfs, ndf_order, contrast_targets);
    nd04_rgb_means = initialize_nd_pair_rgb_export();

    for cc = 1:n_contrastTargets
        % Convert the raw replicate measurements into the plotted mean curve
        % and +/- 1 SD error bars for the current contrast target.
        frame_mean_by_measurement = extract_frame_mean_by_measurement(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        rgb_mean_by_measurement = extract_rgb_mean_by_measurement(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        frame_mean_by_NDF = mean(frame_mean_by_measurement, 3, 'omitnan');
        rgb_mean_by_NDF = mean(rgb_mean_by_measurement, 4, 'omitnan');
        frame_std_by_NDF = std(frame_mean_by_measurement, 0, 3, 'omitnan');
        agc_settings_by_NDF = extract_agc_settings_by_NDF(measurements, cc, n_NDFs, n_settingsLevels, n_measurements);
        curves_by_NDF.contrast_curves{cc} = build_curve_export(contrast_targets(cc), settings_levels, ...
            sorted_ndfs, ndf_order, frame_mean_by_NDF, frame_std_by_NDF, ...
            frame_mean_by_measurement, agc_settings_by_NDF);
        nd04_rgb_means = add_nd_pair_rgb_export_record(nd04_rgb_means, cc, ...
            settings_levels, sorted_ndfs, ndf_order, rgb_mean_by_NDF, frame_mean_by_NDF);

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
            % Plot in sorted NDF order so the panels progress by filter density.
            nn = ndf_order(oo);
            tile_row = floor((oo - 1) / n_groups_per_row);
            tile_group = mod(oo - 1, n_groups_per_row);
            first_tile = tile_row * n_tile_columns + tile_group * n_group_columns + 1;

            responseAxes = nexttile(tiles, first_tile, [1, 2]);
            rgb_response = reshape(rgb_mean_by_NDF(nn, :, :), [n_settingsLevels, 3]);
            plot_raw_response(responseAxes, settings_levels, frame_mean_by_NDF(nn, :), ...
                frame_std_by_NDF(nn, :), rgb_response, sorted_ndfs(oo), [y_min, y_max]);

            settingsAxes = nexttile(tiles, first_tile + 2);
            plot_agc_settings(settingsAxes, agc_settings_by_NDF(nn, :), sorted_ndfs(oo));
        end

        if(options.verbose)
            fprintf('World Camera Frame Mean | AGC target %.3g | plotted %d NDF levels with y max %.3g.\n', ...
                contrast_targets(cc), n_NDFs, y_max);
        end

        drawnow;

        bGammaFig = figure('Name', sprintf('B_Channel_Raw_vs_Per_NDF_Gamma_Fit_AGC_Target_%0.3g', contrast_targets(cc)));
        set(bGammaFig, 'Color', 'w');
        set(bGammaFig, 'Units', 'pixels', 'Position', scaled_figure_position(n_tile_columns + n_groups_per_row, n_tile_rows));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = bGammaFig;

        bGammaTiles = tiledlayout(bGammaFig, n_tile_rows, n_tile_columns, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(bGammaTiles, sprintf('B Channel Raw vs Per-NDF Gamma Fit | AGC Target %.3g', contrast_targets(cc)), ...
            'FontWeight', 'Bold', 'FontSize', 18);

        for oo = 1:numel(ndf_order)
            nn = ndf_order(oo);
            tile_row = floor((oo - 1) / n_groups_per_row);
            tile_group = mod(oo - 1, n_groups_per_row);
            first_tile = tile_row * n_tile_columns + tile_group * n_group_columns + 1;

            bFitAxes = nexttile(bGammaTiles, first_tile, [1, 3]);
            b_response = reshape(rgb_mean_by_NDF(nn, :, 3), [1, n_settingsLevels]);
            b_gamma_fit = fit_b_channel_gamma(settings_levels, b_response, ndf_label(sorted_ndfs(oo)));
            plot_b_channel_raw_vs_gamma_fit(bFitAxes, settings_levels, b_response, ...
                b_gamma_fit, sorted_ndfs(oo), oo == 1);
        end

        drawnow;

        % Add one separate figure where the NDF response curves for this AGC
        % target are overlaid on the same axes.
        superimposedFig = figure('Name', sprintf('Curves_Superimposed_AGC_Target_%0.3g', contrast_targets(cc)));
        set(superimposedFig, 'Color', 'w');
        set(superimposedFig, 'Units', 'pixels', 'Position', scaled_figure_position(6, 2));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = superimposedFig;
        superimposedTiles = tiledlayout(superimposedFig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        superimposedAxes = nexttile(superimposedTiles);
        plot_curves_superimposed(superimposedAxes, settings_levels, frame_mean_by_NDF, ...
            sorted_ndfs, ndf_order, contrast_targets(cc), [y_min, y_max]);
        drawnow;

        settingsShiftFig = figure('Name', sprintf('Gamma_Fitting_Plot_AGC_Target_%0.3g', contrast_targets(cc)));
        set(settingsShiftFig, 'Color', 'w');
        set(settingsShiftFig, 'Units', 'pixels', 'Position', scaled_figure_position(8, 3));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = settingsShiftFig;
        settingsShiftTiles = tiledlayout(settingsShiftFig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(settingsShiftTiles, sprintf('Gamma fitting plot | AGC Target %.3g', contrast_targets(cc)), ...
            'FontWeight', 'Bold', 'FontSize', 18);
        plot_gamma_fitting_shift(settingsShiftTiles, settings_levels, frame_mean_by_NDF, frame_std_by_NDF, ...
            rgb_mean_by_NDF, agc_settings_by_NDF, ...
            sorted_ndfs, ndf_order, contrast_targets(cc), [y_min, y_max]);
        drawnow;

        combinedSingleFig = figure('Name', sprintf('Combined_Shifted_Gamma_Single_Plot_AGC_Target_%0.3g', contrast_targets(cc)));
        set(combinedSingleFig, 'Color', 'w');
        set(combinedSingleFig, 'Units', 'pixels', 'Position', scaled_figure_position(4, 2));
        figure_idx = figure_idx + 1;
        figureHandles{figure_idx} = combinedSingleFig;
        combinedSingleTiles = tiledlayout(combinedSingleFig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(combinedSingleTiles, sprintf('Shifted ND 0 and ND 0.4 Sensor Values | AGC Target %.3g', contrast_targets(cc)), ...
            'FontWeight', 'Bold', 'FontSize', 18);
        combinedSingleAxes = nexttile(combinedSingleTiles);
        plot_combined_shifted_gamma_single(combinedSingleAxes, settings_levels, rgb_mean_by_NDF, ...
            sorted_ndfs, ndf_order, contrast_targets(cc));
        drawnow;
    end

    % Save the exact curve data used by the plots, grouped first by contrast
    % target and then by sorted NDF.
    curves_mat_filename = 'camera_linearity_curves_by_NDF.mat';
    save(curves_mat_filename, 'curves_by_NDF');
    nd04_rgb_mat_filename = 'camera_linearity_ND0_ND0p4_rgb_means.mat';
    save(nd04_rgb_mat_filename, 'nd04_rgb_means');

    if(options.verbose)
        fprintf('Saved camera linearity curves by NDF to %s.\n', curves_mat_filename);
        fprintf('Saved ND 0 / ND 0.4 RGB mean measurements to %s.\n', nd04_rgb_mat_filename);
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

function nd04_rgb_means = initialize_nd_pair_rgb_export()
% Create a minimal export structure for the ND 0 / ND 0.4 means.

    nd04_rgb_means = struct();
    nd04_rgb_means.ND0 = empty_nd_mean_export();
    nd04_rgb_means.ND0p4 = empty_nd_mean_export();
end

function nd04_rgb_means = add_nd_pair_rgb_export_record(nd04_rgb_means, contrast_idx, ...
    settings_levels, sorted_ndfs, ndf_order, rgb_mean_by_NDF, frame_mean_by_NDF)
% Fill the minimal ND 0 / ND 0.4 export from the first contrast target.

    if(contrast_idx ~= 1)
        return;
    end

    target_ndfs = target_gamma_ndfs();
    nd0_rgb = rgb_mean_for_target_ndf(target_ndfs(1), sorted_ndfs, ndf_order, rgb_mean_by_NDF);
    nd04_rgb = rgb_mean_for_target_ndf(target_ndfs(2), sorted_ndfs, ndf_order, rgb_mean_by_NDF);
    nd0_frame = frame_mean_for_target_ndf(target_ndfs(1), sorted_ndfs, ndf_order, frame_mean_by_NDF);
    nd04_frame = frame_mean_for_target_ndf(target_ndfs(2), sorted_ndfs, ndf_order, frame_mean_by_NDF);

    nd04_rgb_means.ND0.settings_levels = settings_levels;
    nd04_rgb_means.ND0.rgb_mean = nd0_rgb;
    nd04_rgb_means.ND0.frame_mean = nd0_frame;

    nd04_rgb_means.ND0p4.settings_levels = settings_levels;
    nd04_rgb_means.ND0p4.rgb_mean = nd04_rgb;
    nd04_rgb_means.ND0p4.frame_mean = nd04_frame;
end

function nd_export = empty_nd_mean_export()
% Empty record for one ND level.

    nd_export = struct( ...
        'settings_levels', [], ...
        'rgb_mean', [], ...
        'frame_mean', []);
end

function rgb_mean = rgb_mean_for_target_ndf(target_ndf, sorted_ndfs, ndf_order, rgb_mean_by_NDF)
% Return a settings-level x RGB matrix for one requested NDF.

    n_settings = size(rgb_mean_by_NDF, 2);
    rgb_mean = nan(n_settings, 3);
    sorted_idx = find(abs(sorted_ndfs - target_ndf) < 1e-8, 1, 'first');

    if(isempty(sorted_idx))
        return;
    end

    nn = ndf_order(sorted_idx);
    rgb_mean = reshape(rgb_mean_by_NDF(nn, :, :), [n_settings, 3]);
end

function frame_mean = frame_mean_for_target_ndf(target_ndf, sorted_ndfs, ndf_order, frame_mean_by_NDF)
% Return a settings-level row vector for one requested NDF.

    n_settings = size(frame_mean_by_NDF, 2);
    frame_mean = nan(1, n_settings);
    sorted_idx = find(abs(sorted_ndfs - target_ndf) < 1e-8, 1, 'first');

    if(isempty(sorted_idx))
        return;
    end

    nn = ndf_order(sorted_idx);
    frame_mean = frame_mean_by_NDF(nn, :);
end

function target_ndfs = target_gamma_ndfs()
% NDF pair used for the shifted gamma comparison and matching RGB export.

    target_ndfs = [0, 0.4];
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

function plot_raw_response(ax, settings_levels, mean_response, error_response, rgb_response, ndf_value, y_limits)
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

    title(ax, ndf_label(ndf_value), 'FontWeight', 'Bold');
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
            'DisplayName', ndf_label(sorted_ndfs(oo)));
    end

    label_min_max_lines(ax, frame_mean_by_NDF(:));
    title(ax, sprintf('Curves Superimposed | AGC Target %.3g', contrast_target), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings Level');
    ylabel(ax, 'Frame Mean');
    ylim(ax, y_limits);
    expand_axes_right_for_line_labels(ax, settings_levels);
    for oo = 1:numel(ndf_order)
        nn = ndf_order(oo);
        label_curve_on_line(ax, settings_levels, frame_mean_by_NDF(nn, :), ...
            ndf_label(sorted_ndfs(oo)), curve_colors(oo, :));
    end
    grid(ax, 'on');
    hold(ax, 'off');
end

function [black_level, white_level] = gamma_response_anchors()
% Fixed B-channel response anchors used for gamma fitting.

    black_level = 16;
    white_level = 254;
end

function [normalized_settings, normalized_response] = normalize_gamma_source_curve(settings_levels, response, black_level, white_level)
% Normalize one B-channel response curve using fixed black/white camera-count anchors.

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

function plot_b_channel_raw_vs_gamma_fit(ax, settings_levels, b_response, gamma_fit, ndf_value, show_legend)
% Plot one NDF B-channel raw curve against its own gamma fit.

    settings_levels = settings_levels(:);
    b_response = b_response(:);
    valid = isfinite(settings_levels) & isfinite(b_response);
    [b_color, fit_color] = b_gamma_plot_colors();

    hold(ax, 'on');

    if(any(valid))
        plot(ax, settings_levels(valid), b_response(valid), '-o', ...
            'Color', b_color, ...
            'MarkerFaceColor', b_color, ...
            'MarkerEdgeColor', b_color, ...
            'MarkerSize', 5, ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Raw B channel');
    end

    fitted_response = evaluate_b_channel_gamma(settings_levels, gamma_fit);
    fit_valid = isfinite(settings_levels) & isfinite(fitted_response);

    if(any(fit_valid))
        plot(ax, settings_levels(fit_valid), fitted_response(fit_valid), '-', ...
            'Color', fit_color, ...
            'LineWidth', 2.0, ...
            'DisplayName', sprintf('%s gamma %.3g | %.0f-%.0f counts', ...
                gamma_fit.label, gamma_fit.gamma, gamma_fit.black_level, gamma_fit.white_level));
    end

    title(ax, ndf_label(ndf_value), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings Level');
    ylabel(ax, 'B Mean');
    format_b_gamma_axes(ax, b_response, fitted_response);
    expand_axes_right_for_line_labels(ax, settings_levels);
    if(any(valid))
        label_curve_on_line(ax, settings_levels(valid), b_response(valid), 'B', b_color, -0.035);
    end
    if(any(fit_valid))
        label_curve_on_line(ax, settings_levels(fit_valid), fitted_response(fit_valid), 'Gamma', fit_color, 0.035);
    end
    grid(ax, 'on');
    if(show_legend)
        legend(ax, 'Location', 'eastoutside', 'FontSize', 8, 'Box', 'off');
    end
    hold(ax, 'off');
end

function [b_color, fit_color] = b_gamma_plot_colors()
% Colors for B-channel data and fitted gamma curves.

    b_color = [0.08, 0.28, 0.85];
    fit_color = [0.92, 0.58, 0.16];
end

function fitted_response = evaluate_b_channel_gamma(settings_levels, gamma_fit)
% Apply a fixed-anchor B-channel gamma curve in camera-count space.

    fitted_response = nan(size(settings_levels));

    if(~isfinite(gamma_fit.gamma) || ~isfinite(gamma_fit.black_level) || ...
            ~isfinite(gamma_fit.white_level) || gamma_fit.white_level <= gamma_fit.black_level)
        return;
    end

    normalized_settings = min(max(settings_levels, 0), inf);
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

function plot_gamma_fitting_shift(tiles, settings_levels, frame_mean_by_NDF, frame_std_by_NDF, rgb_mean_by_NDF, ...
    agc_settings_by_NDF, ...
    sorted_ndfs, ndf_order, contrast_target, frame_y_limits)
% Plot ND 0 / ND 0.4 B-channel gamma fits before and after applying the settings shift.

    target_ndfs = target_gamma_ndfs();
    curve_colors = [
        0.05, 0.25, 0.60
        0.80, 0.25, 0.10
        ];
    gamma_fits = repmat(empty_gamma_fit(), numel(target_ndfs), 1);
    target_matrix_indices = nan(size(target_ndfs));

    for tt = 1:numel(target_ndfs)
        gamma_fits(tt).label = nd_label(target_ndfs(tt));
        sorted_idx = find(abs(sorted_ndfs - target_ndfs(tt)) < 1e-8, 1, 'first');

        if(isempty(sorted_idx))
            continue;
        end

        nn = ndf_order(sorted_idx);
        target_matrix_indices(tt) = nn;
        b_mean = reshape(rgb_mean_by_NDF(nn, :, 3), [1, numel(settings_levels)]);
        gamma_fits(tt) = fit_b_channel_gamma(settings_levels, b_mean, nd_label(target_ndfs(tt)));
    end

    combined_fit = fit_shifted_combined_b_gamma(gamma_fits(1), gamma_fits(2));
    b_y_limits = combined_gamma_response_limits(gamma_fits, combined_fit);
    rawAxes = gobjects(1, numel(target_ndfs));
    settingsAxes = gobjects(1, numel(target_ndfs));

    for tt = 1:numel(target_ndfs)
        settingsAxes(tt) = nexttile(tiles, tt);
        nn = target_matrix_indices(tt);

        if(isfinite(nn))
            plot_agc_settings(settingsAxes(tt), agc_settings_by_NDF(nn, :), target_ndfs(tt));
        else
            show_missing_target_message(settingsAxes(tt), nd_label(target_ndfs(tt)));
        end
    end

    for tt = 1:numel(target_ndfs)
        rawAxes(tt) = nexttile(tiles, 2 + tt);
        nn = target_matrix_indices(tt);

        if(isfinite(nn))
            rgb_response = reshape(rgb_mean_by_NDF(nn, :, :), [numel(settings_levels), 3]);
            plot_raw_response(rawAxes(tt), settings_levels, frame_mean_by_NDF(nn, :), ...
                frame_std_by_NDF(nn, :), rgb_response, target_ndfs(tt), frame_y_limits);
        else
            show_missing_target_message(rawAxes(tt), nd_label(target_ndfs(tt)));
        end
    end

    beforeAxes = nexttile(tiles, 5);
    plot_combined_gamma_before_shift(beforeAxes, gamma_fits, curve_colors, contrast_target, combined_fit, b_y_limits);

    shiftedAxes = nexttile(tiles, 6);
    plot_combined_gamma_after_shift(shiftedAxes, gamma_fits, curve_colors, contrast_target, combined_fit, b_y_limits);

    drawnow;
    add_centered_row_title(settingsAxes, 'Camera settings');
    add_centered_row_title(rawAxes, 'World camera frame mean raw');
    add_centered_row_title([beforeAxes, shiftedAxes], 'Single-function fit with settings shift');
end

function plot_combined_shifted_gamma_single(ax, settings_levels, rgb_mean_by_NDF, sorted_ndfs, ndf_order, contrast_target)
% Plot shifted ND 0 and ND 0.4 sensor values with one shared gamma function.

    [gamma_fits, combined_fit, curve_colors] = build_target_combined_gamma_fit(settings_levels, rgb_mean_by_NDF, ...
        sorted_ndfs, ndf_order);
    y_limits = combined_gamma_response_limits(gamma_fits, combined_fit);

    plot_combined_gamma_after_shift(ax, gamma_fits, curve_colors, contrast_target, combined_fit, y_limits);
    title(ax, sprintf('Settings vs B Sensor Values | S = %.4g | AGC Target %.3g', ...
        combined_fit.settings_shift, contrast_target), 'FontWeight', 'Bold');
    xlabel(ax, 'Settings');
    ylabel(ax, 'B Sensor Value');
end

function [gamma_fits, combined_fit, curve_colors] = build_target_combined_gamma_fit(settings_levels, rgb_mean_by_NDF, sorted_ndfs, ndf_order)
% Build ND 0 / ND 0.4 B-channel gamma fits and their shared shifted model.

    target_ndfs = target_gamma_ndfs();
    curve_colors = [
        0.05, 0.25, 0.60
        0.80, 0.25, 0.10
        ];
    gamma_fits = repmat(empty_gamma_fit(), numel(target_ndfs), 1);

    for tt = 1:numel(target_ndfs)
        gamma_fits(tt).label = nd_label(target_ndfs(tt));
        sorted_idx = find(abs(sorted_ndfs - target_ndfs(tt)) < 1e-8, 1, 'first');

        if(isempty(sorted_idx))
            continue;
        end

        nn = ndf_order(sorted_idx);
        b_mean = reshape(rgb_mean_by_NDF(nn, :, 3), [1, numel(settings_levels)]);
        gamma_fits(tt) = fit_b_channel_gamma(settings_levels, b_mean, nd_label(target_ndfs(tt)));
    end

    combined_fit = fit_shifted_combined_b_gamma(gamma_fits(1), gamma_fits(2));
end

function add_centered_row_title(row_axes, row_title)
% Place one centered title above a row of tiled axes.

    row_axes = row_axes(isgraphics(row_axes, 'axes'));

    if(isempty(row_axes))
        return;
    end

    fig = ancestor(row_axes(1), 'figure');
    positions = nan(numel(row_axes), 4);

    for aa = 1:numel(row_axes)
        old_units = row_axes(aa).Units;
        row_axes(aa).Units = 'normalized';
        positions(aa, :) = row_axes(aa).Position;
        row_axes(aa).Units = old_units;
    end

    left = min(positions(:, 1));
    right = max(positions(:, 1) + positions(:, 3));
    top = max(positions(:, 2) + positions(:, 4));
    width = right - left;
    title_height = 0.025;
    bottom = min(0.975 - title_height, top + 0.006);

    annotation(fig, 'textbox', [left, bottom, width, title_height], ...
        'String', row_title, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'LineStyle', 'none', ...
        'FitBoxToText', 'off');
end

function show_missing_target_message(ax, target_label)
% Fill an expected simplification-panel slot when the target NDF is absent.

    title(ax, target_label, 'FontWeight', 'Bold');
    text(ax, 0.5, 0.5, sprintf('%s not found', target_label), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 12, ...
        'FontWeight', 'bold');
    grid(ax, 'on');
end

function gamma_fit = empty_gamma_fit()
% Construct an empty gamma-fit record.

    gamma_fit = struct( ...
        'label', '', ...
        'settings', [], ...
        'response', [], ...
        'x_fit', [], ...
        'y_fit', [], ...
        'gamma', nan, ...
        'black_level', nan, ...
        'white_level', nan, ...
        'y_min', nan, ...
        'y_max', nan, ...
        'valid', false);
end

function gamma_fit = fit_b_channel_gamma(settings_levels, b_mean, label_text)
% Fit B = black + (white - black) * settings^gamma using fixed anchors.

    gamma_fit = empty_gamma_fit();
    gamma_fit.label = label_text;
    [black_level, white_level] = gamma_response_anchors();

    settings_levels = double(settings_levels(:));
    b_mean = double(b_mean(:));
    valid = isfinite(settings_levels) & isfinite(b_mean) & settings_levels >= 0;

    if(nnz(valid) < 2)
        return;
    end

    settings_levels = settings_levels(valid);
    b_mean = b_mean(valid);
    [settings_levels, sort_order] = sort(settings_levels);
    b_mean = b_mean(sort_order);

    [normalized_settings, normalized_response] = normalize_gamma_source_curve(settings_levels, b_mean, ...
        black_level, white_level);

    if(numel(normalized_settings) < 2)
        return;
    end

    gamma_value = fit_gamma_exponent(normalized_settings, normalized_response);
    x_fit = linspace(min(settings_levels), max(settings_levels), 300);
    fit_struct = struct('gamma', gamma_value, 'black_level', black_level, 'white_level', white_level);
    y_fit = evaluate_b_channel_gamma(x_fit, fit_struct);

    gamma_fit.settings = settings_levels;
    gamma_fit.response = b_mean;
    gamma_fit.x_fit = x_fit;
    gamma_fit.y_fit = y_fit;
    gamma_fit.gamma = gamma_value;
    gamma_fit.black_level = black_level;
    gamma_fit.white_level = white_level;
    gamma_fit.y_min = black_level;
    gamma_fit.y_max = white_level;
    gamma_fit.valid = true;
end

function gamma_value = fit_gamma_exponent(settings_levels, normalized_response)
% Fit response = settings_level^gamma using a positive gamma parameter.

    valid = isfinite(settings_levels) & isfinite(normalized_response) & ...
        settings_levels >= 0 & settings_levels <= 1;
    settings_levels = settings_levels(valid);
    normalized_response = normalized_response(valid);

    if(numel(settings_levels) < 2)
        gamma_value = 1;
        return;
    end

    fit_objective = @(log_gamma) sum((settings_levels .^ exp(log_gamma) - normalized_response) .^ 2, 'omitnan');
    gamma_value = exp(fminsearch(fit_objective, 0, optimset('Display', 'off')));
end

function shifted_settings = inverse_b_channel_gamma(response, gamma_fit)
% Convert B-channel values back to settings values using one gamma fit.

    shifted_settings = nan(size(response));

    if(~gamma_fit.valid || gamma_fit.y_max <= gamma_fit.y_min || ~isfinite(gamma_fit.gamma) || gamma_fit.gamma <= 0)
        return;
    end

    normalized_response = (response - gamma_fit.y_min) ./ (gamma_fit.y_max - gamma_fit.y_min);
    valid = isfinite(normalized_response) & normalized_response >= 0 & normalized_response <= 1;
    shifted_settings(valid) = normalized_response(valid) .^ (1 / gamma_fit.gamma);
end

function shift_value = find_gamma_settings_shift(source_fit, target_fit)
% Find S so source settings + S best match target settings at equal response.

    shift_value = nan;

    if(~source_fit.valid || ~target_fit.valid)
        return;
    end

    y_overlap_min = max(source_fit.y_min, target_fit.y_min);
    y_overlap_max = min(source_fit.y_max, target_fit.y_max);

    if(~isfinite(y_overlap_min) || ~isfinite(y_overlap_max) || y_overlap_max <= y_overlap_min)
        return;
    end

    response_samples = linspace(y_overlap_min, y_overlap_max, 300);
    source_settings = inverse_b_channel_gamma(response_samples, source_fit);
    target_settings = inverse_b_channel_gamma(response_samples, target_fit);
    valid = isfinite(source_settings) & isfinite(target_settings);

    if(~any(valid))
        return;
    end

    shift_value = mean(target_settings(valid) - source_settings(valid), 'omitnan');
end

function combined_fit = fit_shifted_combined_b_gamma(source_fit, target_fit)
% Fit one B-channel gamma function to raw values after shifting source settings by S.

    combined_fit = struct( ...
        'valid', false, ...
        'gamma', nan, ...
        'settings_shift', nan, ...
        'black_level', nan, ...
        'white_level', nan, ...
        'x_fit', [], ...
        'y_fit', []);

    if(~source_fit.valid || ~target_fit.valid)
        return;
    end

    source_settings = source_fit.settings(:);
    target_settings = target_fit.settings(:);
    source_response = source_fit.response(:);
    target_response = target_fit.response(:);
    valid_source = isfinite(source_settings) & isfinite(source_response) & source_settings >= 0;
    valid_target = isfinite(target_settings) & isfinite(target_response) & target_settings >= 0;

    if(nnz(valid_source) < 2 || nnz(valid_target) < 2)
        return;
    end

    source_settings = source_settings(valid_source);
    source_response = source_response(valid_source);
    target_settings = target_settings(valid_target);
    target_response = target_response(valid_target);

    initial_shift = find_gamma_settings_shift(source_fit, target_fit);
    if(~isfinite(initial_shift) || initial_shift <= 0)
        initial_shift = 0.1;
    end
    initial_gamma = mean([source_fit.gamma, target_fit.gamma], 'omitnan');
    if(~isfinite(initial_gamma) || initial_gamma <= 0)
        initial_gamma = 1;
    end
    response_values = [source_response; target_response];
    initial_black_level = min(response_values, [], 'omitnan');
    initial_response_scale = max(response_values, [], 'omitnan') - initial_black_level;
    if(~isfinite(initial_response_scale) || initial_response_scale <= 0)
        initial_response_scale = 1;
    end

    fit_objective = @(params) combined_gamma_objective(params, ...
        source_settings, source_response, target_settings, target_response);
    fitted_params = fminsearch(fit_objective, ...
        [log(initial_gamma), log(initial_shift), initial_black_level, log(initial_response_scale)], ...
        optimset('Display', 'off'));
    gamma_value = exp(fitted_params(1));
    settings_shift = exp(fitted_params(2));
    black_level = fitted_params(3);
    response_scale = exp(fitted_params(4));
    white_level = black_level + response_scale;

    combined_x = [target_settings; source_settings + settings_shift];
    x_fit = linspace(min(combined_x), max(combined_x), 400);
    y_fit = evaluate_b_channel_gamma(x_fit, struct( ...
        'gamma', gamma_value, ...
        'black_level', black_level, ...
        'white_level', white_level));

    combined_fit.valid = true;
    combined_fit.gamma = gamma_value;
    combined_fit.settings_shift = settings_shift;
    combined_fit.black_level = black_level;
    combined_fit.white_level = white_level;
    combined_fit.x_fit = x_fit;
    combined_fit.y_fit = y_fit;
end

function error_value = combined_gamma_objective(params, source_settings, source_response, target_settings, target_response)
% Least-squares objective for one shifted source plus unshifted target gamma in raw sensor units.

    gamma_value = exp(params(1));
    settings_shift = exp(params(2));
    black_level = params(3);
    response_scale = exp(params(4));
    predicted_source = black_level + response_scale .* ((source_settings + settings_shift) .^ gamma_value);
    predicted_target = black_level + response_scale .* (target_settings .^ gamma_value);
    residuals = [predicted_source - source_response; predicted_target - target_response];
    error_value = sum(residuals .^ 2, 'omitnan');
end

function y_limits = combined_gamma_response_limits(gamma_fits, combined_fit)
% Shared y-range for separate points and the combined gamma model.

    response_values = [];

    for gg = 1:numel(gamma_fits)
        if(~gamma_fits(gg).valid)
            continue;
        end

        response_values = [response_values; gamma_fits(gg).response(:)]; %#ok<AGROW>
    end

    if(combined_fit.valid)
        response_values = [response_values; combined_fit.y_fit(:)];
    end

    response_values = response_values(isfinite(response_values));

    if(isempty(response_values))
        y_limits = [0, 300];
        return;
    end

    y_min = min(response_values);
    y_max = max(response_values);

    if(y_max <= y_min)
        y_padding = max(abs(y_max), 1) * 0.05;
    else
        y_padding = 0.08 * (y_max - y_min);
    end

    y_limits = [min(0, y_min - y_padding), y_max + y_padding];
end

function plot_combined_gamma_before_shift(ax, gamma_fits, curve_colors, contrast_target, combined_fit, y_limits)
% Show both ND point sets before applying the fitted source settings shift.

    hold(ax, 'on');
    plot_gamma_points(ax, gamma_fits, curve_colors, false, combined_fit);
    plot_combined_gamma_line(ax, combined_fit);
    title(ax, sprintf('Before shift | AGC Target %.3g', contrast_target), 'FontWeight', 'Bold');
    xlabel(ax, 'CombiLED Settings');
    ylabel(ax, 'B Mean');
    ylim(ax, y_limits);
    expand_combined_gamma_xlim(ax, gamma_fits, combined_fit, false);
    label_combined_gamma_plot(ax, gamma_fits, curve_colors, combined_fit, false);
    grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside', 'FontSize', 9, 'Box', 'off');
    show_missing_gamma_fit_message(ax, gamma_fits);
    hold(ax, 'off');
end

function plot_combined_gamma_after_shift(ax, gamma_fits, curve_colors, contrast_target, combined_fit, y_limits)
% Show shifted ND 0 plus ND 0.4 data explained by one gamma function.

    hold(ax, 'on');
    plot_gamma_points(ax, gamma_fits, curve_colors, true, combined_fit);
    plot_combined_gamma_line(ax, combined_fit);
    if(combined_fit.valid)
        title_text = sprintf('After shift | S = %.4g | AGC Target %.3g', combined_fit.settings_shift, contrast_target);
    else
        title_text = sprintf('After shift | AGC Target %.3g', contrast_target);
    end
    title(ax, title_text, 'FontWeight', 'Bold');
    xlabel(ax, 'CombiLED Settings');
    ylabel(ax, 'B Mean');
    ylim(ax, y_limits);
    expand_combined_gamma_xlim(ax, gamma_fits, combined_fit, true);
    label_combined_gamma_plot(ax, gamma_fits, curve_colors, combined_fit, true);
    grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside', 'FontSize', 9, 'Box', 'off');
    show_missing_gamma_fit_message(ax, gamma_fits);
    if(~combined_fit.valid)
        text(ax, 0.5, 0.5, 'Could not fit combined shifted gamma', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
    hold(ax, 'off');
end

function plot_gamma_points(ax, gamma_fits, curve_colors, apply_shift, combined_fit)
% Plot ND 0 / ND 0.4 B-channel points, optionally shifting ND 0 settings.

    for tt = 1:numel(gamma_fits)
        if(~gamma_fits(tt).valid)
            continue;
        end

        x_values = gamma_fits(tt).settings;
        display_label = gamma_fits(tt).label;
        if(apply_shift && tt == 1 && combined_fit.valid)
            x_values = x_values + combined_fit.settings_shift;
            display_label = sprintf('%s shifted', gamma_fits(tt).label);
        end

        plot(ax, x_values, gamma_fits(tt).response, 'o', ...
            'Color', curve_colors(tt, :), ...
            'MarkerFaceColor', curve_colors(tt, :), ...
            'MarkerEdgeColor', curve_colors(tt, :), ...
            'MarkerSize', 5, ...
            'DisplayName', display_label);
    end
end

function plot_combined_gamma_line(ax, combined_fit)
% Plot the single combined gamma function.

    if(~combined_fit.valid)
        return;
    end

    plot(ax, combined_fit.x_fit, combined_fit.y_fit, '-', ...
        'Color', [0.92, 0.58, 0.16], ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf('Gamma %.3g', combined_fit.gamma));
end

function expand_combined_gamma_xlim(ax, gamma_fits, combined_fit, apply_shift)
% Set x-limits to include shifted settings values and direct labels.

    x_values = [];

    for tt = 1:numel(gamma_fits)
        if(~gamma_fits(tt).valid)
            continue;
        end

        x_this = gamma_fits(tt).settings(:);
        if(apply_shift && tt == 1 && combined_fit.valid)
            x_this = x_this + combined_fit.settings_shift;
        end
        x_values = [x_values; x_this]; %#ok<AGROW>
    end

    if(combined_fit.valid)
        x_values = [x_values; combined_fit.x_fit(:)];
    end

    expand_axes_right_for_line_labels(ax, x_values);
end

function label_combined_gamma_plot(ax, gamma_fits, curve_colors, combined_fit, apply_shift)
% Place separated labels on the combined gamma plots.

    if(numel(gamma_fits) >= 1 && gamma_fits(1).valid)
        x_values = gamma_fits(1).settings;
        label_text = gamma_fits(1).label;
        if(apply_shift && combined_fit.valid)
            x_values = x_values + combined_fit.settings_shift;
            label_text = sprintf('%s +S', gamma_fits(1).label);
        end
        label_curve_at_extreme(ax, x_values, gamma_fits(1).response, label_text, curve_colors(1, :), 'min');
    end

    if(numel(gamma_fits) >= 2 && gamma_fits(2).valid)
        label_curve_at_extreme(ax, gamma_fits(2).settings, gamma_fits(2).response, ...
            gamma_fits(2).label, curve_colors(2, :), 'max');
    end

    if(combined_fit.valid)
        label_curve_on_line(ax, combined_fit.x_fit, combined_fit.y_fit, 'Gamma', [0.92, 0.58, 0.16], 0.04, 0.55);
    end
end

function show_missing_gamma_fit_message(ax, gamma_fits)
% Report missing ND curves directly inside the gamma axes.

    missing_labels = {};

    for tt = 1:numel(gamma_fits)
        if(~gamma_fits(tt).valid)
            if(strlength(gamma_fits(tt).label) == 0)
                missing_labels{end + 1} = sprintf('target curve %d', tt); %#ok<AGROW>
            else
                missing_labels{end + 1} = gamma_fits(tt).label; %#ok<AGROW>
            end
        end
    end

    if(isempty(missing_labels))
        return;
    end

    text(ax, 0.5, 0.08, sprintf('Missing fit: %s', strjoin(missing_labels, ', ')), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2);
end

function plot_agc_settings(ax, agc_settings, ndf_value)
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

    title(ax, sprintf('%s AGC', ndf_label(ndf_value)), 'FontWeight', 'Bold');
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

function label = ndf_label(ndf_value)
% Return the compact label used for NDF curve annotations.

    label = sprintf('NDF %.3g', ndf_value);
end

function label = nd_label(ndf_value)
% Return the compact label used for ND settings-shift annotations.

    label = sprintf('ND %.3g', ndf_value);
end

function expand_axes_right_for_line_labels(ax, x_values)
% Add right-side room so direct line labels are readable.

    x_values = double(x_values(:));
    x_values = x_values(isfinite(x_values));

    if(isempty(x_values))
        return;
    end

    x_min = min(x_values);
    x_max = max(x_values);
    x_span = x_max - x_min;

    if(x_span <= 0)
        x_span = max(abs(x_max), 1);
    end

    xlim(ax, [x_min - 0.03 * x_span, x_max + 0.18 * x_span]);
end

function label_curve_on_line(ax, settings_levels, response_values, label_text, label_color, y_offset_fraction, x_fraction)
% Put a small NDF/ND label directly at the end of a plotted curve.

    if(nargin < 6)
        y_offset_fraction = 0;
    end
    if(nargin < 7)
        x_fraction = 1;
    end

    settings_levels = double(settings_levels(:));
    response_values = double(response_values(:));
    valid = isfinite(settings_levels) & isfinite(response_values);

    if(~any(valid))
        return;
    end

    valid_indices = find(valid);
    x_fraction = min(max(x_fraction, 0), 1);
    label_idx = valid_indices(max(1, round(1 + x_fraction * (numel(valid_indices) - 1))));
    x_limits = xlim(ax);
    y_limits = ylim(ax);
    x_offset = 0.012 * diff(x_limits);
    y_offset = y_offset_fraction * diff(y_limits);

    text(ax, settings_levels(label_idx) + x_offset, response_values(label_idx) + y_offset, label_text, ...
        'Color', label_color, ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'w', ...
        'Margin', 1, ...
        'Clipping', 'off');
end

function label_curve_at_extreme(ax, settings_levels, response_values, label_text, label_color, extreme_name)
% Put a line label near the minimum or maximum y-value on a curve.

    settings_levels = double(settings_levels(:));
    response_values = double(response_values(:));
    valid = isfinite(settings_levels) & isfinite(response_values);

    if(~any(valid))
        return;
    end

    valid_indices = find(valid);

    switch string(extreme_name)
        case "min"
            [~, local_idx] = min(response_values(valid_indices));
            y_offset_fraction = -0.035;
            vertical_alignment = 'top';
        case "max"
            [~, local_idx] = max(response_values(valid_indices));
            y_offset_fraction = 0.035;
            vertical_alignment = 'bottom';
        otherwise
            local_idx = numel(valid_indices);
            y_offset_fraction = 0;
            vertical_alignment = 'middle';
    end

    label_idx = valid_indices(local_idx);
    x_limits = xlim(ax);
    y_limits = ylim(ax);
    x_offset = 0.012 * diff(x_limits);
    y_offset = y_offset_fraction * diff(y_limits);

    text(ax, settings_levels(label_idx) + x_offset, response_values(label_idx) + y_offset, label_text, ...
        'Color', label_color, ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', vertical_alignment, ...
        'BackgroundColor', 'w', ...
        'Margin', 1, ...
        'Clipping', 'off');
end

function label_max_lines_outside_right(ax, setting_labels, setting_maxima, setting_colors)
% Put compact max labels outside both sides of the AGC axes.

    x_limits = xlim(ax);
    x_label_left = x_limits(1) - 0.05 * diff(x_limits);
    x_label_right = x_limits(2) + 0.05 * diff(x_limits);

    for ss = 1:numel(setting_maxima)
        if(~isfinite(setting_maxima(ss)))
            continue;
        end

        if(ss == 1)
            x_label = x_label_left;
            horizontal_alignment = 'right';
        else
            x_label = x_label_right;
            horizontal_alignment = 'left';
        end

        text(ax, x_label, setting_maxima(ss), setting_labels{ss}, ...
            'Color', setting_colors(ss, :), ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', horizontal_alignment, ...
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
    x_label_left = x_limits(1) + 0.02 * diff(x_limits);
    x_label_right = x_limits(2) - 0.02 * diff(x_limits);

    yline(ax, min_value, '--', ...
        'Color', [0.85, 0.10, 0.10], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    yline(ax, max_value, '--', ...
        'Color', [0.15, 0.60, 0.25], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    text(ax, x_label_right, min_value, sprintf('min %.1f', min_value), ...
        'Color', [0.65, 0.15, 0.15], ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
    text(ax, x_label_left, max_value, sprintf('max %.1f', max_value), ...
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

    x_limits = xlim(ax);
    x_span = diff(x_limits);
    x_position = settings_levels(max_idx) + 0.025 * x_span;

    plot(ax, settings_levels(max_idx), max_value, 'p', ...
        'Color', [0.05, 0.35, 0.15], ...
        'MarkerFaceColor', [0.05, 0.35, 0.15], ...
        'MarkerSize', 8, ...
        'HandleVisibility', 'off');
    text(ax, x_position, max_value, ...
        sprintf('max %.1f\nR %.1f  G %.1f  B %.1f', max_value, rgb_values(1), rgb_values(2), rgb_values(3)), ...
        'Color', [0.05, 0.25, 0.10], ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2, ...
        'Clipping', 'off');
end
