function figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
    arguments
        calibration_metadata
        measurements
        options.verbose logical = true;
    end

    figureHandles = {};

    % Defaults for all of the plotting we will do
    set(groot, 'DefaultAxesFontSize', 16);
    set(groot, 'DefaultTextFontSize', 16);
    set(groot, 'DefaultAxesColor', 'w');
    set(groot, 'DefaultFigureColor', 'w');

    [n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements] = size(measurements);
    title_font_size = 16;
    label_font_size = 14;
    tick_font_size = 12;
    legend_font_size = 12;

    settings_levels = calibration_metadata.background_scalars;
    NDFs = calibration_metadata.NDFs;

    cameraLinearityFig = figure('Name', 'World_Camera_Linearity_All_NDFs');
    set(cameraLinearityFig, 'Color', 'w');
    figureHandles{end+1, 1} = cameraLinearityFig;

    responseTiled = tiledlayout(cameraLinearityFig, 1, n_contrastTargets, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(responseTiled, 'World Camera Linearity | All NDFs', 'FontWeight', 'Bold');

    cameraSettingsFig = figure('Name', 'World_Camera_Settings_By_NDF');
    set(cameraSettingsFig, 'Color', 'w');
    figureHandles{end+1, 1} = cameraSettingsFig;

    settingsTiled = tiledlayout(cameraSettingsFig, 1, n_contrastTargets, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(settingsTiled, 'World Camera Settings by NDF', 'FontWeight', 'Bold');

    agc_settings_by_NDF = nan(n_NDFs, n_contrastTargets, 3);

    for cc = 1:n_contrastTargets
        contrast_target = calibration_metadata.contrast_agc_targets(cc);
        responseAx = nexttile(responseTiled);
        hold(responseAx, 'on');

        for nn = 1:n_NDFs
            NDF = NDFs(nn);

            if(options.verbose)
                fprintf("World Camera Linearity | NDF (%d/%d): %.3f | Contrast Target (%d/%d)\n", ...
                    nn, n_NDFs, NDF, cc, n_contrastTargets);
            end

            settings_level_measurements = reshape(measurements(nn, cc, :, :), [n_settingsLevels, n_measurements]);
            settings_levels_averaged = average_settings_levels(settings_level_measurements);
            agc_settings_by_NDF(nn, cc, :) = extract_first_world_settings(settings_level_measurements);

            channel_colors = ndf_channel_colors(nn, n_NDFs);
            plot(responseAx, settings_levels, settings_levels_averaged(:, 1), '-o', 'Color', channel_colors.R, ...
                'LineWidth', 1.5, 'DisplayName', sprintf('R | NDF %.3g', NDF));
            plot(responseAx, settings_levels, settings_levels_averaged(:, 2), '-o', 'Color', channel_colors.G, ...
                'LineWidth', 1.5, 'DisplayName', sprintf('G | NDF %.3g', NDF));
            plot(responseAx, settings_levels, settings_levels_averaged(:, 3), '-o', 'Color', channel_colors.B, ...
                'LineWidth', 1.5, 'DisplayName', sprintf('B | NDF %.3g', NDF));
            plot(responseAx, settings_levels, settings_levels_averaged(:, 4), '-o', 'Color', channel_colors.M, ...
                'LineWidth', 1.5, 'DisplayName', sprintf('Frame Mean | NDF %.3g', NDF));
        end

        title(responseAx, sprintf("Contrast Target: %.2f", contrast_target), 'FontWeight', 'Bold', 'FontSize', title_font_size);
        xlabel(responseAx, 'Settings Level', 'FontSize', label_font_size);
        ylabel(responseAx, 'Mean Response', 'FontSize', label_font_size);
        responseAx.FontSize = tick_font_size;
        ylim(responseAx, [0, 255]);
        grid(responseAx, 'on');
        legend(responseAx, 'Location', 'eastoutside', 'FontSize', legend_font_size);
        hold(responseAx, 'off');

        settingsAx = nexttile(settingsTiled);
        contrast_agc_settings = reshape(agc_settings_by_NDF(:, cc, :), [n_NDFs, 3]);
        plot_camera_settings_by_NDF(settingsAx, NDFs, contrast_agc_settings, contrast_target, ...
            title_font_size, label_font_size, tick_font_size, legend_font_size);
        drawnow;
    end
end


function settings_levels_averaged = average_settings_levels(settings_levels_measurements)
    [n_levels, n_measurements] = size(settings_levels_measurements);
    settings_level_measurements = nan(n_levels, n_measurements, 4);

    for ll = 1:n_levels
        for mm = 1:n_measurements
            measurement = settings_levels_measurements{ll, mm};
            world_values = measurement.W.v;

            validateattributes(world_values, {'numeric'}, {'2d', 'ncols', 4}, ...
                mfilename, 'measurement.W.v');

            % For differentiate_color=true with mean_axes.W=[1, 2], W.v is
            % [n_frames, 4] with columns [R_mean, G_mean, B_mean, frame_mean].
            settings_level_measurements(ll, mm, :) = mean(world_values, 1, 'omitnan');
        end
    end

    settings_levels_averaged = squeeze(mean(settings_level_measurements, 2, 'omitnan'));
end

function world_settings = extract_first_world_settings(settings_levels_measurements)
    world_settings = [nan, nan, nan];

    for ii = 1:numel(settings_levels_measurements)
        measurement = settings_levels_measurements{ii};

        if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'settings'))
            continue;
        end

        settings = measurement.W.settings;
        world_settings = [first_setting_value(settings, 'Again'), ...
                          first_setting_value(settings, 'Dgain'), ...
                          first_setting_value(settings, 'exposure')];
        return;
    end
end

function value = first_setting_value(settings, setting_name)
    value = nan;

    if(~isfield(settings, setting_name))
        return;
    end

    setting_values = settings.(setting_name);
    if(isempty(setting_values))
        return;
    end

    value = setting_values(1);
end

function channel_colors = ndf_channel_colors(ndf_idx, n_NDFs)
    color_scale = ndf_color_scale(ndf_idx, n_NDFs);

    channel_colors.R = tint_color([1.00, 0.00, 0.00], color_scale);
    channel_colors.G = tint_color([0.00, 0.55, 0.00], color_scale);
    channel_colors.B = tint_color([0.00, 0.20, 1.00], color_scale);
    channel_colors.M = tint_color([0.00, 0.00, 0.00], color_scale);
end

function color_scale = ndf_color_scale(ndf_idx, n_NDFs)
    if(n_NDFs <= 1)
        color_scale = 1;
        return;
    end

    color_scale = 0.35 + 0.65 * (ndf_idx - 1) / (n_NDFs - 1);
end

function color = tint_color(base_color, color_scale)
    color = 1 - color_scale * (1 - base_color);
end

function plot_camera_settings_by_NDF(ax, NDFs, settings_values, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
    hold(ax, 'on');

    yyaxis(ax, 'left');
    plot(ax, NDFs, settings_values(:, 1), '-o', 'Color', [0.00, 0.35, 0.75], ...
        'LineWidth', 1.5, 'DisplayName', 'Analogue Gain');
    plot(ax, NDFs, settings_values(:, 2), '-o', 'Color', [0.65, 0.20, 0.75], ...
        'LineWidth', 1.5, 'DisplayName', 'Digital Gain');
    ylabel(ax, 'Gain', 'FontSize', label_font_size);

    yyaxis(ax, 'right');
    plot(ax, NDFs, settings_values(:, 3), '-o', 'Color', [0.95, 0.55, 0.10], ...
        'LineWidth', 1.5, 'DisplayName', 'Exposure');
    ylabel(ax, 'Exposure', 'FontSize', label_font_size);

    title(ax, sprintf("Contrast Target: %.2f", contrast_target), 'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax, 'NDF', 'FontSize', label_font_size);
    ax.FontSize = tick_font_size;
    grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside', 'FontSize', legend_font_size);
    hold(ax, 'off');
end
