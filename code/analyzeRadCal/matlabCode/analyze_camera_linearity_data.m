function figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
% Analyze camera linearity data.
%
% Syntax:
%   figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
%
% Description:
%   This function analyze camera linearity data.
% Inputs:
%   calibration_metadata     - Input used by the function.
%   measurements             - Input used by the function.
%   options                  - Input used by the function.
%
% Outputs:
%   figureHandles            - Output produced by the function.
%
% Examples:
%{
    figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
%}

    arguments
        calibration_metadata
        measurements
        options.verbose logical = true;
    end

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
    channel_labels = {'R Mean', 'G Mean', 'B Mean', 'Frame Mean'};
    channel_names = {'R', 'G', 'B', 'Frame_Mean'};
    ndf_line_colors = high_contrast_ndf_colors(n_NDFs);
    responseTiles = cell(1, numel(channel_labels));
    figureHandles = {};

    for ch = 1:numel(channel_labels)
        responseFig = figure('Name', sprintf('World_Camera_Linearity_%s_All_NDFs', channel_names{ch}));
        set(responseFig, 'Color', 'w');
        figureHandles{end + 1, 1} = responseFig;

        responseTiles{ch} = tiledlayout(responseFig, 1, 2 * n_contrastTargets, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        title(responseTiles{ch}, sprintf('World Camera Linearity and Settings | %s | All NDFs', channel_labels{ch}), ...
            'FontWeight', 'Bold');
    end

    agc_settings_by_NDF = nan(n_NDFs, n_contrastTargets, 3);

    for cc = 1:n_contrastTargets
        contrast_target = calibration_metadata.contrast_agc_targets(cc);
        responseAxes = gobjects(1, numel(channel_labels));
        settingsAxes = gobjects(1, numel(channel_labels));
        response_by_NDF = nan(n_NDFs, n_settingsLevels, numel(channel_labels));

        for ch = 1:numel(channel_labels)
            responseAxes(ch) = nexttile(responseTiles{ch}, 2 * cc - 1);
            settingsAxes(ch) = nexttile(responseTiles{ch}, 2 * cc);
            hold(responseAxes(ch), 'on');
        end

        for nn = 1:n_NDFs
            NDF = NDFs(nn);

            if(options.verbose)
                fprintf("World Camera Linearity | NDF (%d/%d): %.3f | Contrast Target (%d/%d)\n", ...
                    nn, n_NDFs, NDF, cc, n_contrastTargets);
            end

            settings_level_measurements = reshape(measurements(nn, cc, :, :), [n_settingsLevels, n_measurements]);
            settings_levels_averaged = average_settings_levels(settings_level_measurements);
            response_by_NDF(nn, :, :) = reshape(settings_levels_averaged, [1, n_settingsLevels, numel(channel_labels)]);
            agc_settings_by_NDF(nn, cc, :) = extract_first_world_settings(settings_level_measurements);

            for ch = 1:numel(channel_labels)
                plot(responseAxes(ch), settings_levels, settings_levels_averaged(:, ch), '-o', ...
                    'Color', ndf_line_colors(nn, :), ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('NDF %.3g', NDF));
            end
        end

        contrast_agc_settings = reshape(agc_settings_by_NDF(:, cc, :), [n_NDFs, 3]);

        for ch = 1:numel(channel_labels)
            title(responseAxes(ch), ["Raw Data"; sprintf("Contrast Target: %.2f", contrast_target)], ...
                'FontWeight', 'Bold', 'FontSize', title_font_size);
            xlabel(responseAxes(ch), 'Settings Level', 'FontSize', label_font_size);
            ylabel(responseAxes(ch), 'Mean Response', 'FontSize', label_font_size);
            responseAxes(ch).FontSize = tick_font_size;
            ylim(responseAxes(ch), [0, 255]);
            grid(responseAxes(ch), 'on');
            legend(responseAxes(ch), 'Location', 'northeastoutside', 'FontSize', legend_font_size);
            hold(responseAxes(ch), 'off');

            fitResult = fit_agc_response_model(settings_levels, response_by_NDF(:, :, ch), contrast_agc_settings);
            rawFitFig = figure('Name', sprintf('Raw_vs_Fit_%s_Contrast_Target_%0.2f', channel_names{ch}, contrast_target));
            set(rawFitFig, 'Color', 'w');
            set(rawFitFig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.90, 0.85]);
            figureHandles{end + 1, 1} = rawFitFig;
            rawFitTiles = tiledlayout(rawFitFig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
            title(rawFitTiles, sprintf('Raw vs Fit | %s | Contrast AGC Target: %.2f', channel_labels{ch}, contrast_target), ...
                'FontWeight', 'Bold');
            plot_raw_vs_agc_fit_by_NDF(rawFitTiles, settings_levels, response_by_NDF(:, :, ch), ...
                fitResult.direct_response, fitResult.predicted_response, NDFs, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);
            rawFitFig.UserData = fitResult;

            if(options.verbose)
                print_agc_fit_parameter_table(channel_labels{ch}, contrast_target, NDFs, ...
                    contrast_agc_settings, fitResult.direct_parameters, fitResult.agc_predicted_parameters);
            end

            plot_camera_settings_by_NDF(settingsAxes(ch), NDFs, contrast_agc_settings, contrast_target, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);
        end
        drawnow;
    end

    if(options.verbose)
        print_camera_settings_table(NDFs, calibration_metadata.contrast_agc_targets, agc_settings_by_NDF);
    end
end


function settings_levels_averaged = average_settings_levels(settings_levels_measurements)
% Internal helper to average settings levels.
%
% Syntax:
%   settings_levels_averaged = average_settings_levels(settings_levels_measurements)
%
% Description:
%   This local helper function internal helper to average settings levels within its parent workflow.
% Inputs:
%   settings_levels_measurements - Input used by the function.
%
% Outputs:
%   settings_levels_averaged - Output produced by the function.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

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
% Internal helper to extract first world settings.
%
% Syntax:
%   world_settings = extract_first_world_settings(settings_levels_measurements)
%
% Description:
%   This local helper function internal helper to extract first world settings within its parent workflow.
% Inputs:
%   settings_levels_measurements - Input used by the function.
%
% Outputs:
%   world_settings           - Output produced by the function.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

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
% Internal helper to first setting value.
%
% Syntax:
%   value = first_setting_value(settings, setting_name)
%
% Description:
%   This local helper function internal helper to first setting value within its parent workflow.
% Inputs:
%   settings                 - Input used by the function.
%   setting_name             - Input used by the function.
%
% Outputs:
%   value                    - Output produced by the function.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

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

function fitResult = fit_agc_response_model(settings_levels, response_by_NDF, agc_settings)
% Internal helper to model camera response curves from AGC settings.
%
% Syntax:
%   fitResult = fit_agc_response_model(settings_levels, response_by_NDF, agc_settings)
%
% Description:
%   Fits each measured NDF curve with a saturating response equation,
%       y = offset + span * (1 - exp(-rate * x^gamma)),
%   and then predicts those curve parameters from [Again, Dgain, Exposure].
%
% Inputs:
%   settings_levels          - Stimulus setting levels.
%   response_by_NDF          - NDF-by-settings response matrix.
%   agc_settings             - NDF-by-3 matrix [Again, Dgain, Exposure].
%
% Outputs:
%   fitResult                - Struct containing direct and AGC-predicted fits.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    n_NDFs = size(response_by_NDF, 1);
    direct_parameters = nan(n_NDFs, 4);

    for nn = 1:n_NDFs
        direct_parameters(nn, :) = fit_saturating_response(settings_levels, response_by_NDF(nn, :));
    end

    transformed_parameters = [direct_parameters(:, 1), log_positive(direct_parameters(:, 2:4))];
    valid_rows = all(isfinite(transformed_parameters), 2) & all(isfinite(agc_settings), 2) & all(agc_settings > 0, 2);

    if(sum(valid_rows) >= 4)
        design_matrix = agc_design_matrix(agc_settings(valid_rows, :));
        parameter_coefficients = design_matrix \ transformed_parameters(valid_rows, :);
        agc_predicted_transformed = agc_design_matrix(agc_settings) * parameter_coefficients;
    else
        parameter_coefficients = nan(4, 4);
        agc_predicted_transformed = transformed_parameters;
    end

    agc_predicted_parameters = [agc_predicted_transformed(:, 1), exp(agc_predicted_transformed(:, 2:4))];
    direct_response = evaluate_saturating_response(settings_levels, direct_parameters);
    predicted_response = evaluate_saturating_response(settings_levels, agc_predicted_parameters);

    fitResult.direct_parameters = direct_parameters;
    fitResult.agc_predicted_parameters = agc_predicted_parameters;
    fitResult.direct_response = direct_response;
    fitResult.predicted_response = predicted_response;
    fitResult.parameter_coefficients = parameter_coefficients;
    fitResult.response_equation = 'y = offset + span * (1 - exp(-rate * settings_level^gamma))';
    fitResult.parameter_labels = {'offset', 'span', 'rate', 'gamma'};
    fitResult.agc_predictor_labels = {'intercept', 'log10(Again)', 'log10(Dgain)', 'log10(Exposure)'};
end

function parameters = fit_saturating_response(settings_levels, response_values)
% Internal helper to fit one saturating camera-response curve.

    x = settings_levels(:);
    y = response_values(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if(numel(x) < 4)
        parameters = nan(1, 4);
        return;
    end

    x = max(x, 0);
    y = min(max(y, 0), 255);
    offset0 = max(min(y), 0);
    span0 = max(max(y) - offset0, eps);
    rate0 = 1 / max(max(x), eps);
    gamma0 = 1;
    initial_transformed = [offset0, log(span0), log(rate0), log(gamma0)];

    fit_model = @(p, x_fit) p(1) + exp(p(2)) .* ...
        (1 - exp(-exp(p(3)) .* max(x_fit, 0) .^ exp(p(4))));

    objective = @(p) mean((y - fit_model(p, x)).^2, 'omitnan') + ...
        parameter_penalty([p(1), exp(p(2)), exp(p(3)), exp(p(4))]);
    fit_options = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000);
    fitted_transformed = fminsearch(objective, initial_transformed, fit_options);

    parameters = [fitted_transformed(1), exp(fitted_transformed(2)), exp(fitted_transformed(3)), exp(fitted_transformed(4))];
end

function penalty = parameter_penalty(parameters)
% Internal helper to softly discourage implausible response parameters.

    offset = parameters(1);
    span = parameters(2);
    rate = parameters(3);
    gamma = parameters(4);

    penalty = 0;
    penalty = penalty + 1e3 * max(0, -offset).^2;
    penalty = penalty + 1e3 * max(0, offset - 255).^2;
    penalty = penalty + 1e3 * max(0, -span).^2;
    penalty = penalty + 1e3 * max(0, offset + span - 280).^2;
    penalty = penalty + 1e3 * max(0, -rate).^2;
    penalty = penalty + 1e3 * max(0, gamma - 5).^2;
end

function response = evaluate_saturating_response(settings_levels, parameters)
% Internal helper to evaluate saturating response parameters.

    x = max(settings_levels(:)', 0);
    n_curves = size(parameters, 1);
    response = nan(n_curves, numel(x));

    for ii = 1:n_curves
        offset = parameters(ii, 1);
        span = parameters(ii, 2);
        rate = parameters(ii, 3);
        gamma = parameters(ii, 4);

        if(any(~isfinite(parameters(ii, :))))
            continue;
        end

        response(ii, :) = offset + span .* (1 - exp(-rate .* x .^ gamma));
    end
end

function design_matrix = agc_design_matrix(agc_settings)
% Internal helper to build AGC predictors for curve-parameter regression.

    log_agc = log10(max(agc_settings, realmin));
    design_matrix = [ones(size(log_agc, 1), 1), log_agc];
end

function values = log_positive(values)
% Internal helper to safely log positive fit parameters.

    values = log(max(values, realmin));
end

function plot_raw_vs_agc_fit_by_NDF(rawFitTiles, settings_levels, raw_response, direct_fit_response, agc_fit_response, NDFs, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot raw camera responses against AGC-derived fits.

    n_NDFs = size(raw_response, 1);

    for nn = 1:n_NDFs
        ax = nexttile(rawFitTiles);
        hold(ax, 'on');

        plot(ax, settings_levels, raw_response(nn, :), '-o', ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, ...
            'DisplayName', 'Raw Data');
        plot(ax, settings_levels, direct_fit_response(nn, :), '-', ...
            'Color', [0.00, 0.35, 0.85], 'LineWidth', 1.5, ...
            'DisplayName', 'Data Fit');
        plot(ax, settings_levels, agc_fit_response(nn, :), '-', ...
            'Color', [0.85, 0.10, 0.10], 'LineWidth', 1.5, ...
            'DisplayName', 'AGC Fit');

        title(ax, sprintf('NDF %.3g', NDFs(nn)), 'FontWeight', 'Bold', 'FontSize', title_font_size);
        xlabel(ax, 'Settings Level', 'FontSize', label_font_size);
        ylabel(ax, 'Mean Response', 'FontSize', label_font_size);
        ax.FontSize = tick_font_size;
        ylim(ax, [0, 255]);
        grid(ax, 'on');
        legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
        hold(ax, 'off');
    end
end

function colorList = high_contrast_ndf_colors(num_colors)
% Internal helper to create visually separated line colors for NDF plots.
%
% Syntax:
%   colorList = high_contrast_ndf_colors(num_colors)
%
% Description:
%   Returns a high-contrast palette so NDF lines are easy to distinguish
%   within each single-channel world-linearity figure.
% Inputs:
%   num_colors               - Number of colors to return.
%
% Outputs:
%   colorList                - num_colors-by-3 RGB color matrix.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    baseColors = [
        0.0000, 0.4470, 0.7410;  % Blue
        0.8500, 0.3250, 0.0980;  % Vermillion
        0.4660, 0.6740, 0.1880;  % Green
        0.4940, 0.1840, 0.5560;  % Purple
        0.9290, 0.6940, 0.1250;  % Yellow
        0.3010, 0.7450, 0.9330;  % Cyan
        0.6350, 0.0780, 0.1840;  % Maroon
        0.0000, 0.0000, 0.0000;  % Black
        0.0000, 0.6200, 0.4500;  % Bluish green
        0.8000, 0.4700, 0.6500;  % Reddish purple
        0.5800, 0.4000, 0.7400;  % Violet
        0.6000, 0.6000, 0.0000;  % Olive
        ];

    if(num_colors <= size(baseColors, 1))
        colorList = baseColors(1:num_colors, :);
        return;
    end

    extra_colors = max(num_colors - size(baseColors, 1), 0);
    fallbackColors = hsv(extra_colors + 1);
    fallbackColors = fallbackColors(1:extra_colors, :);
    colorList = [baseColors; fallbackColors];
end

function plot_camera_settings_by_NDF(ax, NDFs, settings_values, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot camera settings by ndf.
%
% Syntax:
%   plot_camera_settings_by_NDF(ax, NDFs, settings_values, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   This local helper function internal helper to plot camera settings by ndf within its parent workflow.
% Inputs:
%   ax                       - Input used by the function.
%   NDFs                     - Input used by the function.
%   settings_values          - Input used by the function.
%   contrast_target          - Input used by the function.
%   title_font_size          - Input used by the function.
%   label_font_size          - Input used by the function.
%   tick_font_size           - Input used by the function.
%   legend_font_size         - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

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
    exposureMax = max(settings_values(:, 3), [], 'omitnan');
    if(~isfinite(exposureMax) || exposureMax <= 0)
        exposureMax = 1;
    end
    ylim(ax, [0, exposureMax]);

    title(ax, sprintf("Contrast Target: %.2f", contrast_target), 'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax, 'NDF', 'FontSize', label_font_size);
    ax.FontSize = tick_font_size;
    grid(ax, 'on');
    legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
    hold(ax, 'off');
end

function print_agc_fit_parameter_table(channel_label, contrast_target, NDFs, agc_settings, direct_parameters, agc_predicted_parameters)
% Internal helper to print AGC fit parameters.

    fprintf('\nAGC camera response fit | %s | Contrast Target %.3f\n', channel_label, contrast_target);
    fprintf('%s\n', repmat('-', 1, 161));
    fprintf('%8s  %10s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
        'NDF', 'Again', 'Dgain', 'Exposure', 'Offset', 'Span', 'Rate', 'Gamma', ...
        'AGC Offset', 'AGC Span', 'AGC Rate', 'AGC Gamma');
    fprintf('%s\n', repmat('-', 1, 161));

    for nn = 1:numel(NDFs)
        fprintf('%8.3f  %10.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g\n', ...
            NDFs(nn), agc_settings(nn, 1), agc_settings(nn, 2), agc_settings(nn, 3), ...
            direct_parameters(nn, 1), direct_parameters(nn, 2), direct_parameters(nn, 3), direct_parameters(nn, 4), ...
            agc_predicted_parameters(nn, 1), agc_predicted_parameters(nn, 2), ...
            agc_predicted_parameters(nn, 3), agc_predicted_parameters(nn, 4));
    end

    fprintf('%s\n\n', repmat('-', 1, 161));
end

function print_camera_settings_table(NDFs, contrast_targets, agc_settings_by_NDF)
% Internal helper to print camera settings table.
%
% Syntax:
%   print_camera_settings_table(NDFs, contrast_targets, agc_settings_by_NDF)
%
% Description:
%   This local helper function internal helper to print camera settings table within its parent workflow.
% Inputs:
%   NDFs                     - Input used by the function.
%   contrast_targets         - Input used by the function.
%   agc_settings_by_NDF      - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    fprintf('\nWorld Camera Settings by NDF\n');
    fprintf('%s\n', repmat('-', 1, 75));
    fprintf('%16s  %10s  %12s  %12s  %12s\n', ...
        'Contrast Target', 'NDF', 'Again', 'Dgain', 'Exposure');
    fprintf('%s\n', repmat('-', 1, 75));

    for cc = 1:numel(contrast_targets)
        for nn = 1:numel(NDFs)
            settings_values = reshape(agc_settings_by_NDF(nn, cc, :), [1, 3]);
            fprintf('%16.3f  %10.3f  %12.6g  %12.6g  %12.6g\n', ...
                contrast_targets(cc), NDFs(nn), settings_values(1), settings_values(2), settings_values(3));
        end

        if(cc < numel(contrast_targets))
            fprintf('%s\n', repmat('-', 1, 75));
        end
    end

    fprintf('%s\n\n', repmat('-', 1, 75));
end
