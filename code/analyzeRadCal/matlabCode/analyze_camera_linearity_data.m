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
    channel_labels = {'R Mean', 'G Mean', 'B Mean', 'Frame Mean'};
    channel_names = {'R', 'G', 'B', 'Frame_Mean'};
    agc_line_colors = high_contrast_ndf_colors(n_NDFs);
    responseTiles = cell(1, numel(channel_labels));
    figureHandles = {};

    for ch = 1:numel(channel_labels)
        responseFig = figure('Name', sprintf('World_Camera_Linearity_%s_All_AGC_Conditions', channel_names{ch}));
        set(responseFig, 'Color', 'w');
        figureHandles{end + 1, 1} = responseFig;

        responseTiles{ch} = tiledlayout(responseFig, 1, 2 * n_contrastTargets, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        title(responseTiles{ch}, sprintf('World Camera Linearity and Settings | %s | All AGC Conditions', channel_labels{ch}), ...
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
            if(options.verbose)
                fprintf("World Camera Linearity | AGC source row (%d/%d) | Contrast Target (%d/%d)\n", ...
                    nn, n_NDFs, cc, n_contrastTargets);
            end

            settings_level_measurements = reshape(measurements(nn, cc, :, :), [n_settingsLevels, n_measurements]);
            settings_levels_averaged = average_settings_levels(settings_level_measurements);
            response_by_NDF(nn, :, :) = reshape(settings_levels_averaged, [1, n_settingsLevels, numel(channel_labels)]);
            agc_settings_by_NDF(nn, cc, :) = extract_first_world_settings(settings_level_measurements);

        end

        contrast_agc_settings = reshape(agc_settings_by_NDF(:, cc, :), [n_NDFs, 3]);
        agc_order = sort_agc_settings(contrast_agc_settings);

        for ch = 1:numel(channel_labels)
            for oo = 1:n_NDFs
                nn = agc_order(oo);
                plot(responseAxes(ch), settings_levels, response_by_NDF(nn, :, ch), '-o', ...
                    'Color', agc_line_colors(oo, :), ...
                    'LineWidth', 1.5, 'DisplayName', agc_condition_label(contrast_agc_settings(nn, :)));
            end

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
            plot_raw_vs_agc_fit_by_AGC(rawFitTiles, settings_levels, response_by_NDF(:, :, ch), ...
                fitResult.direct_response, fitResult.predicted_response, contrast_agc_settings, agc_order, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);
            rawFitFig.UserData = fitResult;

            demo_condition_idx = agc_order(ceil(numel(agc_order) / 2));
            gammaDemoFig = figure('Name', sprintf('Gamma_Fit_Demo_%s_Contrast_Target_%0.2f', channel_names{ch}, contrast_target));
            set(gammaDemoFig, 'Color', 'w');
            set(gammaDemoFig, 'Units', 'normalized', 'Position', [0.1, 0.15, 0.75, 0.55]);
            figureHandles{end + 1, 1} = gammaDemoFig;
            plot_gamma_fit_demo(gammaDemoFig, settings_levels, response_by_NDF(demo_condition_idx, :, ch), ...
                fitResult.direct_parameters(demo_condition_idx, :), fitResult.direct_response(demo_condition_idx, :), ...
                fitResult.direct_fit_objects{demo_condition_idx}, contrast_agc_settings(demo_condition_idx, :), ...
                channel_labels{ch}, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size);

            if(strcmp(channel_labels{ch}, 'Frame Mean'))
                [again_group_idx, again_group_fixed] = largest_group_with_fixed_settings(contrast_agc_settings, 1, [2, 3]);
                paramsVsAgainFig = figure('Name', sprintf('Parameters_vs_Analog_Gain_%s_Contrast_Target_%0.2f', channel_names{ch}, contrast_target));
                set(paramsVsAgainFig, 'Color', 'w');
                set(paramsVsAgainFig, 'Units', 'normalized', 'Position', [0.15, 0.10, 0.6, 0.7]);
                figureHandles{end + 1, 1} = paramsVsAgainFig;
                plot_parameters_vs_variable(paramsVsAgainFig, contrast_agc_settings(again_group_idx, 1), ...
                    fitResult.direct_parameters(again_group_idx, :), fitResult.parameter_labels, 'Analog Gain (Again)', ...
                    calibration_metadata.NDFs(again_group_idx), ...
                    sprintf('Fitted Parameters vs Analog Gain | %s | Contrast Target %.2f | Group: Varying Again, Fixed Dgain %.3g & Exposure %.3g', ...
                        channel_labels{ch}, contrast_target, again_group_fixed(1), again_group_fixed(2)), ...
                    title_font_size, label_font_size, tick_font_size);

                [exposure_group_idx, exposure_group_fixed] = largest_group_with_fixed_settings(contrast_agc_settings, 3, [1, 2]);
                paramsVsExposureFig = figure('Name', sprintf('Parameters_vs_Exposure_%s_Contrast_Target_%0.2f', channel_names{ch}, contrast_target));
                set(paramsVsExposureFig, 'Color', 'w');
                set(paramsVsExposureFig, 'Units', 'normalized', 'Position', [0.15, 0.10, 0.6, 0.7]);
                figureHandles{end + 1, 1} = paramsVsExposureFig;
                plot_parameters_vs_variable(paramsVsExposureFig, contrast_agc_settings(exposure_group_idx, 3), ...
                    fitResult.direct_parameters(exposure_group_idx, :), fitResult.parameter_labels, 'Exposure', ...
                    calibration_metadata.NDFs(exposure_group_idx), ...
                    sprintf('Fitted Parameters vs Exposure | %s | Contrast Target %.2f | Group: Fixed Again %.3g & Dgain %.3g, Varying Exposure', ...
                        channel_labels{ch}, contrast_target, exposure_group_fixed(1), exposure_group_fixed(2)), ...
                    title_font_size, label_font_size, tick_font_size);
            end

            if(options.verbose)
                print_agc_fit_parameter_table(channel_labels{ch}, contrast_target, contrast_agc_settings, agc_order, ...
                    fitResult.direct_parameters, fitResult.agc_predicted_parameters);
            end

            plot_camera_settings_by_AGC(settingsAxes(ch), contrast_agc_settings, agc_order, contrast_target, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);
        end
        drawnow;
    end

    if(options.verbose)
        print_camera_settings_table(calibration_metadata.contrast_agc_targets, agc_settings_by_NDF);
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
    found_settings = false;

    for ii = 1:numel(settings_levels_measurements)
        measurement = settings_levels_measurements{ii};

        if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'settings'))
            continue;
        end

        settings = measurement.W.settings;
        measurement_settings = [constant_setting_value(settings, 'Again'), ...
                                constant_setting_value(settings, 'Dgain'), ...
                                constant_setting_value(settings, 'exposure')];

        if(~found_settings)
            world_settings = measurement_settings;
            found_settings = true;
        end
    end
end

function value = constant_setting_value(settings, setting_name)
% Internal helper to validate and return a constant setting value.
%
% Syntax:
%   value = constant_setting_value(settings, setting_name)
%
% Description:
%   This local helper confirms that a setting does not change within a
%   measurement before returning its first value.
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

    setting_values = double(setting_values(:));
    finite_values = setting_values(isfinite(setting_values));

    if(~isempty(finite_values))
        setting_tolerance = 1e-9 * max(1, max(abs(finite_values)));

        if(any(abs(finite_values - finite_values(1)) > setting_tolerance))
            error('analyze_camera_linearity_data:NonconstantAGCSetting', ...
                'AGC setting "%s" changes within a measurement; cannot safely use the first value.', ...
                setting_name);
        end
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
%       y = offset + span * (1 - exp(-rate * x)),
%   to obtain curve parameters at each AGC condition. The relationship
%   between those parameters and the AGC settings [Again, Dgain, Exposure]
%   is then fit with a linear (in log-log space) meta-model. Predictions
%   are synthesized by plugging AGC settings into the meta-model to obtain
%   predicted curve parameters, then evaluating the response equation with
%   those parameters.
%
%   This equation has 3 parameters, not 4: an earlier version also
%   included a "gamma" exponent on x (y = offset + span * (1 -
%   exp(-rate * x^gamma))), but rate and gamma both control the same
%   thing -- how quickly the curve bends into saturation -- so with only
%   ~10 data points per AGC condition, many (rate, gamma) pairs fit
%   almost equally well. That let gamma (and rate) bounce around
%   noisily from one AGC condition to the next instead of varying
%   smoothly. Fixing the exponent removes that redundant degree of
%   freedom, leaving offset, span, and rate as 3 well-identified
%   parameters per AGC condition.
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
    direct_parameters = nan(n_NDFs, 3);
    direct_fit_objects = cell(n_NDFs, 1);

    for nn = 1:n_NDFs
        [direct_parameters(nn, :), direct_fit_objects{nn}] = fit_saturating_response(settings_levels, response_by_NDF(nn, :));
    end

    transformed_parameters = [direct_parameters(:, 1), log_positive(direct_parameters(:, 2:3))];
    valid_rows = all(isfinite(transformed_parameters), 2) & all(isfinite(agc_settings), 2) & all(agc_settings > 0, 2);

    design_matrix = agc_design_matrix(agc_settings);

    if(sum(valid_rows) >= size(design_matrix, 2))
        parameter_coefficients = design_matrix(valid_rows, :) \ transformed_parameters(valid_rows, :);
    else
        parameter_coefficients = nan(size(design_matrix, 2), 3);
    end

    agc_predicted_transformed = design_matrix * parameter_coefficients;
    agc_predicted_parameters = transformed_to_response_parameters(agc_predicted_transformed);
    direct_response = evaluate_saturating_response(settings_levels, direct_parameters);
    predicted_response = evaluate_saturating_response(settings_levels, agc_predicted_parameters);

    fitResult.direct_parameters = direct_parameters;
    fitResult.direct_fit_objects = direct_fit_objects;
    fitResult.agc_predicted_parameters = agc_predicted_parameters;
    fitResult.direct_response = direct_response;
    fitResult.predicted_response = predicted_response;
    fitResult.parameter_coefficients = parameter_coefficients;
    fitResult.meta_model_objective = 'AGC coefficients fit by ordinary least squares regression of curve parameters on AGC settings';
    fitResult.data_fit_method = 'Curve Fitting Toolbox fit/fittype with bounded nonlinear least squares';
    fitResult.response_equation = 'y = offset + span * (1 - exp(-rate * settings_level))';
    fitResult.parameter_labels = {'offset', 'span', 'rate'};
    fitResult.agc_predictor_labels = {'intercept', 'log10(Again)', 'log10(Dgain)', 'log10(Exposure)'};
end

function [parameters, fit_object] = fit_saturating_response(settings_levels, response_values)
% Internal helper to fit one saturating camera-response curve.

    x = settings_levels(:);
    y = response_values(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    fit_object = [];

    if(numel(x) < 3)
        parameters = nan(1, 3);
        return;
    end

    x = max(x, 0);
    y = min(max(y, 0), 255);
    offset0 = max(min(y), 0);
    span0 = max(max(y) - offset0, eps);
    rate0 = 1 / max(max(x), eps);

    response_fit_type = fittype( ...
        'offset + span * (1 - exp(-rate * x))', ...
        'independent', 'x', ...
        'coefficients', {'offset', 'span', 'rate'});
    fit_options = fitoptions(response_fit_type);
    fit_options.Display = 'Off';
    fit_options.StartPoint = [offset0, span0, rate0];
    fit_options.Lower = [0, 0, 0];
    fit_options.Upper = [255, 280, Inf];

    fit_object = fit(x, y, response_fit_type, fit_options);
    parameters = [fit_object.offset, fit_object.span, fit_object.rate];
end

function response_parameters = transformed_to_response_parameters(transformed_parameters)
% Internal helper to convert transformed meta-model outputs to response parameters.

    response_parameters = [transformed_parameters(:, 1), exp(transformed_parameters(:, 2:3))];
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

        if(any(~isfinite(parameters(ii, :))))
            continue;
        end

        response(ii, :) = offset + span .* (1 - exp(-rate .* x));
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

function plot_raw_vs_agc_fit_by_AGC(rawFitTiles, settings_levels, raw_response, direct_fit_response, agc_fit_response, agc_settings, agc_order, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot raw camera responses against AGC-derived fits.

    n_NDFs = size(raw_response, 1);

    for oo = 1:n_NDFs
        nn = agc_order(oo);
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

        title(ax, agc_condition_label(agc_settings(nn, :)), 'FontWeight', 'Bold', 'FontSize', title_font_size);
        xlabel(ax, 'Settings Level', 'FontSize', label_font_size);
        ylabel(ax, 'Mean Response', 'FontSize', label_font_size);
        ax.FontSize = tick_font_size;
        ylim(ax, [0, 255]);
        grid(ax, 'on');
        legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
        hold(ax, 'off');
    end
end

function plot_gamma_fit_demo(figHandle, settings_levels, raw_response, fitted_parameters, direct_fit_response, fit_object, agc_settings, channel_label, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to walk through one AGC condition's gamma-curve fit.
%
% Description:
%   Left panel shows the raw data, the fitted curve, and the fitted
%   [offset, span, rate] parameters for a single AGC condition, alongside
%   the response equation
%       y = offset + span * (1 - exp(-rate * x)).
%   Right panel is a sanity check: the Curve Fitting Toolbox fit object
%   itself is evaluated on a fine grid and compared against the same
%   equation typed out by hand from the fitted parameters, to confirm
%   the two agree.

    offset = fitted_parameters(1);
    span = fitted_parameters(2);
    rate = fitted_parameters(3);

    x = settings_levels(:)';
    x_fine = linspace(min(x), max(x), 200);
    manual_fit_fine = offset + span .* (1 - exp(-rate .* x_fine));
    toolbox_fit_fine = fit_object(x_fine(:))';
    max_abs_diff = max(abs(manual_fit_fine - toolbox_fit_fine));

    tiles = tiledlayout(figHandle, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiles, sprintf('Gamma Response Fit Demo | %s | Contrast Target %.2f | AGC: %s', ...
        channel_label, contrast_target, agc_condition_label(agc_settings)), 'FontWeight', 'Bold');

    ax1 = nexttile(tiles);
    hold(ax1, 'on');
    plot(ax1, x, raw_response, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Raw Data');
    plot(ax1, x, direct_fit_response, '-', 'Color', [0.00, 0.35, 0.85], 'LineWidth', 2, 'DisplayName', 'Data Fit');
    title(ax1, {'Fit and Parameters', 'y = offset + span (1 - exp(-rate \cdot x))'}, ...
        'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax1, 'Settings Level', 'FontSize', label_font_size);
    ylabel(ax1, 'Mean Response', 'FontSize', label_font_size);
    ax1.FontSize = tick_font_size;
    ylim(ax1, [0, 255]);
    grid(ax1, 'on');
    legend(ax1, 'Location', 'southeast', 'FontSize', legend_font_size);
    text(ax1, 0.03, 0.97, sprintf('offset = %.2f\nspan = %.2f\nrate = %.3g', offset, span, rate), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', label_font_size, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold(ax1, 'off');

    ax2 = nexttile(tiles);
    hold(ax2, 'on');
    plot(ax2, x, raw_response, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Raw Data');
    plot(ax2, x_fine, toolbox_fit_fine, '-', 'Color', [0.00, 0.35, 0.85], 'LineWidth', 4, 'DisplayName', 'Curve Fit Object (Toolbox)');
    plot(ax2, x_fine, manual_fit_fine, '--', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2, 'DisplayName', 'Recomputed by Hand from Fitted Parameters');
    title(ax2, {'Sanity Check', 'y = offset + span (1 - exp(-rate \cdot x))', ...
        sprintf('max |formula - fit| = %.2e', max_abs_diff)}, ...
        'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax2, 'Settings Level', 'FontSize', label_font_size);
    ylabel(ax2, 'Mean Response', 'FontSize', label_font_size);
    ax2.FontSize = tick_font_size;
    ylim(ax2, [0, 255]);
    grid(ax2, 'on');
    legend(ax2, 'Location', 'southeast', 'FontSize', legend_font_size);
    text(ax2, 0.03, 0.97, sprintf('offset = %.2f\nspan = %.2f\nrate = %.3g', offset, span, rate), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', label_font_size, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold(ax2, 'off');
end

function [group_idx, fixed_values] = largest_group_with_fixed_settings(agc_settings, vary_column, fixed_columns)
% Internal helper to find the AGC conditions that hold fixed_columns
% constant while spanning the widest set of distinct values in vary_column.
%
% Description:
%   Groups AGC conditions by their (tolerance-matched) values in
%   fixed_columns, then returns the largest such group, sorted by its
%   value in vary_column. This isolates, e.g., "same Dgain and Exposure,
%   only Again varies" so a parameter-vs-AGC plot reflects a single
%   controlled variable rather than three AGC settings changing at once.

    fixed_settings = agc_settings(:, fixed_columns);
    [~, ~, group_labels] = uniquetol(fixed_settings, 1e-6, 'ByRows', true);
    group_sizes = accumarray(group_labels, 1);
    [~, best_group] = max(group_sizes);
    group_idx = find(group_labels == best_group);
    [~, sort_order] = sort(agc_settings(group_idx, vary_column));
    group_idx = group_idx(sort_order);
    fixed_values = fixed_settings(group_idx(1), :);
end

function plot_parameters_vs_variable(figHandle, x_values, direct_parameters, parameter_labels, x_label, ndf_values, plot_title, title_font_size, label_font_size, tick_font_size)
% Internal helper to plot each fitted curve parameter against one AGC variable.
%
% Description:
%   One subplot per curve parameter (offset, span, rate), showing that
%   parameter's directly-fit value against a single AGC setting (e.g.
%   Again or Exposure) for a set of conditions where the other AGC
%   settings were held fixed. Each point is labeled with its NDF value,
%   and the group size is reported in the overall title.

    n_parameters = numel(parameter_labels);
    n_points = numel(x_values);

    tiles = tiledlayout(figHandle, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiles, sprintf('%s | n = %d points', plot_title, n_points), 'FontWeight', 'Bold');

    for pp = 1:n_parameters
        ax = nexttile(tiles);
        hold(ax, 'on');
        plot(ax, x_values, direct_parameters(:, pp), '-o', ...
            'Color', [0.00, 0.35, 0.85], 'MarkerFaceColor', [0.00, 0.35, 0.85], 'MarkerSize', 7, 'LineWidth', 1.5);

        for kk = 1:n_points
            text(ax, x_values(kk), direct_parameters(kk, pp), sprintf('  NDF %.3g', ndf_values(kk)), ...
                'FontSize', tick_font_size, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end

        title(ax, parameter_labels{pp}, 'FontWeight', 'Bold', 'FontSize', title_font_size);
        xlabel(ax, x_label, 'FontSize', label_font_size);
        ylabel(ax, parameter_labels{pp}, 'FontSize', label_font_size);
        ax.FontSize = tick_font_size;
        ax.XScale = 'log';
        grid(ax, 'on');
        hold(ax, 'off');
    end
end

function agc_order = sort_agc_settings(agc_settings)
% Internal helper to order response curves by measured AGC values.

    agc_scale = agc_exposure_gain_product(agc_settings);
    sort_values = [agc_scale(:), agc_settings];
    sort_values(~isfinite(sort_values)) = inf;
    [~, agc_order] = sortrows(sort_values);
end

function agc_scale = agc_exposure_gain_product(agc_settings)
% Internal helper to combine gain and exposure into one ordering value.

    agc_scale = prod(agc_settings, 2);
end

function label = agc_condition_label(agc_values)
% Internal helper to label plots by AGC values.

    label = sprintf('Again %.3g | Dgain %.3g | Exp %.3g', ...
        agc_values(1), agc_values(2), agc_values(3));
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

function plot_camera_settings_by_AGC(ax, settings_values, agc_order, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot camera settings by sorted AGC condition.
%
% Syntax:
%   plot_camera_settings_by_AGC(ax, settings_values, agc_order, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   This local helper function plots camera settings after sorting curves by
%   their AGC values.
% Inputs:
%   ax                       - Input used by the function.
%   settings_values          - Input used by the function.
%   agc_order                - Row order sorted by AGC values.
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
    sorted_settings = settings_values(agc_order, :);
    agc_condition_idx = 1:size(sorted_settings, 1);

    yyaxis(ax, 'left');
    plot(ax, agc_condition_idx, sorted_settings(:, 1), '-o', 'Color', [0.00, 0.35, 0.75], ...
        'LineWidth', 1.5, 'DisplayName', 'Analogue Gain');
    plot(ax, agc_condition_idx, sorted_settings(:, 2), '-o', 'Color', [0.65, 0.20, 0.75], ...
        'LineWidth', 1.5, 'DisplayName', 'Digital Gain');
    ylabel(ax, 'Gain', 'FontSize', label_font_size);

    yyaxis(ax, 'right');
    plot(ax, agc_condition_idx, sorted_settings(:, 3), '-o', 'Color', [0.95, 0.55, 0.10], ...
        'LineWidth', 1.5, 'DisplayName', 'Exposure');
    ylabel(ax, 'Exposure', 'FontSize', label_font_size);
    exposureMax = max(sorted_settings(:, 3), [], 'omitnan');
    if(~isfinite(exposureMax) || exposureMax <= 0)
        exposureMax = 1;
    end
    ylim(ax, [0, exposureMax]);

    title(ax, sprintf("Contrast Target: %.2f", contrast_target), 'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax, 'AGC Condition (sorted)', 'FontSize', label_font_size);
    xticks(ax, agc_condition_idx);
    ax.FontSize = tick_font_size;
    grid(ax, 'on');
    legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
    hold(ax, 'off');
end

function print_agc_fit_parameter_table(channel_label, contrast_target, agc_settings, agc_order, direct_parameters, agc_predicted_parameters)
% Internal helper to print AGC fit parameters.

    fprintf('\nAGC camera response fit | %s | Contrast Target %.3f\n', channel_label, contrast_target);
    fprintf('%s\n', repmat('-', 1, 146));
    fprintf('%8s  %12s  %10s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
        'AGC Idx', 'AGC Product', 'Again', 'Dgain', 'Exposure', 'Offset', 'Span', 'Rate', ...
        'AGC Offset', 'AGC Span', 'AGC Rate');
    fprintf('%s\n', repmat('-', 1, 146));

    for oo = 1:numel(agc_order)
        nn = agc_order(oo);
        fprintf('%8d  %12.5g  %10.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g  %12.5g\n', ...
            oo, agc_exposure_gain_product(agc_settings(nn, :)), ...
            agc_settings(nn, 1), agc_settings(nn, 2), agc_settings(nn, 3), ...
            direct_parameters(nn, 1), direct_parameters(nn, 2), direct_parameters(nn, 3), ...
            agc_predicted_parameters(nn, 1), agc_predicted_parameters(nn, 2), agc_predicted_parameters(nn, 3));
    end

    fprintf('%s\n\n', repmat('-', 1, 146));
end

function print_camera_settings_table(contrast_targets, agc_settings_by_NDF)
% Internal helper to print camera settings table.
%
% Syntax:
%   print_camera_settings_table(contrast_targets, agc_settings_by_NDF)
%
% Description:
%   This local helper function internal helper to print camera settings table within its parent workflow.
% Inputs:
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

    fprintf('\nWorld Camera Settings by AGC Condition\n');
    fprintf('%s\n', repmat('-', 1, 75));
    fprintf('%16s  %10s  %12s  %12s  %12s  %12s\n', ...
        'Contrast Target', 'AGC Idx', 'AGC Product', 'Again', 'Dgain', 'Exposure');
    fprintf('%s\n', repmat('-', 1, 75));

    for cc = 1:numel(contrast_targets)
        settings_values_for_contrast = reshape(agc_settings_by_NDF(:, cc, :), [size(agc_settings_by_NDF, 1), 3]);
        agc_order = sort_agc_settings(settings_values_for_contrast);

        for oo = 1:numel(agc_order)
            nn = agc_order(oo);
            settings_values = reshape(agc_settings_by_NDF(nn, cc, :), [1, 3]);
            fprintf('%16.3f  %10d  %12.6g  %12.6g  %12.6g  %12.6g\n', ...
                contrast_targets(cc), oo, agc_exposure_gain_product(settings_values), ...
                settings_values(1), settings_values(2), settings_values(3));
        end

        if(cc < numel(contrast_targets))
            fprintf('%s\n', repmat('-', 1, 75));
        end
    end

    fprintf('%s\n\n', repmat('-', 1, 75));
end
