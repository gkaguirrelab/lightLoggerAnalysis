function figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
% Analyze world-camera linearity data across AGC (gain/exposure) conditions.
%
% Syntax:
%   figureHandles = analyze_camera_linearity_data(calibration_metadata, measurements, options)
%
% Description:
%   For every combination of contrast target and color channel, this
%   function:
%     1. Extracts the raw mean-pixel response curve (response vs.
%        stimulus "settings level") for each AGC condition (NDF row),
%        along with the camera's [Again, Dgain, Exposure] hardware
%        settings for that condition.
%     2. Fits a saturating response equation to each condition's curve
%        only while analogue gain is not maxed out, then fits a
%        meta-model relating those curve parameters to analogue gain
%        (see fit_agc_response_model).
%     3. Produces a family of diagnostic figures: raw data overlays, raw
%        vs. fit comparisons, a worked fitting example with a sanity
%        check, analogue-gain parameter plots, and camera-settings
%        summaries.
%   All figures are collected and returned via figureHandles, and (when
%   options.verbose is true) text summaries are also printed to the
%   console.
%
% Inputs:
%   calibration_metadata     - Struct describing the calibration run. Must contain:
%                                 .background_scalars    - stimulus settings levels swept per AGC condition (x-axis).
%                                 .contrast_agc_targets   - contrast target value for each contrast-target index.
%                                 .NDFs                   - NDF (neutral density filter) value for each AGC condition/row.
%   measurements             - [n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements] cell array of
%                               raw measurement structs, each with fields .W.v (world-camera pixel means) and
%                               .W.settings (AGC hardware settings for that measurement).
%   options.verbose          - Logical; if true, print progress and summary tables to the console. Default true.
%
% Outputs:
%   figureHandles            - Cell array of all figure handles created by this analysis.
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

    % Give every figure created below a consistent, presentation-ready
    % appearance without repeating the same 'set' calls everywhere.
    configure_default_plot_appearance();

    [n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements] = size(measurements);
    title_font_size = 16;
    label_font_size = 14;
    tick_font_size = 12;
    legend_font_size = 12;

    settings_levels = calibration_metadata.background_scalars;
    channel_labels = {'R Mean', 'G Mean', 'B Mean', 'Frame Mean'};
    channel_names = {'R', 'G', 'B', 'Frame_Mean'};
    agc_line_colors = high_contrast_ndf_colors(n_NDFs);
    figureHandles = {};

    % One summary figure per channel, each holding a raw-data/settings
    % tile pair for every contrast target. The tiles are populated below
    % as each contrast target's data is processed.
    [figureHandles, responseTiles] = create_channel_summary_figures(figureHandles, channel_labels, channel_names, n_contrastTargets);

    % Accumulates every contrast target's AGC settings so the
    % end-of-run summary table (print_camera_settings_table) can report
    % across all of them at once.
    agc_settings_by_NDF = nan(n_NDFs, n_contrastTargets, 3);

    for cc = 1:n_contrastTargets
        contrast_target = calibration_metadata.contrast_agc_targets(cc);

        % Pull this contrast target's raw response curves and AGC
        % settings out of the measurements array, one row per NDF.
        [response_by_NDF, contrast_agc_settings] = extract_contrast_target_data( ...
            measurements, cc, n_contrastTargets, n_NDFs, n_settingsLevels, n_measurements, numel(channel_labels), options.verbose);
        agc_settings_by_NDF(:, cc, :) = reshape(contrast_agc_settings, [n_NDFs, 1, 3]);

        % Order AGC conditions from lowest to highest overall
        % gain*exposure. This order is reused for plot line colors,
        % legends, and printed tables so a given condition is easy to
        % track across every figure.
        agc_order = sort_agc_settings(contrast_agc_settings);

        for ch = 1:numel(channel_labels)
            % Grab this channel's pre-reserved tile pair from the
            % per-channel summary figure created above.
            responseAxes = nexttile(responseTiles{ch}, 2 * cc - 1);
            settingsAxes = nexttile(responseTiles{ch}, 2 * cc);

            plot_raw_response_by_AGC(responseAxes, settings_levels, response_by_NDF(:, :, ch), ...
                contrast_agc_settings, agc_order, agc_line_colors, contrast_target, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);

            % Fit the response model and generate every diagnostic figure
            % for this channel/contrast-target combination.
            figureHandles = analyze_and_plot_channel_fit(figureHandles, settings_levels, response_by_NDF(:, :, ch), ...
                contrast_agc_settings, agc_order, calibration_metadata.NDFs, channel_labels{ch}, channel_names{ch}, ...
                contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size, options.verbose);

            plot_camera_settings_by_AGC(settingsAxes, contrast_agc_settings, agc_order, contrast_target, ...
                title_font_size, label_font_size, tick_font_size, legend_font_size);
        end

        drawnow;
    end

    if(options.verbose)
        print_camera_settings_table(calibration_metadata.contrast_agc_targets, agc_settings_by_NDF);
    end
end


%% ========================================================================
%  Plot/figure setup helpers
%  ========================================================================

function configure_default_plot_appearance()
% Internal helper to set root graphics defaults shared by every figure.
%
% Syntax:
%   configure_default_plot_appearance()
%
% Description:
%   Sets a handful of 'Default...' properties on the graphics root (
%   groot) so every figure, axes, and text object created afterward in
%   this analysis inherits a consistent font size and white background,
%   without repeating the same formatting calls at every plotting site.
%
% Inputs:
%   None.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    set(groot, 'DefaultAxesFontSize', 16);
    set(groot, 'DefaultTextFontSize', 16);
    set(groot, 'DefaultAxesColor', 'w');
    set(groot, 'DefaultFigureColor', 'w');
end

function [figureHandles, responseTiles] = create_channel_summary_figures(figureHandles, channel_labels, channel_names, n_contrastTargets)
% Internal helper to create the one-figure-per-channel linearity summary
% figures, each pre-laid-out with a raw-data/settings tile pair for every
% contrast target.
%
% Syntax:
%   [figureHandles, responseTiles] = create_channel_summary_figures(figureHandles, channel_labels, channel_names, n_contrastTargets)
%
% Description:
%   Creates numel(channel_labels) figures, each containing a
%   1-by-(2*n_contrastTargets) tiled layout: one "raw data" tile and one
%   "camera settings" tile per contrast target. The tiles are created
%   empty here; the caller fills them in later via nexttile once each
%   contrast target's data has been extracted.
%
% Inputs:
%   figureHandles      - Existing figure-handle accumulator (cell array) to append to.
%   channel_labels     - Display names for each channel (e.g. 'R Mean').
%   channel_names      - Filesystem/figure-name-safe versions of channel_labels.
%   n_contrastTargets  - Number of contrast targets (controls tiled-layout width).
%
% Outputs:
%   figureHandles  - Updated figure-handle accumulator, with the new figures appended.
%   responseTiles  - Cell array of tiledlayout handles, one per channel.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    responseTiles = cell(1, numel(channel_labels));

    for ch = 1:numel(channel_labels)
        responseFig = figure('Name', sprintf('World_Camera_Linearity_%s_All_AGC_Conditions', channel_names{ch}));
        set(responseFig, 'Color', 'w');
        figureHandles{end + 1, 1} = responseFig;

        responseTiles{ch} = tiledlayout(responseFig, 1, 2 * n_contrastTargets, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        title(responseTiles{ch}, sprintf('World Camera Linearity and Settings | %s | All AGC Conditions', channel_labels{ch}), ...
            'FontWeight', 'Bold');
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
%   within each single-channel world-linearity figure. Uses a fixed,
%   hand-picked colorblind-friendly palette for the first 12 colors, then
%   falls back to an HSV-spaced palette if more colors are needed.
%
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

    % Fast path: the hand-picked palette already covers num_colors.
    if(num_colors <= size(baseColors, 1))
        colorList = baseColors(1:num_colors, :);
        return;
    end

    % Otherwise, pad out with additional colors spaced around the HSV wheel.
    extra_colors = max(num_colors - size(baseColors, 1), 0);
    fallbackColors = hsv(extra_colors + 1);
    fallbackColors = fallbackColors(1:extra_colors, :);
    colorList = [baseColors; fallbackColors];
end


%% ========================================================================
%  Raw measurement extraction helpers
%  ========================================================================

function [response_by_NDF, contrast_agc_settings] = extract_contrast_target_data(measurements, cc, n_contrastTargets, n_NDFs, n_settingsLevels, n_measurements, n_channels, verbose)
% Internal helper to pull one contrast target's raw responses and AGC
% settings out of the measurements array, one row per NDF condition.
%
% Syntax:
%   [response_by_NDF, contrast_agc_settings] = extract_contrast_target_data(measurements, cc, n_contrastTargets, n_NDFs, n_settingsLevels, n_measurements, n_channels, verbose)
%
% Description:
%   For a single contrast-target index cc, loops over every NDF (AGC)
%   condition, averages its repeated measurements at each settings level
%   (see average_settings_levels), and reads off the AGC hardware
%   settings that were in effect during that sweep (see
%   extract_first_world_settings). The averaged 8-bit responses are
%   normalized by 255 so all downstream fitting and plotting work on a
%   [0, 1] response scale. The result is one averaged, normalized
%   response curve and one [Again, Dgain, Exposure] triplet per NDF
%   condition.
%
% Inputs:
%   measurements       - [n_NDFs, n_contrastTargets, n_settingsLevels, n_measurements] cell array of raw measurement structs.
%   cc                 - Index of the contrast target to extract.
%   n_contrastTargets  - Total number of contrast targets (used only for the verbose progress message).
%   n_NDFs             - Number of NDF (AGC) conditions.
%   n_settingsLevels   - Number of stimulus settings levels swept per condition.
%   n_measurements     - Number of repeated measurements per settings level.
%   n_channels         - Number of color channels produced by average_settings_levels (currently 4: R/G/B/Frame mean).
%   verbose            - Whether to print a per-NDF progress message.
%
% Outputs:
%   response_by_NDF        - [n_NDFs, n_settingsLevels, n_channels] averaged response curves,
%                             normalized to [0, 1] (raw 8-bit values divided by 255).
%   contrast_agc_settings  - [n_NDFs, 3] matrix of [Again, Dgain, Exposure], one row per NDF.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    response_by_NDF = nan(n_NDFs, n_settingsLevels, n_channels);
    contrast_agc_settings = nan(n_NDFs, 3);

    for nn = 1:n_NDFs
        if(verbose)
            fprintf("World Camera Linearity | AGC source row (%d/%d) | Contrast Target (%d/%d)\n", ...
                nn, n_NDFs, cc, n_contrastTargets);
        end

        % This NDF condition's settings-level sweep, as a
        % [n_settingsLevels, n_measurements] cell of raw measurement structs.
        settings_level_measurements = reshape(measurements(nn, cc, :, :), [n_settingsLevels, n_measurements]);

        % Average repeated measurements down to one response curve per
        % channel, then normalize the 8-bit values to [0, 1].
        settings_levels_averaged = average_settings_levels(settings_level_measurements);
        response_by_NDF(nn, :, :) = reshape(settings_levels_averaged, [1, n_settingsLevels, n_channels]) ./ 255;

        % The AGC (Again/Dgain/Exposure) hardware settings are constant
        % across a single NDF condition's sweep, so read them once.
        contrast_agc_settings(nn, :) = extract_first_world_settings(settings_level_measurements);
    end
end

function settings_levels_averaged = average_settings_levels(settings_levels_measurements)
% Internal helper to average repeated measurements at each settings level.
%
% Syntax:
%   settings_levels_averaged = average_settings_levels(settings_levels_measurements)
%
% Description:
%   For each settings level, averages the per-frame world-camera pixel
%   means across frames (within a measurement) and then across the
%   repeated measurements at that settings level, producing one
%   [R, G, B, frame mean] value per settings level.
%
% Inputs:
%   settings_levels_measurements - [n_levels, n_measurements] cell array of raw measurement structs,
%                                   each with field .W.v = [n_frames, 4] pixel means.
%
% Outputs:
%   settings_levels_averaged - [n_levels, 4] matrix of averaged [R, G, B, frame mean] values.
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
            % Average across frames first, within this one measurement.
            settings_level_measurements(ll, mm, :) = mean(world_values, 1, 'omitnan');
        end
    end

    % Now average across the repeated measurements at each settings level.
    settings_levels_averaged = squeeze(mean(settings_level_measurements, 2, 'omitnan'));
end

function world_settings = extract_first_world_settings(settings_levels_measurements)
% Internal helper to read the constant AGC settings for one NDF condition.
%
% Syntax:
%   world_settings = extract_first_world_settings(settings_levels_measurements)
%
% Description:
%   Scans the measurements for the first one that has recorded camera
%   settings, and returns its [Again, Dgain, exposure] values. Since AGC
%   settings are held fixed across an entire NDF condition's settings-level
%   sweep, only the first available measurement needs to be read (and
%   constant_setting_value verifies it really is constant within that
%   measurement).
%
% Inputs:
%   settings_levels_measurements - Cell array of raw measurement structs for one NDF condition
%                                   (any shape; iterated linearly).
%
% Outputs:
%   world_settings           - [1, 3] vector of [Again, Dgain, exposure]; NaN entries if not found.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    world_settings = [nan, nan, nan];
    found_settings = false;

    for ii = 1:numel(settings_levels_measurements)
        measurement = settings_levels_measurements{ii};

        % Skip measurements that don't carry camera settings at all.
        if(isempty(measurement) || ~isfield(measurement, 'W') || ~isfield(measurement.W, 'settings'))
            continue;
        end

        settings = measurement.W.settings;
        measurement_settings = [constant_setting_value(settings, 'Again'), ...
                                constant_setting_value(settings, 'Dgain'), ...
                                constant_setting_value(settings, 'exposure')];

        % Only the first measurement with settings is needed; later ones
        % are still visited (to let constant_setting_value validate them)
        % but do not overwrite the first result.
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
%   measurement before returning its first value. If the setting varies
%   by more than a small numerical tolerance, it raises an error rather
%   than silently returning a potentially wrong value.
%
% Inputs:
%   settings                 - Struct of recorded camera settings for one measurement.
%   setting_name              - Field name of the setting to read (e.g. 'Again').
%
% Outputs:
%   value                    - The setting's (constant) value, or NaN if the field is missing/empty.
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
        % Allow for floating-point noise, scaled to the magnitude of the values.
        setting_tolerance = 1e-9 * max(1, max(abs(finite_values)));

        if(any(abs(finite_values - finite_values(1)) > setting_tolerance))
            error('analyze_camera_linearity_data:NonconstantAGCSetting', ...
                'AGC setting "%s" changes within a measurement; cannot safely use the first value.', ...
                setting_name);
        end
    end

    value = setting_values(1);
end


%% ========================================================================
%  AGC condition ordering and grouping utilities
%  ========================================================================

function agc_order = sort_agc_settings(agc_settings)
% Internal helper to order AGC conditions from lowest to highest overall gain.
%
% Syntax:
%   agc_order = sort_agc_settings(agc_settings)
%
% Description:
%   Ranks AGC conditions primarily by their combined gain*exposure
%   product (see agc_exposure_gain_product), breaking ties using the raw
%   [Again, Dgain, Exposure] values. Non-finite entries are sorted last.
%   This ordering is reused everywhere a consistent "low gain -> high
%   gain" condition order is needed (line colors, legends, tables).
%
% Inputs:
%   agc_settings              - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%
% Outputs:
%   agc_order                 - [n_NDFs, 1] row indices into agc_settings, sorted ascending by gain.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    agc_scale = agc_exposure_gain_product(agc_settings);
    sort_values = [agc_scale(:), agc_settings];
    sort_values(~isfinite(sort_values)) = inf;
    [~, agc_order] = sortrows(sort_values);
end

function agc_scale = agc_exposure_gain_product(agc_settings)
% Internal helper to combine gain and exposure into one ordering value.
%
% Syntax:
%   agc_scale = agc_exposure_gain_product(agc_settings)
%
% Description:
%   Multiplies [Again, Dgain, Exposure] together row-wise, giving a
%   single scalar per AGC condition that increases with overall
%   sensitivity/brightness. Used purely for sorting and display, not for
%   any modeling calculation.
%
% Inputs:
%   agc_settings              - [n, 3] matrix of [Again, Dgain, Exposure] (or a single [1, 3] row).
%
% Outputs:
%   agc_scale                 - [n, 1] row-wise product.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    agc_scale = prod(agc_settings, 2);
end

function label = agc_condition_label(agc_values)
% Internal helper to label plots by AGC values.
%
% Syntax:
%   label = agc_condition_label(agc_values)
%
% Description:
%   Formats one AGC condition's [Again, Dgain, Exposure] triplet into a
%   short human-readable string for use in legends and subplot titles.
%
% Inputs:
%   agc_values                - [1, 3] vector of [Again, Dgain, Exposure].
%
% Outputs:
%   label                     - Formatted string, e.g. 'Again 3.05 | Dgain 1 | Exp 8290'.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    label = sprintf('Again %.3g | Dgain %.3g | Exp %.3g', ...
        agc_values(1), agc_values(2), agc_values(3));
end

%% ========================================================================
%  Response curve model: per-condition fits and the AGC meta-model
%  ========================================================================

function fitResult = fit_agc_response_model(settings_levels, response_by_NDF, agc_settings, ~)
% Internal helper to model camera response curves from analogue gain.
%
% Syntax:
%   fitResult = fit_agc_response_model(settings_levels, response_by_NDF, agc_settings, ndf_values)
%
% Description:
%   Fits each measured NDF curve, expressed as a normalized response in
%   [0, 1] (raw 8-bit value / 255), with a saturating response equation
%   whose endpoints are pinned -- y = 0 at x = 0, plateau of 1:
%       y = (1 - exp(-(x / x0)^gamma))^shape.
%   This yields three curve parameters (x0, gamma, shape) at each fitted
%   condition. Only conditions where analogue gain has not reached its
%   maximum are fit; rows at max analogue gain are outside this
%   analogue-gain-only model and are reported as N/A.
%
%   Why the endpoints are pinned: earlier versions fit a free offset
%   (floor) and span (ceiling) as well. Empirically the fitted span was
%   250-255/255 in every condition (the camera output genuinely spans
%   the 8-bit range), so those two extra parameters mostly absorbed
%   noise -- offset in particular traded off against gamma at the dim
%   end of the curve, making both bounce across conditions. Pinning the
%   floor to 0 and the plateau to 1 leaves x0, gamma, and shape to model
%   where the curve rises, how sharply it accelerates, and how much the
%   lower limb bends before saturation.
%
%   The meta-model uses analogue gain only: each log-parameter is fit as
%   a three-coefficient quadratic polynomial in Again over the non-maxed
%   analogue gain regime. Digital gain and exposure are deliberately not
%   used here.
%   Digital gain and exposure are deliberately not used here.
%
% Inputs:
%   settings_levels          - Stimulus setting levels (shared x-axis for every AGC condition).
%   response_by_NDF          - [n_NDFs, n_settingsLevels] normalized ([0, 1]) response matrix,
%                               one row per AGC condition.
%   agc_settings             - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   ndf_values               - [n_NDFs, 1] NDF filter value for each condition (log10 attenuation).
%
% Outputs:
%   fitResult                - Struct containing:
%                                 .direct_parameters        - [n_NDFs, 3] per-condition [x0, gamma, shape].
%                                 .direct_fit_objects       - Cell array of the underlying cfit objects.
%                                 .agc_predicted_parameters - [n_NDFs, 3] parameters predicted by the analogue-gain model.
%                                 .direct_response          - Response curves evaluated from direct_parameters.
%                                 .predicted_response       - Response curves evaluated from agc_predicted_parameters.
%                                 .parameter_coefficients   - Polynomial coefficients for [log(x0), log(gamma), log(shape)].
%                                 .parameter_labels         - Names of the curve parameters.
%                                 .agc_predictor_labels     - Name of the analogue-gain predictor.
%                                 (plus descriptive metadata fields; see assignments below).
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    % Step 1: fit each analogue-gain-eligible condition's response curve
    % independently to obtain [x0, gamma, shape]. Maxed-out analogue gain
    % rows are intentionally left as NaN/N/A.
    n_NDFs = size(response_by_NDF, 1);
    direct_parameters = nan(n_NDFs, 3);
    direct_fit_objects = cell(n_NDFs, 1);
    analogue_gain_fit_mask = analog_gain_fit_mask(agc_settings);

    for nn = 1:n_NDFs
        if(~analogue_gain_fit_mask(nn))
            continue;
        end

        [direct_parameters(nn, :), direct_fit_objects{nn}] = fit_saturating_response(settings_levels, response_by_NDF(nn, :));
    end

    % Step 2: use an analogue-gain-only meta-model to predict curve
    % parameters inside the non-maxed analogue gain regime.
    [agc_predicted_parameters, meta_model] = predict_response_parameters_from_analog_gain(agc_settings, direct_parameters, analogue_gain_fit_mask);

    % Step 3: use the meta-model to predict curve parameters from
    % analogue gain, then synthesize both the "direct fit" and
    % "AGC-predicted" response curves for downstream plotting.
    direct_response = evaluate_saturating_response(settings_levels, direct_parameters);
    predicted_response = evaluate_saturating_response(settings_levels, agc_predicted_parameters);

    % Step 4: package everything (including descriptive metadata) for
    % downstream plotting, printing, and UserData attachment.
    fitResult.direct_parameters = direct_parameters;
    fitResult.direct_fit_objects = direct_fit_objects;
    fitResult.agc_predicted_parameters = agc_predicted_parameters;
    fitResult.analogue_gain_fit_mask = analogue_gain_fit_mask;
    fitResult.raw_response = response_by_NDF;
    fitResult.direct_response = direct_response;
    fitResult.predicted_response = predicted_response;
    fitResult.parameter_coefficients = meta_model.polynomial_coefficients;
    fitResult.initial_parameter_coefficients = [];
    fitResult.meta_model = meta_model;
    fitResult.meta_model_objective = 'Polynomial fit of log curve parameters over analogue gain, restricted to rows where analogue gain is not maxed out';
    fitResult.data_fit_method = 'Multi-start fminsearch on log parameters with the response constrained to [0, 1] by the pinned equation';
    fitResult.response_equation = 'y = (1 - exp(-(settings_level / x0)^gamma))^shape, y normalized to [0, 1]';
    fitResult.parameter_labels = {'x0', 'gamma', 'shape'};
    fitResult.agc_predictor_labels = {'Again (analogue gain only; maxed-out rows N/A)'};
end

function [parameters, fit_object] = fit_saturating_response(settings_levels, response_values)
% Internal helper to fit one saturating camera-response curve.
%
% Syntax:
%   [parameters, fit_object] = fit_saturating_response(settings_levels, response_values)
%
% Description:
%   Fits y = (1 - exp(-(x / x0)^gamma))^shape to one AGC condition's
%   (settings_levels, normalized response) pairs via bounded nonlinear
%   least squares, after discarding any non-finite points. The floor (0
%   at x = 0) and plateau (1) of the curve are pinned, so only x0
%   (scale), gamma (rise steepness), and shape (lower-limb curvature)
%   are free.
%
%   The starting point for the nonlinear fit is obtained from the exact
%   linearization of the model: applying the complementary log-log
%   transform to the unsaturated points gives an initial [x0, gamma]
%   estimate for the unpowered Weibull:
%       log(-log(1 - y)) = gamma * log(x) - gamma * log(x0),
%   which is a straight line in log(x). That start point is expanded
%   into several deterministic shape starts and the best least-squares
%   result is retained.
%
% Inputs:
%   settings_levels           - Stimulus settings levels (x values).
%   response_values           - Measured normalized ([0, 1]) response at each settings level.
%
% Outputs:
%   parameters                - [1, 3] fitted [x0, gamma, shape]; all NaN if too few valid points.
%   fit_object                - Struct with the fitted parameters and SSE (empty if fit not attempted).
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    % Drop any non-finite (missing/invalid) data points before fitting.
    x = settings_levels(:);
    y = response_values(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    fit_object = [];

    % Need more points than free parameters to fit meaningfully.
    if(numel(x) < 4)
        parameters = nan(1, 3);
        return;
    end

    % Clip to the physically valid ranges (non-negative x, normalized y)
    % before fitting.
    x = max(x, 0);
    y = min(max(y, 0), 1);

    % Starting point from the complementary log-log linearization, using
    % only points on the informative rising portion of the curve (the
    % transform diverges as y -> 1 and is undefined at y = 0 or x = 0).
    x0_0 = max(x) / 2;
    gamma0 = 1.2;
    rising = (y > 0.005) & (y < 0.98) & (x > 0);

    if(sum(rising) >= 2)
        linear_fit = polyfit(log(x(rising)), log(-log(1 - y(rising))), 1);
        linearized_gamma = linear_fit(1);
        linearized_x0 = exp(-linear_fit(2) / linear_fit(1));

        if(isfinite(linearized_gamma) && linearized_gamma > 0.1 && linearized_gamma < 5)
            gamma0 = linearized_gamma;
        end
        if(isfinite(linearized_x0) && linearized_x0 > 0)
            x0_0 = linearized_x0;
        end
    end

    % fminsearch is unconstrained, so optimize log-parameters and clamp
    % only through the response equation's positive parameters.
    [lower_bounds, upper_bounds] = response_parameter_bounds();
    starts = [
        x0_0, gamma0, 1.0;
        max(x) / 4, 0.8, 0.8;
        max(x) / 2, 1.5, 1.5;
        max(x), 2.5, 0.6;
        max(x) / 3, 4.0, 2.0;
        ];

    fit_options = optimset('Display', 'off', ...
                           'MaxFunEvals', 50000, ...
                           'MaxIter', 50000, ...
                           'TolFun', 1e-12, ...
                           'TolX', 1e-12);
    best_sse = inf;
    best_parameters = nan(1, 3);

    for ss = 1:size(starts, 1)
        start_parameters = min(max(starts(ss, :), lower_bounds), upper_bounds);
        start_transformed = log(start_parameters);
        objective = @(transformed_parameters) response_sse_from_transformed_parameters( ...
            transformed_parameters, x, y, lower_bounds, upper_bounds);
        [fitted_transformed, sse] = fminsearch(objective, start_transformed, fit_options);

        if(sse < best_sse)
            best_sse = sse;
            best_parameters = min(max(exp(fitted_transformed), lower_bounds), upper_bounds);
        end
    end

    parameters = best_parameters;
    fit_object = struct('x0', parameters(1), ...
                        'gamma', parameters(2), ...
                        'shape', parameters(3), ...
                        'sse', best_sse);
end

function sse = response_sse_from_transformed_parameters(transformed_parameters, x, y, lower_bounds, upper_bounds)
% Internal helper to evaluate direct-fit SSE from log-parameters.

    parameters = min(max(exp(transformed_parameters(:)'), lower_bounds), upper_bounds);
    fitted = evaluate_saturating_response(x, parameters);
    residuals = fitted(:) - y(:);
    sse = sum(residuals .^ 2, 'omitnan');
end

function fit_mask = analog_gain_fit_mask(agc_settings)
% Internal helper to identify rows before analogue gain maxes out.

    analog_gain = agc_settings(:, 1);
    finite_gain = analog_gain(isfinite(analog_gain));

    if(isempty(finite_gain))
        fit_mask = false(size(analog_gain));
        return;
    end

    max_analog_gain = max(finite_gain);
    gain_tolerance = 1e-9 * max(1, abs(max_analog_gain));
    fit_mask = isfinite(analog_gain) & (analog_gain < max_analog_gain - gain_tolerance);
end

function [predicted_parameters, meta_model] = predict_response_parameters_from_analog_gain(agc_settings, direct_parameters, fit_mask)
% Internal helper to map analogue gain to response parameters by polynomial fit.
%
% Syntax:
%   [predicted_parameters, meta_model] = predict_response_parameters_from_analog_gain(agc_settings, direct_parameters, fit_mask)
%
% Description:
%   Uses only analogue gain (Again) as the meta-model coordinate. Rows
%   where analogue gain is maxed out are left as NaN because digital gain
%   dominates those conditions and is outside this model by design. The
%   fitted parameters are log-transformed before the polynomial fit so
%   synthesized response parameters remain positive.

    analog_gain = agc_settings(:, 1);
    transformed_parameters = log_positive(direct_parameters);
    valid_rows = fit_mask(:) & isfinite(analog_gain) & all(isfinite(transformed_parameters), 2);

    predicted_parameters = nan(size(direct_parameters));
    meta_model = struct('kind', 'polynomial', ...
                        'coordinate_label', 'Again', ...
                        'coordinate', analog_gain, ...
                        'fit_mask', fit_mask, ...
                        'polynomial_degree', nan, ...
                        'polynomial_coefficients', [], ...
                        'unique_coordinate', [], ...
                        'unique_transformed_parameters', []);

    if(sum(valid_rows) < 2)
        return;
    end

    valid_coordinate = analog_gain(valid_rows);
    valid_parameters = transformed_parameters(valid_rows, :);
    [unique_coordinate, ~, unique_idx] = unique(valid_coordinate);
    unique_parameters = nan(numel(unique_coordinate), size(valid_parameters, 2));

    for uu = 1:numel(unique_coordinate)
        unique_parameters(uu, :) = mean(valid_parameters(unique_idx == uu, :), 1, 'omitnan');
    end

    polynomial_degree = min(2, numel(unique_coordinate) - 1);
    polynomial_coefficients = nan(polynomial_degree + 1, size(unique_parameters, 2));
    predicted_transformed = nan(size(transformed_parameters));

    for pp = 1:size(transformed_parameters, 2)
        polynomial_coefficients(:, pp) = polyfit(unique_coordinate, unique_parameters(:, pp), polynomial_degree)';
        predicted_transformed(fit_mask, pp) = polyval(polynomial_coefficients(:, pp)', analog_gain(fit_mask));
    end

    predicted_parameters = transformed_to_response_parameters(predicted_transformed);
    predicted_parameters(~fit_mask, :) = nan;
    meta_model.polynomial_degree = polynomial_degree;
    meta_model.polynomial_coefficients = polynomial_coefficients;
    meta_model.unique_coordinate = unique_coordinate;
    meta_model.unique_transformed_parameters = unique_parameters;
end

function response_parameters = transformed_to_response_parameters(transformed_parameters)
% Internal helper to convert transformed meta-model outputs to response parameters.
%
% Syntax:
%   response_parameters = transformed_to_response_parameters(transformed_parameters)
%
% Description:
%   Inverts the log transform applied to the curve parameters before the
%   analogue-gain polynomial meta-model (see fit_agc_response_model).
%
% Inputs:
%   transformed_parameters    - [n, 3] matrix of [log(x0), log(gamma), log(shape)].
%
% Outputs:
%   response_parameters       - [n, 3] matrix of [x0, gamma, shape] in their natural (non-log) units.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    [lower_bounds, upper_bounds] = response_parameter_bounds();
    response_parameters = min(max(exp(transformed_parameters), lower_bounds), upper_bounds);
end

function [lower_bounds, upper_bounds] = response_parameter_bounds()
% Internal helper to keep response parameters in their fitted ranges.

    lower_bounds = [1e-4, 0.05, 0.05];
    upper_bounds = [10, 12, 12];
end

function response = evaluate_saturating_response(settings_levels, parameters)
% Internal helper to evaluate saturating response parameters.
%
% Syntax:
%   response = evaluate_saturating_response(settings_levels, parameters)
%
% Description:
%   Evaluates the normalized response y = (1 - exp(-(x / x0)^gamma))^shape at
%   every settings level, for one or more parameter rows (e.g. one per
%   AGC condition). Rows containing any non-finite parameter produce a
%   row of NaN rather than an error.
%
% Inputs:
%   settings_levels           - Stimulus settings levels (x values), shared across all rows.
%   parameters                - [n_curves, 3] matrix of [x0, gamma, shape], one row per curve.
%
% Outputs:
%   response                  - [n_curves, numel(settings_levels)] evaluated normalized ([0, 1]) response curves.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    x = max(settings_levels(:)', 0);
    n_curves = size(parameters, 1);
    response = nan(n_curves, numel(x));

    for ii = 1:n_curves
        x0 = parameters(ii, 1);
        gamma = parameters(ii, 2);
        shape = parameters(ii, 3);

        if(any(~isfinite(parameters(ii, :))))
            continue;
        end

        response(ii, :) = (1 - exp(-(x ./ x0) .^ gamma)) .^ shape;
        response(ii, :) = min(max(response(ii, :), 0), 1);
    end
end

function values = log_positive(values)
% Internal helper to safely log positive fit parameters.
%
% Syntax:
%   values = log_positive(values)
%
% Description:
%   Applies log() after clamping to a small positive floor (realmin), so
%   that zero or negative parameter values (which should not occur for
%   x0/gamma, but could from a degenerate fit) do not produce -Inf
%   or complex results.
%
% Inputs:
%   values                    - Array of values expected to be positive.
%
% Outputs:
%   values                    - log() of the floored input, same size as input.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    values = log(max(values, realmin));
end


%% ========================================================================
%  Per-channel fit + diagnostic-figure orchestration
%  ========================================================================

function figureHandles = analyze_and_plot_channel_fit(figureHandles, settings_levels, response_for_channel, agc_settings, agc_order, ndf_values, channel_label, channel_name, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size, verbose)
% Internal helper to fit the response model for one channel/contrast
% target and generate all of its diagnostic figures.
%
% Syntax:
%   figureHandles = analyze_and_plot_channel_fit(figureHandles, settings_levels, response_for_channel, agc_settings, agc_order, ndf_values, channel_label, channel_name, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size, verbose)
%
% Description:
%   Ties together the curve-fitting model (fit_agc_response_model) with
%   its diagnostic plots for a single channel and contrast target:
%     1. Fit the saturating response model at every AGC condition and fit
%        the AGC-to-parameter meta-model.
%     2. Plot raw data vs. direct fit vs. AGC-predicted fit, one tile per
%        AGC condition.
%     3. Walk through one representative AGC condition's fit in detail,
%        including a sanity check of the fitted equation.
%     4. For the Frame Mean channel only, plot how the fitted parameters
%        vary across the isolated "vary Again" and "vary Exposure" AGC
%        groups (see largest_group_with_fixed_settings).
%     5. Optionally print a text table of fitted parameters.
%
% Inputs:
%   figureHandles          - Existing figure-handle accumulator (cell array) to append to.
%   settings_levels        - Stimulus settings levels (x-axis for the response curves).
%   response_for_channel   - [n_NDFs, n_settingsLevels] raw response matrix for this channel.
%   agc_settings           - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   agc_order              - AGC-sorted row order (see sort_agc_settings).
%   ndf_values             - NDF value for each row of agc_settings, used for point labels.
%   channel_label          - Display name for this channel (e.g. 'Frame Mean').
%   channel_name           - Filesystem/figure-name-safe version of channel_label.
%   contrast_target        - Contrast target value, used in titles/figure names.
%   title_font_size        - Input used by the function.
%   label_font_size        - Input used by the function.
%   tick_font_size         - Input used by the function.
%   legend_font_size       - Input used by the function.
%   verbose                - Whether to print the fitted-parameter table.
%
% Outputs:
%   figureHandles          - Updated figure-handle accumulator, with the new figures appended.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    % Step 1: fit the per-condition response curves and the AGC meta-model.
    fitResult = fit_agc_response_model(settings_levels, response_for_channel, agc_settings, ndf_values);

    % Step 2: raw data vs. direct fit vs. AGC-predicted fit, one tile per AGC condition.
    rawFitFig = figure('Name', sprintf('Raw_vs_Fit_%s_Contrast_Target_%0.2f', channel_name, contrast_target));
    set(rawFitFig, 'Color', 'w');
    set(rawFitFig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.90, 0.85]);
    figureHandles{end + 1, 1} = rawFitFig;
    rawFitTiles = tiledlayout(rawFitFig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(rawFitTiles, sprintf('Raw vs Fit | %s | Contrast AGC Target: %.2f', channel_label, contrast_target), ...
        'FontWeight', 'Bold');
    plot_raw_vs_agc_fit_by_AGC(rawFitTiles, settings_levels, response_for_channel, ...
        fitResult.direct_response, fitResult.predicted_response, agc_settings, agc_order, ...
        fitResult.analogue_gain_fit_mask, title_font_size, label_font_size, tick_font_size, legend_font_size);
    % Attach the full fit result to the figure so it can be inspected
    % later (parameters, coefficients, equation strings, etc.) without
    % needing to re-run the fit.
    rawFitFig.UserData = fitResult;

    % Step 3: walk through one representative fitted analogue-gain
    % condition, including a sanity check of the fitted equation.
    fitted_condition_idx = find(fitResult.analogue_gain_fit_mask);
    if(~isempty(fitted_condition_idx))
        demo_condition_idx = fitted_condition_idx(ceil(numel(fitted_condition_idx) / 2));
        gammaDemoFig = figure('Name', sprintf('Gamma_Fit_Demo_%s_Contrast_Target_%0.2f', channel_name, contrast_target));
        set(gammaDemoFig, 'Color', 'w');
        set(gammaDemoFig, 'Units', 'normalized', 'Position', [0.1, 0.15, 0.75, 0.55]);
        figureHandles{end + 1, 1} = gammaDemoFig;
        plot_gamma_fit_demo(gammaDemoFig, settings_levels, response_for_channel(demo_condition_idx, :), ...
            fitResult.direct_parameters(demo_condition_idx, :), fitResult.direct_response(demo_condition_idx, :), ...
            fitResult.direct_fit_objects{demo_condition_idx}, agc_settings(demo_condition_idx, :), ...
            channel_label, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size);
    end

    % Step 4: parameter-vs-analogue-gain plots, restricted to the Frame
    % Mean channel (the channel used for the analogue-gain model in practice).
    if(strcmp(channel_label, 'Frame Mean'))
        again_group_idx = find(fitResult.analogue_gain_fit_mask);
        paramsVsAgainFig = figure('Name', sprintf('Parameters_vs_Analogue_Gain_%s_Contrast_Target_%0.2f', channel_name, contrast_target));
        set(paramsVsAgainFig, 'Color', 'w');
        set(paramsVsAgainFig, 'Units', 'normalized', 'Position', [0.15, 0.10, 0.6, 0.7]);
        figureHandles{end + 1, 1} = paramsVsAgainFig;
        plot_parameters_vs_variable(paramsVsAgainFig, agc_settings(again_group_idx, 1), ...
            fitResult.direct_parameters(again_group_idx, :), fitResult.agc_predicted_parameters(again_group_idx, :), ...
            fitResult.parameter_labels, 'Analogue Gain (Again)', ...
            ndf_values(again_group_idx), ...
            sprintf('Fitted Parameters vs Analogue Gain | %s | Contrast Target %.2f | Analogue gain not maxed out', ...
                channel_label, contrast_target), ...
            title_font_size, label_font_size, tick_font_size, legend_font_size);
    end

    % Step 5: optional text summary of fitted parameters for this channel/contrast target.
    if(verbose)
        print_agc_fit_parameter_table(channel_label, contrast_target, agc_settings, agc_order, ...
            fitResult.direct_parameters, fitResult.agc_predicted_parameters);
    end
end


%% ========================================================================
%  Plotting helpers
%  ========================================================================

function plot_raw_response_by_AGC(ax, settings_levels, response_for_channel, agc_settings, agc_order, agc_line_colors, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot every AGC condition's raw response curve on one
% axes, colored and ordered by AGC value.
%
% Syntax:
%   plot_raw_response_by_AGC(ax, settings_levels, response_for_channel, agc_settings, agc_order, agc_line_colors, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   Draws one line per AGC (NDF) condition -- mean response vs. settings
%   level -- on the given axes, using agc_order to both pick the plotting
%   order and assign a consistent line color, so a given condition is
%   easy to compare across figures.
%
% Inputs:
%   ax                    - Target axes.
%   settings_levels       - Stimulus settings levels (x-axis).
%   response_for_channel  - [n_NDFs, n_settingsLevels] raw response matrix for one channel.
%   agc_settings          - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   agc_order             - AGC-sorted row order (see sort_agc_settings).
%   agc_line_colors       - [n_NDFs, 3] RGB color per plotted line, indexed by sorted position.
%   contrast_target       - Contrast target value, used in the title only.
%   title_font_size       - Input used by the function.
%   label_font_size       - Input used by the function.
%   tick_font_size        - Input used by the function.
%   legend_font_size      - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    n_NDFs = numel(agc_order);
    hold(ax, 'on');

    for oo = 1:n_NDFs
        nn = agc_order(oo);
        plot(ax, settings_levels, response_for_channel(nn, :), '-o', ...
            'Color', agc_line_colors(oo, :), ...
            'LineWidth', 1.5, 'DisplayName', agc_condition_label(agc_settings(nn, :)));
    end

    title(ax, ["Raw Data"; sprintf("Contrast Target: %.2f", contrast_target)], ...
        'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax, 'Settings Level', 'FontSize', label_font_size);
    ylabel(ax, 'Mean Response (normalized)', 'FontSize', label_font_size);
    ax.FontSize = tick_font_size;
    ylim(ax, [0, 1]);
    grid(ax, 'on');
    legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
    hold(ax, 'off');
end

function plot_raw_vs_agc_fit_by_AGC(rawFitTiles, settings_levels, raw_response, direct_fit_response, agc_fit_response, agc_settings, agc_order, fit_mask, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot raw camera responses against AGC-derived fits.
%
% Syntax:
%   plot_raw_vs_agc_fit_by_AGC(rawFitTiles, settings_levels, raw_response, direct_fit_response, agc_fit_response, agc_settings, agc_order, fit_mask, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   For each AGC condition (in agc_order), draws one tile containing:
%   the raw data, and, for analogue-gain-eligible rows, the curve
%   directly fit to that condition's own data ("Data Fit") plus the
%   analogue-gain-derived curve ("AGC Fit"). Rows where analogue gain has
%   maxed out are titled N/A and show raw data only.
%
% Inputs:
%   rawFitTiles           - Parent tiledlayout to add one tile to per AGC condition.
%   settings_levels       - Stimulus settings levels (x-axis).
%   raw_response          - [n_NDFs, n_settingsLevels] raw response matrix.
%   direct_fit_response   - [n_NDFs, n_settingsLevels] curve fit directly to each condition's own data.
%   agc_fit_response      - [n_NDFs, n_settingsLevels] curve predicted from AGC settings via the meta-model.
%   agc_settings          - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   agc_order             - AGC-sorted row order (see sort_agc_settings); also the tile order.
%   fit_mask              - Logical vector. True for rows included in the analogue-gain-only fit.
%   title_font_size       - Input used by the function.
%   label_font_size       - Input used by the function.
%   tick_font_size        - Input used by the function.
%   legend_font_size      - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    n_NDFs = size(raw_response, 1);

    for oo = 1:n_NDFs
        nn = agc_order(oo);
        ax = nexttile(rawFitTiles);
        hold(ax, 'on');

        plot(ax, settings_levels, raw_response(nn, :), '-o', ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, ...
            'DisplayName', 'Raw Data');

        if(fit_mask(nn))
            plot(ax, settings_levels, direct_fit_response(nn, :), '-', ...
                'Color', [0.00, 0.35, 0.85], 'LineWidth', 1.5, ...
                'DisplayName', 'Data Fit');
            plot(ax, settings_levels, agc_fit_response(nn, :), '-', ...
                'Color', [0.85, 0.10, 0.10], 'LineWidth', 1.5, ...
                'DisplayName', 'AGC Fit');
            title_text = agc_condition_label(agc_settings(nn, :));
        else
            title_text = ["N/A"; agc_condition_label(agc_settings(nn, :))];
        end

        title(ax, title_text, 'FontWeight', 'Bold', 'FontSize', title_font_size);
        xlabel(ax, 'Settings Level', 'FontSize', label_font_size);
        ylabel(ax, 'Mean Response (normalized)', 'FontSize', label_font_size);
        ax.FontSize = tick_font_size;
        ylim(ax, [0, 1]);
        grid(ax, 'on');
        legend(ax, 'Location', 'northeastoutside', 'FontSize', legend_font_size);
        hold(ax, 'off');
    end
end

function plot_gamma_fit_demo(figHandle, settings_levels, raw_response, fitted_parameters, direct_fit_response, ~, agc_settings, channel_label, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to walk through one AGC condition's response-curve fit.
%
% Syntax:
%   plot_gamma_fit_demo(figHandle, settings_levels, raw_response, fitted_parameters, direct_fit_response, fit_object, agc_settings, channel_label, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   Left panel shows the raw data, the fitted curve, and the fitted
%   [x0, gamma, shape] parameters for a single AGC condition, alongside the
%   normalized response equation
%       y = (1 - exp(-(x / x0)^gamma))^shape,
%   whose endpoints -- y = 0 at x = 0 and a plateau of 1 -- are pinned
%   rather than fit.
%   Right panel is a sanity check: the fitted values on the measurement
%   grid are compared against the same equation typed out by hand from
%   the fitted parameters, to confirm the two agree (max_abs_diff should
%   be ~0).
%
% Inputs:
%   figHandle             - Parent figure to lay the two panels out on.
%   settings_levels       - Stimulus settings levels (x-axis) for this AGC condition.
%   raw_response          - [1, n_settingsLevels] normalized ([0, 1]) response for this AGC condition.
%   fitted_parameters     - [1, 3] fitted [x0, gamma, shape] for this AGC condition.
%   direct_fit_response   - [1, n_settingsLevels] curve evaluated from fitted_parameters, at settings_levels.
%   fit_object            - Struct with diagnostic details for this AGC condition.
%   agc_settings          - [1, 3] vector of [Again, Dgain, Exposure] for this AGC condition.
%   channel_label         - Display name for this channel, used in the title.
%   contrast_target       - Contrast target value, used in the title.
%   title_font_size       - Input used by the function.
%   label_font_size       - Input used by the function.
%   tick_font_size        - Input used by the function.
%   legend_font_size      - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    x0 = fitted_parameters(1);
    gamma = fitted_parameters(2);
    shape = fitted_parameters(3);

    % Evaluate the fitted equation on the measurement grid and a fine x grid.
    x = settings_levels(:)';
    x_fine = linspace(min(x), max(x), 200);
    manual_fit_fine = evaluate_saturating_response(x_fine, fitted_parameters);
    manual_fit_at_measurements = evaluate_saturating_response(x, fitted_parameters);
    max_abs_diff = max(abs(manual_fit_at_measurements - direct_fit_response), [], 'omitnan');

    tiles = tiledlayout(figHandle, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiles, sprintf('Gamma Response Fit Demo | %s | Contrast Target %.2f | AGC: %s', ...
        channel_label, contrast_target, agc_condition_label(agc_settings)), 'FontWeight', 'Bold');

    % Left panel: raw data, fitted curve, and the fitted parameters.
    ax1 = nexttile(tiles);
    hold(ax1, 'on');
    plot(ax1, x, raw_response, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Raw Data');
    plot(ax1, x, direct_fit_response, '-', 'Color', [0.00, 0.35, 0.85], 'LineWidth', 2, 'DisplayName', 'Data Fit');
    title(ax1, {'Fit and Parameters', 'y = (1 - exp(-(x / x_0)^{gamma}))^{shape}'}, ...
        'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax1, 'Settings Level', 'FontSize', label_font_size);
    ylabel(ax1, 'Mean Response (normalized)', 'FontSize', label_font_size);
    ax1.FontSize = tick_font_size;
    ylim(ax1, [0, 1]);
    grid(ax1, 'on');
    legend(ax1, 'Location', 'southeast', 'FontSize', legend_font_size);
    text(ax1, 0.03, 0.97, sprintf('x0 = %.3g\ngamma = %.2f\nshape = %.2f\n(y(0) = 0, plateau = 1 pinned)', x0, gamma, shape), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', label_font_size, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold(ax1, 'off');

    % Right panel: fitted values vs. hand-computed formula, to sanity-check
    % that the equation as written matches what was actually plotted.
    ax2 = nexttile(tiles);
    hold(ax2, 'on');
    plot(ax2, x, raw_response, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'Raw Data');
    plot(ax2, x, direct_fit_response, 'o', 'Color', [0.00, 0.35, 0.85], 'MarkerFaceColor', [0.00, 0.35, 0.85], 'MarkerSize', 5, 'DisplayName', 'Data Fit');
    plot(ax2, x_fine, manual_fit_fine, '--', 'Color', [0.85, 0.10, 0.10], 'LineWidth', 2, 'DisplayName', 'Recomputed by Hand from Fitted Parameters');
    title(ax2, {'Sanity Check', 'y = (1 - exp(-(x / x_0)^{gamma}))^{shape}', ...
        sprintf('max |formula - fit| = %.2e', max_abs_diff)}, ...
        'FontWeight', 'Bold', 'FontSize', title_font_size);
    xlabel(ax2, 'Settings Level', 'FontSize', label_font_size);
    ylabel(ax2, 'Mean Response (normalized)', 'FontSize', label_font_size);
    ax2.FontSize = tick_font_size;
    ylim(ax2, [0, 1]);
    grid(ax2, 'on');
    legend(ax2, 'Location', 'southeast', 'FontSize', legend_font_size);
    text(ax2, 0.03, 0.97, sprintf('x0 = %.3g\ngamma = %.2f\nshape = %.2f', x0, gamma, shape), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', label_font_size, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold(ax2, 'off');
end

function plot_parameters_vs_variable(figHandle, x_values, direct_parameters, predicted_parameters, parameter_labels, x_label, ndf_values, plot_title, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot each fitted curve parameter against one AGC variable.
%
% Syntax:
%   plot_parameters_vs_variable(figHandle, x_values, direct_parameters, predicted_parameters, parameter_labels, x_label, ndf_values, plot_title, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   One subplot per curve parameter (x0, gamma), showing that
%   parameter's directly-fit value against a single AGC setting (e.g.
%   Again or Exposure) for a set of conditions where the other AGC
%   settings were held fixed (see largest_group_with_fixed_settings).
%   The meta-model's predicted parameter values for the same conditions
%   are overlaid, so per-condition fit scatter can be judged against the
%   smooth relationship the meta-model captures. Each point is labeled
%   with its NDF value -- important because the NDF (light level) also
%   varies within these groups, and the AGC only partially compensates
%   for it, so a jagged trend against the AGC setting alone can reflect
%   NDF variation rather than fit noise. The group size is reported in
%   the overall title.
%
% Inputs:
%   figHandle             - Parent figure to lay the subplots out on.
%   x_values              - AGC setting values for the isolated group (the plot's x-axis).
%   direct_parameters     - [n_points, n_parameters] directly-fit parameters for the isolated group.
%   predicted_parameters  - [n_points, n_parameters] meta-model-predicted parameters for the same conditions.
%   parameter_labels      - Cell array of parameter names, one per column of direct_parameters.
%   x_label               - X-axis label describing the isolated AGC setting.
%   ndf_values            - NDF value for each point, used as a per-point text label.
%   plot_title            - Overall title describing which AGC variable/group this plot isolates.
%   title_font_size       - Input used by the function.
%   label_font_size       - Input used by the function.
%   tick_font_size        - Input used by the function.
%   legend_font_size      - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    n_parameters = numel(parameter_labels);
    n_points = numel(x_values);
    [x_values, sort_order] = sort(x_values(:));
    direct_parameters = direct_parameters(sort_order, :);
    predicted_parameters = predicted_parameters(sort_order, :);
    ndf_values = ndf_values(sort_order);

    tiles = tiledlayout(figHandle, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiles, sprintf('%s | n = %d points', plot_title, n_points), 'FontWeight', 'Bold');

    for pp = 1:n_parameters
        ax = nexttile(tiles);
        hold(ax, 'on');
        plot(ax, x_values, direct_parameters(:, pp), '-o', ...
            'Color', [0.00, 0.35, 0.85], 'MarkerFaceColor', [0.00, 0.35, 0.85], ...
            'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'Direct fit');
        plot(ax, x_values, predicted_parameters(:, pp), '--s', ...
            'Color', [0.85, 0.10, 0.10], 'MarkerFaceColor', 'none', ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Meta-model prediction');

        % Label each point with its NDF value so a given physical
        % condition can be traced back from this plot.
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
        legend(ax, 'Location', 'best', 'FontSize', legend_font_size);
        hold(ax, 'off');
    end
end

function plot_camera_settings_by_AGC(ax, settings_values, agc_order, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
% Internal helper to plot camera settings by sorted AGC condition.
%
% Syntax:
%   plot_camera_settings_by_AGC(ax, settings_values, agc_order, contrast_target, title_font_size, label_font_size, tick_font_size, legend_font_size)
%
% Description:
%   Plots Analog Gain and Digital Gain (left y-axis) and Exposure (right
%   y-axis) against AGC condition index, after sorting conditions by
%   agc_order, so the reader can see how the AGC algorithm traded off
%   gain vs. exposure across the sweep of NDF conditions.
%
% Inputs:
%   ax                       - Target axes.
%   settings_values          - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   agc_order                - Row order sorted by AGC value (see sort_agc_settings).
%   contrast_target          - Contrast target value, used in the title only.
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


%% ========================================================================
%  Console reporting helpers
%  ========================================================================

function print_agc_fit_parameter_table(channel_label, contrast_target, agc_settings, agc_order, direct_parameters, agc_predicted_parameters)
% Internal helper to print a text table of fitted AGC response parameters.
%
% Syntax:
%   print_agc_fit_parameter_table(channel_label, contrast_target, agc_settings, agc_order, direct_parameters, agc_predicted_parameters)
%
% Description:
%   Prints, for one channel/contrast target, a row per AGC condition
%   (sorted by agc_order) with its AGC settings, its directly-fit curve
%   parameters, and the parameters predicted by the AGC meta-model, so
%   the two can be compared at a glance.
%
% Inputs:
%   channel_label              - Display name for this channel, used in the header.
%   contrast_target            - Contrast target value, used in the header.
%   agc_settings               - [n_NDFs, 3] matrix of [Again, Dgain, Exposure].
%   agc_order                  - AGC-sorted row order (see sort_agc_settings); also the print order.
%   direct_parameters          - [n_NDFs, 3] directly-fit [x0, gamma, shape].
%   agc_predicted_parameters   - [n_NDFs, 3] [x0, gamma, shape] predicted by the meta-model.
%
% Outputs:
%   None (writes to the console via fprintf).
%
% Examples:
%{
    % See analyze_camera_linearity_data.m for usage context.
%}

    fprintf('\nAGC camera response fit | %s | Contrast Target %.3f\n', channel_label, contrast_target);
    fprintf('%s\n', repmat('-', 1, 146));
    fprintf('%8s  %12s  %10s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
        'AGC Idx', 'AGC Product', 'Again', 'Dgain', 'Exposure', 'X0', 'Gamma', ...
        'Shape', 'Meta X0', 'Meta Gamma', 'Meta Shape');
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
% Internal helper to print a text table of AGC settings across the whole run.
%
% Syntax:
%   print_camera_settings_table(contrast_targets, agc_settings_by_NDF)
%
% Description:
%   For every contrast target, prints one row per AGC condition (sorted
%   by AGC value) with its [Again, Dgain, Exposure] settings, giving a
%   compact end-of-run summary of how the AGC algorithm behaved across
%   the whole calibration.
%
% Inputs:
%   contrast_targets           - Vector of contrast target values.
%   agc_settings_by_NDF        - [n_NDFs, n_contrastTargets, 3] matrix of [Again, Dgain, Exposure].
%
% Outputs:
%   None (writes to the console via fprintf).
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
