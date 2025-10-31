% Analyze linearity calibration data collected from the MS
function  analyze_ms_linearity_data(calibration_metadata, measurements, opts)
% Analyze the results of an ms linearity light logger calibration measurement (post-conversion)
%
% Syntax:
%  analyze_ms_linearity_data(calibration_metadata, measurements)
%
% Description:
%   Analyze linearity calibration data collected from the MS.
%   Illustrate several plots showing the raw sensor counts
%   from all channels of each chip at various settings/NDF
%   levels, as well as the linearity of these in comparison
%   to the predicted counts at a given NDF level.
%
% Inputs:
%   calibration_metadata  - Struct. ms_linearity substruct
%                           of the light logger metadata
%                           struct containing only
%                           ms_linearity related metadata
%
%   measurements         - Cell. Parsed recordings made from
%                          the light logger and converted to MATLAB
%                          type.
%   plotSettingLevel     - Bool. Whether to plot the counts by settigns
%   level
%
% Outputs:
%
%   NONE
%   Given the parsed and converted metadata for an ms linearity
%   calibration measurement, analyze the data and plot the MS
%   linearity across NDF levels.
%
% Inputs:
%   calibration_metadata        - Struct. Converted metadata for
%                                 the ms linearity reading
%
%   measurements                - Cell. The parsed + converted
%                                 ms linearity readings
%
%   opts                        - struct for options.
%                                    -plotSettingLevel - bool for whether to
%                                    plot all the mod settings for each
%                                    channel.
%                                   -plotAllNDF - bool for whether to plot
%                                   all NDFs(true) or just 0-5(false)
%{
    path_to_experiment = "/example/path"; 
    converted_light_logger_data = convert_light_logger_calibration_data(path_to_experiment, true, true true, true); 
    analyze_ms_linearity_data(converted_light_logger_data.metdata.ms_linearity, converted_light_logger_data.readings.ms_linearity);
%}
arguments
    calibration_metadata; % Struct representing the metadata for the ms_linearity calibration measurement
    measurements; % Parsed and converted recordings from the light logger
    opts.plotSettingLevel logical = false;
    opts.plotAllNDF logical = true;
    opts.plotIllum logical = false;
    opts.save_illum_to_MS logical = false;
end

% Save the path to CombiExperiments. We will use this as a relative
% path to find other files
combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');

% Load the minispect SPDs
spectral_sensitivity_map = containers.Map({'ASM7341', 'TSL2591'},...
    {fullfile(combiExperiments_path,'data','ASM7341_spectralSensitivity.mat');
    fullfile(combiExperiments_path,'data','TSL2591_spectralSensitivity.mat')...
    }...
    );

% Create a map for the filters used to select good indices from the resulting curves for each chip
as_chip_point_filter = @(x, y) and(and(~isinf(y), ~isinf(x)), y >= 0.25); % AS chip we want to exclude points in the mud
ts_chip_point_filter= @(x, y) and(and(~isinf(y), ~isinf(x)), y < max(y)); % TS chip we want to exclude points that are saturated
goodIdxFilterMap = containers.Map({'ASM7341', 'TSL2591'},...
    {as_chip_point_filter, ts_chip_point_filter}...
    );


% Create a map for the limits for the chips' associated curves
lim_map = containers.Map({'ASM7341', 'TSL2591'},...
    {[-1, 5], [-6, 6]}...
    );

% Initialize a map between chips and the number of channels that they have
n_channels_map = containers.Map({'ASM7341', 'TSL2591'},...
    {size(measurements{1,1}.M.v.AS, 2),...
    size(measurements{1,1}.M.v.TS, 2)...
    }...
    );


% Get the background that was modified by the settings
% scalar and shown to the MS
background = calibration_metadata.background;

% Retrieve the list of background scalars
background_scalars = calibration_metadata.background_scalars;

% Determine the number of settings levels that were exposed
n_settings_levels = numel(background_scalars);

% Save a map between the predicted and measured counts
% across NDF levels
predicted_measured_map = containers.Map( ...
    {'ASM7341', 'TSL2591'}, ...
    { ...
    {zeros(0, n_channels_map('ASM7341')), zeros(0, n_channels_map('ASM7341'))}, ...
    {zeros(0, n_channels_map('TSL2591')), zeros(0, n_channels_map('TSL2591'))} ...
    } ...
    );

% Initialize a matrix of starting/ending indices
% for each NDF for each chip
NDF_start_end_map = containers.Map( ...
    {'ASM7341', 'TSL2591'}, ...
    { zeros(numel(calibration_metadata.NDFs), 2), ...
    zeros(numel(calibration_metadata.NDFs), 2) ...
    } ...
    );

% Make a list of colors for each ND level for the conjoined plot
colorList = [
    0.6350, 0.0780, 0.1840   % Red
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330;  % Light Blue
    0, 0.4470, 0.7410;   % Blue
    0.4940, 0.1840, 0.5560;  % Purple
    ];

% add path to isetBIO CIE luminous efficiency function
addpath('~/Documents/MATLAB/toolboxes/Psychtoolbox-3/Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles');
load("T_CIE_Y2.mat");
wls_CIE_Y2 = SToWls(S_CIE_Y2); % convert to wavelength

% First, let's iterate over the chips
chips = keys(spectral_sensitivity_map);
for cc = 1:numel(chips)

    % Retrieve the name of the chip we are analyzing
    chip = chips{cc};

    % Grab the channels of the chip we are fitting
    n_detector_channels = n_channels_map(chip);

    % Retrieve the limits for this chip
    limits = lim_map(chip);

    if opts.plotSettingLevel
        % Create a tiledlayout figure we will use to show the counts by settings
        % level for this chip across NDF levels
        [rows, cols] = find_min_figsize(numel(calibration_metadata.NDFs));
        figure;
        counts_by_NDF_tiled = tiledlayout(rows, cols);
        title(counts_by_NDF_tiled, sprintf("Counts by NDF Level | C: %s", chip), 'FontWeight', 'Bold');

        % Create a tabbed figure with one tab per NDF level.
        % Each tab will have all of a given chip's channels on it for that NDF
        uif = uifigure('Name', sprintf("%s Channel Linearity Per NDF", chip));
        tab_group = uitabgroup(uif);
    end

    % Initialize a cell array that will hold the measured/predicted value by NDF
    measured_predicted_by_NDF = {};

    % Iterate over the NDF levels
    for nn = 1:numel(calibration_metadata.NDFs)

        % Retrieve the current NDF
        NDF = calibration_metadata.NDFs(nn);

        if opts.plotSettingLevel
            % Create a new tab for this NDF for measured/predicted plot
            tab = uitab(tab_group, 'Title', sprintf("NDF %.2f", NDF));

            % Define a tiled layout that will go in this tab
            [rows, cols] = find_min_figsize(n_channels_map(chip));
            measured_predicted_tiled = tiledlayout(tab, rows, cols);
        end

        % Retrieve the cal file for this NDF
        cal = calibration_metadata.cal_files{nn};

        % Get the source from the cal file
        sourceS = [380 2 201];
        sourceP_abs = cal.processedData.P_device;

        % Retrieve the wavelengths
        wls = SToWls(sourceS);

        % Reformat minispect SPDs
        minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS);

        % Initialize detector counts
        sum_detector_counts = 0;
        predictedCounts = 0;

        % Find associated detectorP_rel
        detectorP_rel = minipspectP_rels_map(chip);

        % Extract detector counts
        detector_counts = extract_detector_counts(nn, measurements, chip);

        % Average across measurements and readings
        detector_counts = squeeze(mean(detector_counts, [2, 3]));

        if opts.plotSettingLevel
            % Plot detector counts across settings levels
            counts_ax = nexttile(counts_by_NDF_tiled);
            title(counts_ax, sprintf("Averaged Counts by Settings Level | NDF %.3f", NDF));
            hold(counts_ax, 'on');

            for ch = 1:n_detector_channels
                plot(counts_ax, background_scalars, detector_counts(:, ch), "-x", 'DisplayName', sprintf("Ch%d", ch));
            end

            xlabel(counts_ax, "Settings Level");
            ylabel(counts_ax, "Averged Count");
            legend(counts_ax, 'Location', 'best');
            hold(counts_ax, 'off');
        end

        % Initialize predicted counts and illuminance
        sphereSPDs = nan(n_settings_levels, sourceS(3));
        predictedCounts = nan(n_settings_levels, n_detector_channels);
        T_CIE_Y2_resamp = interp1(wls_CIE_Y2, T_CIE_Y2, wls, 'linear', 0);

        for ss = 1:n_settings_levels
            source_settings = background * background_scalars(ss);
            sphereSPDs(ss,:) = ( (sourceP_abs*source_settings') / sourceS(2) );

            sphereIrrad(ss,:) = sphereSPDs(ss,:) .* pi;
            irradXCIE(ss,:) = sphereIrrad(ss,:) .* T_CIE_Y2_resamp';
            integratedIrradXCIE(ss) = sum(irradXCIE(ss,:));
            illum(ss) = integratedIrradXCIE(ss) * 683;

            predictedCounts(ss,:) = sphereSPDs(ss,:) * detectorP_rel;
        end

        measured = detector_counts;
        predicted = predictedCounts;

        if opts.plotSettingLevel
            for ch = 1:n_detector_channels
                measured_predicted_ax = nexttile(measured_predicted_tiled);

                plot(measured_predicted_ax, log10((predicted(:, ch)).*pi), log10(measured(:, ch)), '-x', 'DisplayName', 'Data');
                hold(measured_predicted_ax, 'on');

                p = polyfit(log10(predicted(:, ch)), log10(measured(:, ch)), 1);
                fitY = polyval(p, log10(predicted(:, ch)));
                plot(measured_predicted_ax, log10(predicted(:, ch)).*pi, fitY, '-r', 'DisplayName', 'Fit');

                plot(measured_predicted_ax, [limits(1), limits(2)], [limits(1), limits(2)], ':k', "DisplayName", "ReferenceLine");

                xlim(measured_predicted_ax, limits);
                xlabel(measured_predicted_ax, sprintf('%s predicted irradiance [log]', chip));
                ylim(measured_predicted_ax, limits);
                ylabel(measured_predicted_ax, sprintf('%s measured counts [log]', chip));

                legend(measured_predicted_ax, 'Location','best');
                title(measured_predicted_ax, sprintf('channel %d, [slope intercept] = %2.2f, %2.2f', ch, p));
                hold(measured_predicted_ax, 'off');
            end
        end

        measured_predicted_by_NDF{nn} = {measured, predicted};
        PR670_Illum_by_NDF(cc,nn,:) = illum;

    end % NDF loop

    % Plot linearity across all NDF levels for a given chip
    [rows, cols] = find_min_figsize(n_detector_channels);
    % channeels we care about, coeffs(slope, intercept)
    if cc ==1
        illum_to_MS = nan((n_detector_channels - 1), 2);
    end

    for ch = 1:n_detector_channels
        figure;
        across_NDF_channel_ax = axes;
        hold(across_NDF_channel_ax, 'on');

        if opts.plotAllNDF
            n_ndfs_to_plot = numel(calibration_metadata.NDFs);
        else
            n_ndfs_to_plot = 5;
        end

        if ~opts.plotIllum
            for nn = 1:n_ndfs_to_plot
                NDF_measured_predicted = measured_predicted_by_NDF{nn};
                measured = NDF_measured_predicted{1};
                predicted = NDF_measured_predicted{2} * pi;

                h = scatter(across_NDF_channel_ax, ...
                    log10(predicted(:, ch)), log10(measured(:, ch)), ...
                    'o', 'MarkerFaceColor', colorList(nn,:), ...
                    'DisplayName', sprintf("NDF%.1g (%.2f lux)", calibration_metadata.NDFs(nn), round(PR670_Illum_by_NDF(cc,nn,5), 2, "significant")));

                h.MarkerFaceAlpha = 0.4;
                h.MarkerEdgeAlpha = 0.4;
            end

            plot(across_NDF_channel_ax, [limits(1), limits(2)], [limits(1), limits(2)], ':k', "DisplayName", "IdentityLine");

            title(across_NDF_channel_ax, sprintf("Channel %d", ch));
            xlim(across_NDF_channel_ax, limits);
            xlabel(across_NDF_channel_ax, sprintf('%s predicted irradiance [log]', chip));
            ylim(across_NDF_channel_ax, limits);
            ylabel(across_NDF_channel_ax, sprintf('%s measured counts [log]', chip));

            legend(across_NDF_channel_ax, 'Location', 'bestoutside');
            set(gca, 'box', 'off', 'color', 'none');
            set(gcf, 'color', 'w');
        else
            for nn = 1:n_ndfs_to_plot
                NDF_measured_predicted = measured_predicted_by_NDF{nn};
                measured = NDF_measured_predicted{1};
                measured_all_NDF_this_chip(nn,:) = measured(:, ch);

                p = scatter(across_NDF_channel_ax, ...
                    squeeze(log10(PR670_Illum_by_NDF(cc,nn,:))), ...
                    squeeze(log10(measured(:, ch))), ...
                    'o', 'MarkerFaceColor', colorList(nn,:), ...
                    'DisplayName', sprintf("NDF%.1g (%.2f lux)", calibration_metadata.NDFs(nn), round(PR670_Illum_by_NDF(cc,nn,5), 2, "significant")));

                p.MarkerFaceAlpha = 0.4;
                p.MarkerEdgeAlpha = 0.4;
            end

            set(across_NDF_channel_ax, 'YLim', limits);
            set(across_NDF_channel_ax, 'XLimMode', 'auto');
            xlabel(across_NDF_channel_ax, 'log Illuminance [lux]');
            ylabel(across_NDF_channel_ax, 'log measured sensor count');
        end

        hold(across_NDF_channel_ax, 'off');

        x = squeeze(log10(PR670_Illum_by_NDF(cc,1:n_ndfs_to_plot,:)));
        y = squeeze(log10(measured_all_NDF_this_chip(1:n_ndfs_to_plot,:)));

        x = x(:);
        y = y(:);
        valid_idx = isfinite(x) & isfinite(y);
        x = x(valid_idx);
        y = y(valid_idx);

        coeffs = polyfit(x, y, 1);
        %store the coefficients for chip one, all channels except IR
        if cc ==1 && ch<10
            illum_to_MS(ch,:) = coeffs;
        end

        x_fit = linspace(min(x), max(x), 100);
        y_fit = polyval(coeffs, x_fit);

        hold(across_NDF_channel_ax, 'on');
        plot(across_NDF_channel_ax, x_fit, y_fit, 'k:', 'LineWidth', 1, "DisplayName", "Fit Line");
        legend(across_NDF_channel_ax, 'Location', 'bestoutside');
    end % channel loop
end % chip loop
if opts.save_illum_to_MS
    filename = [combiExperiments_path, '/data/PR670_illum_to_MS_fits'];
    save(filename, "illum_to_MS");
end
end % function loop

% Local function to reformat the minispect SPDs to be in the space of
% the source SPDs
function minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS)
% Next, let's iterate over the chips at this NDF level
% and reformat the minispect SPDs to be in the space of
% the source SPDs
minipspectP_rels_map = containers.Map();
chips = keys(spectral_sensitivity_map);
for cc = 1:numel(chips)
    % Retrieve the current chip
    chip = chips{cc};

    miniSpectSPDPath = spectral_sensitivity_map(chip);
    load(miniSpectSPDPath,'T');
    minispectS = WlsToS(T.wl);
    minispectP_rel = T{:,2:end};

    detectorP_rel = [];
    for jj = 1:size(minispectP_rel,2)
        detectorP_rel(:,jj) = interp1(SToWls(minispectS),minispectP_rel(:,jj),SToWls(sourceS));
    end

    % Save the new SPDs
    minipspectP_rels_map(chip) = detectorP_rel;
end

return ;
end

% Local function to find the min square figsize required to plot data
function [rows, cols] = find_min_figsize(num_plots)
% Iterate over the ints between 1 and num_plots
for ii = 1:num_plots
    rows = ii;
    cols = ii;

    % Determine if we have reached the target
    if(rows * cols >= num_plots)
        return ;
    end

end

end

% Local function to extract solely the values of a given MS
% sensor from the measurements
function counts_mat = extract_detector_counts(NDF_num, measurements, chip)
% Extract some information about the experiment
[num_NDF_levels, num_settings_levels, n_measures] = size(measurements);

% First, let's find the max number of readings we have and the number
% of channels that we have so that we can allocate a matrix
min_num_readings = inf;
n_channels = 0;
for ss = 1:num_settings_levels
    for nn = 1:n_measures

        % Extract this measurement struct
        measurement = measurements{NDF_num, ss, nn};
        measurement_counts = 0;

        % Add the number of readings
        if(chip == "ASM7341")
            measurement_counts = measurement.M.v.AS;


            % Save the total number of channels if we have not already
            if(n_channels == 0)
                n_channels = size(measurement.M.v.AS, 2);
            end

        else
            measurement_counts = measurement.M.v.TS;

            % Save the total number of channels if we have not already
            if(n_channels == 0)
                n_channels = size(measurement.M.v.TS, 2);
            end

        end

        % Calculate the number of readings
        num_readings = size(measurement_counts, 1);
        if(num_readings == 0)
            % Output a warning, because if this happens, it's evil and
            % scary and we should probably fix the light logger to make
            % this not happen, and also should know the data is synthetic.
            warning("NDF: %d | Settings Level: %d | Measurement: %d has no readings. Generating synthetic datapoint from average", NDF_num-1, ss, nn);

            % Collect the readings that are not missing.

            % First, create a set of the number of measures
            measurement_index_set = unique(1:n_measures);

            % Now, remove the definitely bad one from this set
            good_candidates_idx = setdiff(measurement_index_set, nn);
            good_candidates = measurements{NDF_num, ss, good_candidates_idx};

            % However, if our definitely bad one is number 1,
            % let's check to make sure the candidates are good themselves.
            % Missing 1 datapoint is okay, more than 1, we should error


            % First, let's go over the candidates and find the minimum
            % number of examples that they share between themselves
            min_shared_candidate_readings = inf;
            for candidate_idx_idx = 1:numel(good_candidates_idx)
                % Retrieve the index of the candidate
                candidate_idx = good_candidates_idx(candidate_idx_idx);

                % Retrieve the candidate measurement
                candidate_measurement = measurements{NDF_num, ss, candidate_idx};

                % Initialize variable to hold good candidates readings
                candidate_counts = 0;

                % Retrieve the readings for the chip
                if(chip == "ASM7341")
                    candidate_counts = candidate_measurement.M.v.AS;

                else
                    candidate_counts = candidate_measurement.M.v.TS;
                end

                % If another example has no readings, this is really bad
                % so let's error
                if(size(candidate_counts, 1) == 0)
                    error("Recording for settings level %d has multiple measurements without readings", ss);
                end

                % Otherwise, calculate the new minimum shared size
                min_shared_candidate_readings = min(min_shared_candidate_readings, size(candidate_counts, 1));

            end

            % Initialize synthetic example. We will fill this in
            % with the first good candidate and then sum the rest to this.
            % Then, we will element wise divide to make the average.
            synthetic_example = [];

            % Iterate over the candidates
            for candidate_idx_idx = 1:numel(good_candidates_idx)
                % Retrieve the index of the candidate
                candidate_idx = good_candidates_idx(candidate_idx_idx);

                % Initialize variable to hold good candidates readings
                candidate_counts = 0;

                candidate_measurement = measurements{NDF_num, ss, candidate_idx};

                % Retrieve the readings for the chip
                if(chip == "ASM7341")
                    candidate_counts = candidate_measurement.M.v.AS;

                else
                    candidate_counts = candidate_measurement.M.v.TS;
                end

                % If we are on the first index, simply save the matrix
                if(candidate_idx_idx == 1)
                    synthetic_example = candidate_counts(1:min_shared_candidate_readings, :);
                    continue ;
                end

                % Otherwise, we can start constructing the average between
                % the candidates. So, let's add the readings from this
                % candidate to the growing list
                synthetic_example = synthetic_example + candidate_counts(1:min_shared_candidate_readings, :);

            end

            % Now, take the average of all channels for the synthetic
            % example
            synthetic_example  = synthetic_example ./ numel(good_candidates_idx);

            % Assign the synthetic example to this datapoint
            measurement_counts = synthetic_example;

            % Now, we need to save this back into the measurements array
            % because we iterate over it again
            if(chip == "ASM7341")
                measurement.M.v.AS = measurement_counts;
            else
                measurement.M.v.TS = measurement_counts;
            end

            % Resave the edited measurement.
            measurements{NDF_num, ss, nn} = measurement;


        end

        % Recalculate the number of readings, now that we have
        % potentially filled in with a synthetic example
        num_readings = size(measurement_counts, 1);

        % Calculate the min number of readings
        min_num_readings = min(min_num_readings, num_readings);
    end
end


% Initialize a matrix to store the values for this NDF
counts_mat = nan(num_settings_levels, n_measures, min_num_readings, n_channels);

% Next, we will go over each measurement and extract the channels
% for the desired chip
for ss = 1:num_settings_levels
    for nn = 1:n_measures
        % Extract this measurement struct
        measurement = measurements{NDF_num, ss, nn};

        % Retrieve the readings from the MS
        readings = 0;
        if(chip == 'ASM7341')
            readings = measurement.M.v.AS;
        else
            readings = measurement.M.v.TS;
        end

        % Insert these readings into the matrix
        counts_mat(ss, nn, 1:min_num_readings, :) = readings(1:min_num_readings, :);
    end

end
return ;

end
