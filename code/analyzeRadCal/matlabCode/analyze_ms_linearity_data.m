% Analyze linearity calibration data collected from the MS
function figureHandles = analyze_ms_linearity_data(calibration_metadata, measurements, opts)
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
    path_to_converted_data = fullfile(getpref('lightLoggerAnalysis', 'dropboxBaseDir'), "FLIC_data/LightLoggerRadCal/W1P1M1/MSLinearityMeasurements/MS_linearity_data.mat");
    converted_data = load(path_to_converted_data).MS_linearity_data;

    % MEAN RADIANCE PLOTS 
    analyze_ms_linearity_data(converted_data.metadata.ms_linearity, converted_data.readings.ms_linearity);

    % ILLUMINANCE PLOTS 
    analyze_ms_linearity_data(converted_data.metadata.ms_linearity, converted_data.readings.ms_linearity);
%}
arguments
    calibration_metadata; % Struct representing the metadata for the ms_linearity calibration measurement
    measurements; % Parsed and converted recordings from the light logger
    opts.output_path = false;
end

n_ndfs_to_plot = 4;

figureHandles = {};

% Defaults for all of the plotting we will do
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultAxesColor', 'w');
set(groot, 'DefaultFigureColor', 'w');

% Save the path to CombiExperiments. We will use this as a relative
% path to find other files
combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');

% Locate repository data relative to this file so the analysis does not
% depend on a machine-specific light_logger_analysis_path preference.
analysisDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(analysisDirectory)));
asSensitivityFile = fullfile( ...
    projectRoot, 'data', 'ASM7341_spectralSensitivity.mat');
assert(isfile(asSensitivityFile), ...
    'analyze_ms_linearity_data:MissingASSensitivity', ...
    'AS7341 spectral-sensitivity file not found: %s', ...
    asSensitivityFile);

% Only the AS7341 sensitivity calibration is stored in this repository.
spectral_sensitivity_map = containers.Map( ...
    {'ASM7341'}, {asSensitivityFile});

% Create a map for the filters used to select good indices from the resulting curves for each chip
as_chip_point_filter = @(x, y) and(and(~isinf(y), ~isinf(x)), y >= 0.25); % AS chip we want to exclude points in the mud
goodIdxFilterMap = containers.Map( ...
    {'ASM7341'}, {as_chip_point_filter});


% Create a map for the limits for the chips' associated curves
lim_map = containers.Map({'ASM7341'}, {[-1, 5]});

% Initialize a map between chips and the number of channels that they have
n_channels_map = containers.Map( ...
    {'ASM7341'}, {size(measurements{1,1}.M.v.AS, 2)});


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
    {'ASM7341'}, ...
    {{zeros(0, n_channels_map('ASM7341')), ...
    zeros(0, n_channels_map('ASM7341'))}});

% Initialize a matrix of starting/ending indices
% for each NDF for each chip
NDF_start_end_map = containers.Map( ...
    {'ASM7341'}, ...
    {zeros(numel(calibration_metadata.NDFs), 2)});

% Make a list of colors for each ND level for the conjoined plot.
colorList = ndf_color_list(numel(calibration_metadata.NDFs));

% Locate the Psychtoolbox CIE luminous-efficiency data.
matlabDirectory = fileparts(fileparts(projectRoot));
cieFilePath = fullfile(matlabDirectory, 'toolboxes', ...
    'Psychtoolbox-3', 'Psychtoolbox', ...
    'PsychColorimetricData', 'PsychColorimetricMatFiles', ...
    'T_CIE_Y2.mat');
assert(isfile(cieFilePath), ...
    'analyze_ms_linearity_data:MissingCIEData', ...
    'CIE luminous-efficiency file not found: %s', cieFilePath);
cieData = load(cieFilePath, 'T_CIE_Y2', 'S_CIE_Y2');
T_CIE_Y2 = cieData.T_CIE_Y2;
S_CIE_Y2 = cieData.S_CIE_Y2;
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

    % Initialize a cell array that will hold the measured/predicted value by NDF
    countsAndCalculatedRadiance = {};

    % Iterate over the NDF levels
    for nn = 1:numel(calibration_metadata.NDFs)

        % Retrieve the current NDF
        NDF = calibration_metadata.NDFs(nn);

        % Retrieve the cal file for this NDF
        cal = calibration_metadata.cal_files{nn};

        % Get the source from the cal file
        sourceS = cal.rawData.S;
        sourceP_abs = cal.processedData.P_device;

        % Retrieve the wavelengths
        wls = SToWls(sourceS);

        % Reformat minispect SPDs
        minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS);

        % Initialize detector counts
        sum_detector_counts = 0;
        calculatedRadiance = 0;

        % Find associated detectorP_rel
        detectorP_rel = minipspectP_rels_map(chip);

        % Extract detector counts
        detectorCounts = extract_detector_counts(nn, measurements, chip);

        % Average across measurements and readings
        detectorCounts = squeeze(mean(detectorCounts, [2, 3]));

        % Initialize predicted counts and illuminance
        sphereSPDs = nan(n_settings_levels, sourceS(3));
        calculatedRadiance = nan(n_settings_levels, n_detector_channels);
        T_CIE_Y2_resamp = interp1(wls_CIE_Y2, T_CIE_Y2, wls, 'linear', 0);

        for ss = 1:n_settings_levels
            source_settings = background * background_scalars(ss);
            sphereSPDs(ss,:) = sourceP_abs*source_settings';
            calculatedRadiance(ss,:) = sphereSPDs(ss,:) * detectorP_rel;
        end

        countsAndCalculatedRadiance{nn} = {detectorCounts, calculatedRadiance};

    end % NDF loop

    % Plot linearity across all NDF levels for a given chip
    [rows, cols] = find_min_figsize(n_detector_channels);
    % channels we care about, coeffs(slope, intercept)
    if cc ==1
        illum_to_MS = nan((n_detector_channels - 1), 2);
    end

    % Store one log-log linear fit per detector channel. Each row is a
    % channel and the two columns are [slope, intercept].
    logLogFitCoefficients = nan(n_detector_channels, 2);

    % Exclude measurements near the detector floor and the 16-bit ceiling
    % from the fit. These thresholds are applied to the original counts,
    % even though the regression itself is performed in log10 coordinates.
    fitCountLo = 10;
    fitCountHigh = 2^16 - 11;

    for ch = 1:n_detector_channels
        acrossNdfFig = figure('Name', sprintf("MS_Linearity_%s_Channel_%d", chip, ch));
        across_NDF_channel_ax = axes;
        hold(across_NDF_channel_ax, 'on');
        figureHandles{end+1,1} = acrossNdfFig;

        % Accumulate all plotted NDF levels so a single model can be fit to
        % the complete log-log response of this channel.
        measured_all_NDF_this_chip = nan(n_ndfs_to_plot, n_settings_levels);
        predicted_all_NDF_this_chip = nan(n_ndfs_to_plot, n_settings_levels);

        for nn = 1:n_ndfs_to_plot
            detectorCounts = countsAndCalculatedRadiance{nn}{1};
            calculatedRadiance = countsAndCalculatedRadiance{nn}{2};

            % Save this NDF's predicted radiance and measured counts in
            % log10 units. Rows identify NDFs and columns identify stimulus
            % settings.
            predicted_all_NDF_this_chip(nn,:) = log10(calculatedRadiance(:, ch)).';
            measured_all_NDF_this_chip(nn,:) = log10(detectorCounts(:, ch)).';

            h = scatter(across_NDF_channel_ax, ...
                predicted_all_NDF_this_chip(nn,:), measured_all_NDF_this_chip(nn,:), ...
                'o', 'MarkerFaceColor', colorList(nn,:), ...
                'DisplayName', sprintf('NDF-%g', calibration_metadata.NDFs(nn)));

            %'DisplayName', sprintf("NDF%.1g (%.2f lux)", calibration_metadata.NDFs(nn), round(PR670_Illum_by_NDF(cc,nn,5), 2, "significant")));

            h.MarkerFaceAlpha = 0.4;
            h.MarkerEdgeAlpha = 0.4;
        end

        % Combine the NDF-by-setting matrices into one pair of vectors for
        % this channel's regression.
        fitX = predicted_all_NDF_this_chip(:);
        fitY = measured_all_NDF_this_chip(:);

        % Keep finite values whose measured counts fall inside the desired
        % raw-count range. Comparing against log10 thresholds is equivalent
        % to applying the thresholds before the log transform.
        validFitIdx = isfinite(fitX) & isfinite(fitY) & ...
            fitY >= log10(fitCountLo) & fitY <= log10(fitCountHigh);

        % Fit log10(counts) = slope*log10(radiance) + intercept, then draw
        % the fitted segment only across the radiance values used to fit it.
        logLogFitCoefficients(ch,:) = polyfit(fitX(validFitIdx), fitY(validFitIdx), 1);
        fitXRange = [min(fitX(validFitIdx)), max(fitX(validFitIdx))];
        plot(across_NDF_channel_ax, fitXRange, polyval(logLogFitCoefficients(ch,:), fitXRange), ...
            '-k', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Linear fit'));

        plot(across_NDF_channel_ax, [limits(1), limits(2)], [limits(1), limits(2)], ':k', "DisplayName", "IdentityLine");

        title(across_NDF_channel_ax, sprintf("Channel %d", ch));
        xlabel(across_NDF_channel_ax, sprintf('%s predicted radiance [log W/m^2/sr]', chip));
        ylabel(across_NDF_channel_ax, sprintf('%s measured counts [log]', chip));

        legend(across_NDF_channel_ax, 'Location', 'bestoutside');
        set(gca, 'box', 'off', 'color', 'none');
        set(gcf, 'color', 'w');

        hold(across_NDF_channel_ax, 'off');


    end % channel loop

    % Persist the channel fits and their raw-count inclusion thresholds for
    % downstream conversion code.
    fitObj.coeff = logLogFitCoefficients;

    save(fullfile(projectRoot, 'derived', 'MSLinearityLogLogFits.mat'), ...
        'fitObj');
end % chip loop

end % function loop

function colorList = ndf_color_list(num_NDFs)
% Internal helper to create one plotting color per NDF level.
%
% Syntax:
%   colorList = ndf_color_list(num_NDFs)
%
% Description:
%   Returns the historical 7-color NDF palette for small calibrations and
%   interpolates through that palette when a calibration has more NDF levels.

baseColors = [
    0.6350, 0.0780, 0.1840   % Red
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330;  % Light Blue
    0, 0.4470, 0.7410;       % Blue
    0.4940, 0.1840, 0.5560;  % Purple
    ];

if(num_NDFs <= size(baseColors, 1))
    colorList = baseColors(1:num_NDFs, :);
    return;
end

basePositions = linspace(1, num_NDFs, size(baseColors, 1));
requestedPositions = 1:num_NDFs;
colorList = interp1(basePositions, baseColors, requestedPositions, 'linear');
end

% Local function to reformat the minispect SPDs to be in the space of
% the source SPDs
function minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS)
% Internal helper to reformat spds.
%
% Syntax:
%   minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS)
%
% Description:
%   This local helper function internal helper to reformat spds within its parent workflow.
% Inputs:
%   spectral_sensitivity_map - Input used by the function.
%   sourceS                  - Input used by the function.
%
% Outputs:
%   minipspectP_rels_map     - Output produced by the function.
%
% Examples:
%{
        % See analyze_ms_linearity_data.m for usage context.
%}

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
% Internal helper to find min figsize.
%
% Syntax:
%   rows, cols = find_min_figsize(num_plots)
%
% Description:
%   This local helper function internal helper to find min figsize within its parent workflow.
% Inputs:
%   num_plots                - Input used by the function.
%
% Outputs:
%   rows                     - Output produced by the function.
%   cols                     - Output produced by the function.
%
% Examples:
%{
        % See analyze_ms_linearity_data.m for usage context.
%}

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
% Internal helper to extract detector counts.
%
% Syntax:
%   counts_mat = extract_detector_counts(NDF_num, measurements, chip)
%
% Description:
%   This local helper function internal helper to extract detector counts within its parent workflow.
% Inputs:
%   NDF_num                  - Input used by the function.
%   measurements             - Input used by the function.
%   chip                     - Input used by the function.
%
% Outputs:
%   counts_mat               - Output produced by the function.
%
% Examples:
%{
    % See analyze_ms_linearity_data.m for usage context.
%}

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
        measurement_counts = measurement.M.v.AS;

        % Save the total number of channels if we have not already
        if(n_channels == 0)
            n_channels = size(measurement.M.v.AS, 2);
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
                candidate_counts = candidate_measurement.M.v.AS;

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
                candidate_counts = candidate_measurement.M.v.AS;

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
            measurement.M.v.AS = measurement_counts;

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
        readings = measurement.M.v.AS;

        % Insert these readings into the matrix
        counts_mat(ss, nn, 1:min_num_readings, :) = readings(1:min_num_readings, :);
    end

end
return ;

end
