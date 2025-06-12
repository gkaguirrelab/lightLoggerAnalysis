function msLinearityFigHandles = analyze_ms_linearity_data(ms_linearity_CalibrationData,sorted_measurements,cal,NDF,whichChannelsToPlot,msLinearityFigHandles)
% Collect calibration data for the linearity of the light-sensing
% chips of the light logger (AS/TS)
%
% Syntax:
%   success = collect_ms_linearity_data(label, dropbox_savedir, CombiLED, background, settings_scalars, n_seconds, bluetooth_central)
%
% Description:
%   Collect measurements on the linearity of the sensors in the MS
%   by exposing it, connected via the light logger's Bluetooth
%   functionality, to linearly increasing static settings on the CombiLED,
%   and recording short videos at each step. These videos are then
%   uploaded to DropBox.
%
% Inputs:
%   label                 - String. Optional notes to attach
%                           to the front of the saved filenames
%
%   dropbox_savedir       - String. The path to the directory
%                           on DropBox where the files will be
%                           output
%
%   CombiLED              - Object. Serial controller for the CombiLED
%
%   background            - Vector. Represents the initial settings
%                           of the CombiLED that we will then scale.
%
%   settings_scalars      - Vector. Represents the linearly spaced
%                           scalars in order that we will loop over
%                           and apply to the background.
%
%   n_seconds             - Int. The duration of the recording at each
%                           settings level
%
%   bluetooth_central     - Object. Represents a Python module (bluetooth_control.py)
%                           with functions for bluetooth communication.
%
% Outputs:
%
%   success               - Int. Returns 1 on success, 0 on failure.
%
% Examples:
%{
    % Point the path to a copy of the converted, calibration data
    lightLoggerCalDataFilePath = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_data/LightLogger_RadCal/W1P2M3/NDF_0/LightLoggerCalibrationData.mat';
    % Load this
    load(lightLoggerCalDataFilePath,'LightLoggerCalibrationData');
    % Extract the fields we need
    CalibrationData = LightLoggerCalibrationData.CalibrationData.ms_linearity;
    cal = LightLoggerCalibrationData.CalibrationData.cal;
    NDF = LightLoggerCalibrationData.CalibrationData.ndf;
    sorted_measurements = LightLoggerCalibrationData.ParsedReadings.ms_linearity;
    % Call the analysis routine
    analyze_ms_linearity_data(CalibrationData,sorted_measurements,cal,NDF);
%}

    % Save the path to CombiExperiments. We will use this as a relative
    % path to find other files
    combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');

    % Get the source from the cal file, as we need this to resample the
    % detector spectral sensitivity functions
    sourceS = cal.rawData.S; % This is the baseCal adjusted for transmitence

    % Extract information regarding the light source that was used to
    % calibrate the minispect
    sourceP_abs = cal.processedData.P_device;

    % Retrieve the wavelengths
    wls = SToWls(sourceS);

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
        {[-1, 5], [-1, 6]}...
        );

    % Initialize a map between chips and the number of channels that they have
    n_channels_map = containers.Map({'ASM7341', 'TSL2591'},...
        {size(sorted_measurements{1,1}.M.v.AS, 2),...
        size(sorted_measurements{1,1}.M.v.TS, 2)...
        }...
        );


    % For each chip, reformat the minispect SPDs to be in the space of the
    % sourceSPDs
    minipspectP_rels_map = containers.Map();
    chips = keys(lim_map);

    for ii = 1:numel(chips)
        miniSpectSPDPath = spectral_sensitivity_map(chips{ii});
        load(miniSpectSPDPath,'T');
        minispectS = WlsToS(T.wl);
        minispectP_rel = T{:,2:end};

        detectorP_rel = [];
        for jj = 1:size(minispectP_rel,2)
            detectorP_rel(:,jj) = interp1(SToWls(minispectS),minispectP_rel(:,jj),SToWls(sourceS));
        end

        minipspectP_rels_map(chips{ii}) = detectorP_rel;
    end

    % Get the background that was modified by the settings
    % scalar and shown to the MS
    background = ms_linearity_CalibrationData.background;

    % Retrieve the list of background scalars
    background_scalars = ms_linearity_CalibrationData.background_scalars;

    % The number of different settings that were shown to the MS
    n_primary_steps = numel(background_scalars);

    % Iterate over the chips
    for cc = 1:numel(chips)

        % Retrieve the name of the chip we are analyzing
        chip = chips{cc};

        % Initialize a summation variable for the detector counts
        % as we are going to average them over the reps
        sum_detector_counts = 0;

        % Initialize the predictedCounts variable to some value.
        % this is to have it in scope for use later.
        predictedCounts = 0;

        % Find the associated detectorP_rel for this chip
        detectorP_rel = minipspectP_rels_map(chip);

        % Grab the channels of the chip we are fitting;
        n_detector_channels = n_channels_map(chip);

        % Extract all of the counts from the sorted measurements
        % This gives you a settings x measurement x readings x channels
        % matrix
        detector_counts = extract_detector_counts(sorted_measurements, chip);

        % Now, let's take the mean across measurements (the number of times)
        % we exposed this settings level as well as the mean across readings
        % the number of readings at each measurement. We may later be
        % interested in the variability across the set of measures.
        % This gives you a mat of settings level x count
        detector_counts = squeeze(mean(detector_counts, [2, 3]));

        % Plot the detector counts across settings levels
        figure ;
        title(sprintf("Averaged Counts by Settings Level | NDF %.3f | C: %s", NDF, chip));
        hold on ;
        for ch = 1:n_detector_channels
            plot(background_scalars, detector_counts(:, ch), "-x", 'DisplayName', sprintf("CH%d", ch));
        end

        % Label the plot
        xlabel("Settings Level");
        ylabel("Averged Count");

        % Show the legend for the plot
        legend show ;

        % Initialize some variables to hold loop results
        sphereSPDs = nan(n_primary_steps, sourceS(3));
        predictedCounts = nan(n_primary_steps, n_detector_channels);

        % Iterate over the primary steps, that is, the number of scalars
        for ss = 1:n_primary_steps
            % Get the source settings by multiplying background
            % by the scalar value at primaryStep kk
            source_settings = background * background_scalars(ss);

            % Derive the sphereSPD for this step in units of W/m2/sr/nm. We divide
            % by the nanometer sampling given in S to cast the units as nm, as
            % opposed to (e.g.) per 2 nm.
            sphereSPDs(ss,:) = ( (sourceP_abs*source_settings')/sourceS(2) );

            % Derive the prediction of the relative counts based upon the sphereSPD
            % and the minispectP_rel.
            predictedCounts(ss,:) = sphereSPDs(ss,:)*detectorP_rel;
        end % End n_primary steps

        % Plot this chip's measured vs predicted counts

        % Find the limits for this chip
        limits = lim_map(chip);

        % Retrieve the filter function used to exclude points
        % from fitting LBF
        goodIdxFilter = goodIdxFilterMap(chip);

        % Retrieve the measured/predicted counts for this chip
        measured = detector_counts;
        predicted = predictedCounts;

        % Loop across the channels of the chip and show predicted vs measured
        for ch = 1:n_detector_channels

            if any(whichChannelsToPlot{cc} == ch)
                if isempty(msLinearityFigHandles)
                    msLinearityFigHandles = cell(numel(chips),n_detector_channels);
                end
                if isempty(msLinearityFigHandles{cc,ch})
                    msLinearityFigHandles{cc,ch} = figure();
                else
                    figure(msLinearityFigHandles{cc,ch});
                end

                % Draw a reference line
                plot([limits(1),limits(2)],[limits(1),limits(2)],':k', "DisplayName", "ReferenceLine");
                hold on;

                % Log transform the measured and predicted counts
                predicted_channel_readings = log10(predicted(:,ch));
                measured_channel_readings = log10(measured(:, ch));

                % Fit a linear model, but only to the "good" points (i.e., those
                % that are finite and not at the ceiling or floor. We also exclude
                % the points measured using the ND6 filter, as we do not have an
                % independent measure of the spectral transmittance of these.
                p = polyfit(predicted_channel_readings, measured_channel_readings,1);
                fitY = polyval(p, predicted_channel_readings);

                % Plot the original points
                plot(predicted_channel_readings, measured_channel_readings, '-x', 'DisplayName', 'Measured');
                plot(predicted_channel_readings, fitY, '-k', 'LineWidth', 1.5, 'DisplayName', "Fit");

                % Pretty up the plot
                xlim(limits);
                xlabel(sprintf('%s predicted counts [log]', chip));
                ylim(limits);
                ylabel(sprintf('%s measured counts [log]', chip));

                legend('Location','southeast');

                title(sprintf('channel %d, [slope intercept] = %2.2f, %2.2f',ch,p));

            end % Plot this channel

        end % End channels

    end % end chips

end

% Local function to extract solely the values of a given MS
% sensor from the sorted measurement measurements
function counts_mat = extract_detector_counts(sorted_measurements, chip)
% Extract some information about the experiment
[num_settings_levels, n_measures] = size(sorted_measurements);

% Next, we will go over each measurement and extract the channels
% for the desired chip
for ss = 1:num_settings_levels
    for nn = 1:n_measures
        % Extract this measurement struct
        measurement = sorted_measurements{ss, nn};

        % Next, we will extract the desired chips' readings
        if(chip == "ASM7341")
            counts_mat(ss, nn, :, :) = measurement.M.v.AS;
            continue ;
        end
        counts_mat(ss, nn, :, :) = measurement.M.v.TS;
    end
end

return ;

end