function [MS2illum_lux, diagnostics] = msCounts2Illuminance(ms_counts, options)
% Convert multispectral sensor counts to illuminance.
%
% Syntax:
%   MS2illum_lux = msCounts2Illuminance(ms_counts, options)
%   [MS2illum_lux, diagnostics] = msCounts2Illuminance(ms_counts, options)
%
% Description:
%   Converts the first eight filtered spectral channels to independent lux
%   estimates using PR670-derived calibration coefficients, then averages
%   those estimates. Diagnostic mode additionally recreates the legacy
%   nine-channel result for direct numerical comparison.
% Inputs:
%   ms_counts                - Input used by the function.
%   options                  - Name-value options. Set diagnostics to true
%                              to calculate an 8-versus-9-channel comparison,
%                              and verbose to true to print a short preview.
%
% Outputs:
%   MS2illum_lux             - Eight-channel illuminance estimate in lux.
%   diagnostics              - Raw counts, per-channel lux estimates, and
%                              8-versus-9-channel comparison results.
%
% Examples:
%{
    [MS2illum_lux, diagnostics] = msCounts2Illuminance( ...
        ms_counts, diagnostics=true, verbose=true)
%}

    arguments 
        ms_counts 
        options.force_recalc = false; 
        options.diagnostics (1, 1) logical = false;
        options.verbose (1, 1) logical = false;
        options.diagnostic_rows (1, 1) double {mustBeInteger, mustBePositive} = 5;
    end 

    % Define these variables as persistent so that we can save time 
    % and not have to load them in if we call them repeatedly 
    persistent illum_to_MS
    if(isempty(illum_to_MS) || options.force_recalc)
        % load coefficients for converting MS chip values to illuminance
        dropBoxBaseDir = getpref('combiExperiments','dropboxBaseDir');
        calDir = '/FLIC_data/LightLoggerRadCal/W1P1M1/MSLinearityMeasurements';
        MS_illumfile = [dropBoxBaseDir,calDir, '/PR670_illum_to_MS_fits.mat'];
        illum_to_MS = load(MS_illumfile, 'illum_to_MS').illum_to_MS;
    end 

    % Always use the eight filtered spectral channels. MATLAB columns 1-8
    % correspond to Python channels AS_0-AS_7; Clear and NIR are excluded.
    nChannels = 8;
    if size(ms_counts, 2) < nChannels
        error('ms_counts must contain at least 8 channels. Found %d.', size(ms_counts, 2));
    end
    if size(illum_to_MS, 1) < nChannels
        error('illum_to_MS must contain at least 8 channel fits. Found %d.', size(illum_to_MS, 1));
    end
    selected_counts = ms_counts(:, 1:nChannels);

    % Log conversion is valid only for positive values in the channels that
    % actually contribute to the illuminance estimate.
    if any(selected_counts(:) <= 0)
        error('The first 8 ms_counts channels must be strictly positive to apply log10.');
    end
    ms_log = log10(selected_counts);

    % Preallocate one per-channel illuminance estimate for every sample.
    nSamples = size(ms_log, 1);
    MS2illum_lux_by_channel = nan(nSamples, nChannels);

    % Invert the regression per channel
    for ch = 1:nChannels
        m = illum_to_MS(ch, 1);
        b = illum_to_MS(ch, 2);

        % Invert the linear relationship
        log10_illum = (ms_log(:, ch) - b) ./ m;

        % Convert back to lux
        MS2illum_lux_by_channel(:, ch) = 10.^log10_illum;
    end
    % Average the independently calibrated channels. This remains the
    % production return value regardless of whether diagnostics are enabled.
    MS2illum_lux = mean(MS2illum_lux_by_channel,2, 'omitnan');

    % Package the production calculation so callers can inspect it without
    % changing the primary return type or parsing console output.
    diagnostics = struct(...
        'available', false, ...
        'raw_counts', selected_counts, ...
        'illuminance_by_channel', MS2illum_lux_by_channel, ...
        'lux_8_channels', MS2illum_lux, ...
        'lux_9_channels', [], ...
        'difference_8_minus_9_lux', [], ...
        'maximum_absolute_difference_lux', NaN, ...
        'maximum_relative_difference', NaN, ...
        'exactly_identical', false, ...
        'numerically_identical', false, ...
        'comparison_tolerance_lux', NaN);

    % Recreate the legacy calculation only when requested. The ninth
    % channel's raw magnitude is not compared directly because every channel
    % has its own calibration slope and intercept before averaging.
    if options.diagnostics
        if size(ms_counts, 2) < 9 || size(illum_to_MS, 1) < 9
            warning(['Nine-channel diagnostics require at least 9 count ', ...
                'channels and 9 calibration rows.']);
            return
        end

        ninth_channel_counts = ms_counts(:, 9);
        if any(ninth_channel_counts <= 0)
            warning(['Nine-channel diagnostics were skipped because the ', ...
                'ninth channel contains nonpositive counts.']);
            return
        end

        ninth_channel_log10_illum = (...
            log10(ninth_channel_counts) - illum_to_MS(9, 2)) ...
            ./ illum_to_MS(9, 1);
        ninth_channel_illuminance = 10.^ninth_channel_log10_illum;
        nine_channel_illuminance_by_channel = [ ...
            MS2illum_lux_by_channel, ninth_channel_illuminance];
        lux_9_channels = mean(...
            nine_channel_illuminance_by_channel, 2, 'omitnan');
        difference_lux = MS2illum_lux - lux_9_channels;
        absolute_difference_lux = abs(difference_lux);
        relative_difference = absolute_difference_lux ...
            ./ max(abs(lux_9_channels), eps);
        comparison_scale = max([1; abs(MS2illum_lux(:)); abs(lux_9_channels(:))]);
        comparison_tolerance_lux = 1e-10 * comparison_scale;

        diagnostics.available = true;
        diagnostics.raw_counts = ms_counts(:, 1:9);
        diagnostics.illuminance_by_channel = nine_channel_illuminance_by_channel;
        diagnostics.lux_9_channels = lux_9_channels;
        diagnostics.difference_8_minus_9_lux = difference_lux;
        diagnostics.maximum_absolute_difference_lux = max(absolute_difference_lux);
        diagnostics.maximum_relative_difference = max(relative_difference);
        diagnostics.exactly_identical = isequaln(MS2illum_lux, lux_9_channels);
        diagnostics.numerically_identical = all(...
            absolute_difference_lux <= comparison_tolerance_lux);
        diagnostics.comparison_tolerance_lux = comparison_tolerance_lux;

        if options.verbose
            % Print only a small sample of the arrays, followed by
            % whole-recording statistics, so diagnostics remain readable.
            nRowsToDisplay = min(options.diagnostic_rows, nSamples);
            fprintf('\nMS COUNTS TO ILLUMINANCE DIAGNOSTICS\n');
            fprintf('=====================================\n');
            fprintf('Raw counts, first %d row(s), channels 1-9:\n', nRowsToDisplay);
            disp(diagnostics.raw_counts(1:nRowsToDisplay, :));
            fprintf('Calibrated lux by channel, first %d row(s):\n', nRowsToDisplay);
            disp(diagnostics.illuminance_by_channel(1:nRowsToDisplay, :));

            comparisonTable = table(...
                (1:nRowsToDisplay)', ...
                MS2illum_lux(1:nRowsToDisplay), ...
                lux_9_channels(1:nRowsToDisplay), ...
                difference_lux(1:nRowsToDisplay), ...
                'VariableNames', {'Sample', 'ReturnedLux8', 'LegacyLux9', ...
                'Lux8MinusLux9'});
            fprintf('Returned and legacy lux, first %d row(s):\n', nRowsToDisplay);
            disp(comparisonTable);
            fprintf('Maximum absolute difference: %.12g lux\n', ...
                diagnostics.maximum_absolute_difference_lux);
            fprintf('Maximum relative difference: %.12g\n', ...
                diagnostics.maximum_relative_difference);
            fprintf('Exactly identical: %s\n', ...
                string(diagnostics.exactly_identical));
            fprintf('Numerically identical at tolerance %.12g lux: %s\n\n', ...
                diagnostics.comparison_tolerance_lux, ...
                string(diagnostics.numerically_identical));
        end
    end
end
