function  temporal_offsets_secs = analyze_phase_fit_data(CalibrationData, sorted_measurements)

% Specify the number of measurements to discard from the start of the pupil
% recording (as there is a start-up effect on the signal here)
nPupilFramesToDiscard = 5;

% Retrieve the frequencies and contrast levels used to make the measurements
% as well as the number of measurements made at each level (this should be 1 here)
contrast_levels = CalibrationData.contrast_levels;
frequencies = CalibrationData.frequencies;
n_measures = CalibrationData.n_measures;

% We will assume that the first measurement (frequency) is the slow frequency
% and the second measurement (frequency) is the high frequency
assert(numel(contrast_levels) == 1 && numel(frequencies) == 2);
assert(numel(sorted_measurements) == numel(contrast_levels) * numel(frequencies) * n_measures);

% Now, we will iterate over the measurements. We will do this to take the average phase
temporal_offsets_secs = containers.Map({'W-AS', 'W-TS', 'W-P'}, { zeros([n_measures, 1]), zeros([n_measures, 1]), zeros([n_measures, 1]) });

for nn = 1:n_measures
    % Retrieve the sensor readings for this measurement
    slow_measurement = sorted_measurements{1, 1, nn};
    fast_measurement = sorted_measurements{1, 2, nn};

    % Adjust the start times to be relative to world camera
    startTime = slow_measurement.W.t(1);
    slow_measurement.W.t = slow_measurement.W.t - startTime;
    slow_measurement.P.t = slow_measurement.P.t - startTime;
    slow_measurement.M.t = slow_measurement.M.t - startTime;
    slow_measurement.S.t = slow_measurement.S.t - startTime;

        startTime = fast_measurement.W.t(1);
    fast_measurement.W.t = fast_measurement.W.t - startTime;
    fast_measurement.P.t = fast_measurement.P.t - startTime;
    fast_measurement.M.t = fast_measurement.M.t - startTime;
    fast_measurement.S.t = fast_measurement.S.t - startTime;
    
    % Calculate the phase difference of the MS chips to the world camera
    % from the slow measurement. Calculate and store the temporal offset
    % between sensors in units of seconds.
    [world_ms_AS_phase_diff, world_fit_slow_AS, AS_fit] = calculate_phase_offset(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), ...
        slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), ...
        frequencies(1));

    world_ms_AS_temporal_offset_secs = temporal_offsets_secs('W-AS');
    world_ms_AS_temporal_offset_secs(nn) = (world_ms_AS_phase_diff / (2 * pi)) / frequencies(1);
    temporal_offsets_secs('W-AS') = world_ms_AS_temporal_offset_secs;

    [world_ms_TS_phase_diff, world_fit_slow_TS, TS_fit] = calculate_phase_offset(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), ...
        slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.TS(:, 1)), ...
        frequencies(1));

    world_ms_TS_temporal_offset_secs = temporal_offsets_secs('W-TS');
    world_ms_TS_temporal_offset_secs(nn) = (world_ms_TS_phase_diff / (2 * pi)) / frequencies(1);
    temporal_offsets_secs('W-TS') = world_ms_TS_temporal_offset_secs;

    % Calculate the phase difference between the world and pupil cameras from the fast measurements
    [world_pupil_phase_diff, world_fit_fast, pupil_fit] = calculate_phase_offset(fast_measurement.W.t, convert_to_contrast(fast_measurement.W.v), ...
        fast_measurement.P.t(nPupilFramesToDiscard:end), convert_to_contrast(fast_measurement.P.v(nPupilFramesToDiscard:end)), ...
        frequencies(2));

    world_pupil_temporal_offset_secs = temporal_offsets_secs('W-P');
    world_pupil_temporal_offset_secs(nn) = (world_pupil_phase_diff / (2 * pi)) / frequencies(2);
    temporal_offsets_secs('W-P') = world_pupil_temporal_offset_secs;

    figure ;
    t = tiledlayout(2,1);
    h = sgtitle(sprintf("Sensor Fits | M: %d", nn));
    h.FontWeight = 'bold';

    % First tile we will plot the slow video world camera and its fit
    nexttile;
    title(sprintf("World and MS-AS | F: %.3f", frequencies(1)));
    hold on;
    yyaxis left
    plot(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), '.k','MarkerSize',10, 'DisplayName', 'World Measured');
    plot(slow_measurement.W.t, world_fit_slow_TS, '-k', 'DisplayName', 'World Fit');
    yyaxis right
    plot(slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.r','MarkerSize',20, 'DisplayName', 'AS Measured');
    plot(slow_measurement.M.t, AS_fit, '-r', 'DisplayName', 'AS Fit');
    legend show ;

    % Next tile we will print the slow video TS chip and its fit
    % nexttile ;
    % title(sprintf("MS-TS Fit | F: %.3f", frequencies(1)));
    % hold on;
    % plot(slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.TS(:, 1)), '.', 'DisplayName', 'Measured');
    % plot(slow_measurement.M.t, TS_fit, '.', 'DisplayName', 'Fit');
    % 
    % legend show ;

    % Next tile we will print the Fast video world camera and its fit
    nexttile;
    title(sprintf("World and pupil | F: %.3f", frequencies(2)));
    hold on;
    yyaxis left
    plot(fast_measurement.W.t, convert_to_contrast(fast_measurement.W.v), '.k','MarkerSize',10, 'DisplayName', 'World Measured');
    plot(fast_measurement.W.t, world_fit_fast, '-k', 'DisplayName', 'World Fit');
    yyaxis right
    plot(fast_measurement.P.t, convert_to_contrast(fast_measurement.P.v), '.b','MarkerSize',10, 'DisplayName', 'Pupil Measured');
    plot(fast_measurement.P.t(nPupilFramesToDiscard:end), pupil_fit, '-b', 'DisplayName', 'Pupil Fit');
    xlim([2 3]);
    legend show ;


    % Plot the sensors before and after alignment
    figure ;
    t = tiledlayout(2, 2);
    h = sgtitle(sprintf("Sensor Phase Alignment | M: %d", nn));
    h.FontWeight = 'bold';

    % First, plot the sensors before the alignment
    nexttile;
    title(sprintf("%.3f | Before Alignment", frequencies(1)));
    hold on ;
    yyaxis left
    plot(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), '.', 'DisplayName', 'World');
    yyaxis right
    plot(slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.', 'DisplayName', 'MS-AS');
    plot(slow_measurement.M.t, convert_to_contrast(slow_measurement.M.v.TS(:, 1)), '.', 'DisplayName', 'MS-TS');
    xlabel("Time [s]");
    ylabel("Contrast");

    % Show the legend for this plot
    legend show;

    % First, plot the sensors after the alignment
    nexttile;
    title(sprintf("%.3f | After Alignment", frequencies(1)));
    hold on ;
    yyaxis left
    plot(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), '.', 'DisplayName', 'World');
    yyaxis right
    plot(slow_measurement.M.t+world_ms_AS_temporal_offset_secs(nn), convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.', 'DisplayName', 'MS-AS');
    plot(slow_measurement.M.t+world_ms_TS_temporal_offset_secs(nn), convert_to_contrast(slow_measurement.M.v.TS(:, 1)), '.', 'DisplayName', 'MS-TS');
    xlabel("Time [s]");
    ylabel("Contrast");

    % Show the legend for this plot
    legend show;

    % First, plot the sensors before the alignment
    nexttile;
    title(sprintf("%.3f | Before Alignment", frequencies(2)));
    hold on ;
    yyaxis left
    plot(fast_measurement.W.t, convert_to_contrast(fast_measurement.W.v), '.', 'DisplayName', 'World');
    yyaxis right
    plot(fast_measurement.P.t(nPupilFramesToDiscard:end), convert_to_contrast(fast_measurement.P.v(nPupilFramesToDiscard:end)), '.', 'DisplayName', 'Pupil');
    xlabel("Time [s]");
    ylabel("Contrast");

    % Show the legend for this plot
    legend show;

    % First, plot the sensors after the alignment
    nexttile;
    title(sprintf("%.3f | After Alignment", frequencies(2)));
    hold on ;
    yyaxis left
    plot(fast_measurement.W.t, convert_to_contrast(fast_measurement.W.v), '.', 'DisplayName', 'World');
    yyaxis right
    plot(fast_measurement.P.t(nPupilFramesToDiscard:end)+world_pupil_temporal_offset_secs(nn), convert_to_contrast(fast_measurement.P.v(nPupilFramesToDiscard:end)), '.', 'DisplayName', 'Pupil');
    xlabel("Time [s]");
    ylabel("Contrast");

    % Show the legend for this plot
    legend show;

end


end


%% LOCAL FUNCTIONS

% Calculate the phase offset between the sensors and plot them before and after adjustment
function [phase_offset, A_fit, B_fit] = calculate_phase_offset(sensorA_t, sensorA_v, sensorB_t, sensorB_v, frequency)
% First, fit the two waves independently
[A_r2, A_amplitude, A_phase, A_fit] = fourierRegression( sensorA_v, sensorA_t, frequency );
[B_r2, B_amplitude, B_phase, B_fit] = fourierRegression( sensorB_v, sensorB_t, frequency );

% Calculate the phase offset
phase_offset = A_phase - B_phase;

% Wrap to the pi/2 domain
if phase_offset > pi/2
    phase_offset = phase_offset - pi;
end

return ;

end

% Convert raw data in arbirtary units to contrast units
function contrast_v = convert_to_contrast(v)
% Take the mean of the input vector
m = mean(v);

% Convert to contrast
contrast_v = (v - m) / m;
return ;

end