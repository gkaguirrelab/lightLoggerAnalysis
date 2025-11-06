function analyze_phase_fit_data(calibration_metadata, measurements)
% Analyze the results of a phase fitting light logger calibration measurement (post-conversion)
%
% Syntax:
%  analyze_phase_fit_data(calibration_metadata, measurements)
%
% Description:
%   Given the parsed and converted metadata for a phase fitting  
%   calibration measurement, analyze the data and calculate the phase 
%   offset in seconds between pairs of sensors. 
%   
% Inputs:
%   calibration_metadata        - Struct. Converted metadata for 
%                                 the phase fitting reading 
%   
%   measurements                - Cell. The parsed + converted 
%                                 phase fitting readings
%                              
%
% Examples:
%{
    path_to_experiment = "/example/path"; 
    converted_light_logger_data = convert_light_logger_calibration_data(path_to_experiment, true, true true, true); 
    analyze_ms_linearity_data(converted_light_logger_data.metdata.phase_fitting, converted_light_logger_data.readings.phase_fitting);
%}
    arguments 
        calibration_metadata; % The parsed + converted metadata for the phase alignment measuremnet 
        measurements; % The parsed + converted metadata for the phase alignment measuremnet 
    end 

    % Specify the number of measurements to discard from the start of the pupil
    % recording (as there is a start-up effect on the signal here)
    nPupilFramesToDiscard = 120;

    % Retrieve the frequencies and contrast levels used to make the measurements
    % as well as the number of measurements made at each level (this should be 1 here)
    NDFs = calibration_metadata.NDFs; 
    contrast_levels = calibration_metadata.contrast_levels;
    frequencies = calibration_metadata.frequencies;
    n_measures = calibration_metadata.n_measures;

    % First, we will iterate over the NDFs
    for nn = 1:numel(NDFs)
        % Retrieve the current NDF
        NDF = NDFs(nn); 

        % Iterate over the contrast levels 
        for cc = 1:numel(contrast_levels)
            % Retrieve the current contrast level
            contrast = contrast_levels(cc); 

            % Iterate over the frequencies 
            for ff = 1:numel(frequencies)
                % Retrieve the current frequency 
                frequency = frequencies(ff); 

                % Initialize a variable to save the phase offsets
                % for our two sensors over the course of the n measures
                temporal_offset_seconds_by_frequency = zeros(2, n_measures); 
                ms_temporal_offset_mat_row = 1; 
                pupil_temporal_offset_mat_row = 2; 

                % Iterate over the measurments at this frequency 
                for mm = 1:n_measures
                    % Retrieve the measurement 
                    measurement = measurements{nn, cc, ff, mm}; 

                    %%%%%%{ Extract World Measurement %}%%%%%%
                    world_measurement = measurement.W; 
                    world_t = world_measurement.t; 
                    world_v = world_measurement.v; 
                    world_v_contrast = (world_v - mean(world_v)) / mean(world_v); 
                    
                    %%%%%%{ Extract Pupil Measurement %}%%%%%%
                    pupil_measurement = measurement.P; 
                    pupil_t = pupil_measurement.t; 
                    pupil_v = pupil_measurement.v; 
                    pupil_v_contrast = (pupil_v - mean(pupil_v)) / mean(pupil_v); 

                    %%%%%%{ Extract MS_AS Measurement %}%%%%%%
                    MS_AS_measurement = measurement.M; 
                    MS_AS_t = MS_AS_measurement.t.AS; 
                    MS_AS_v = MS_AS_measurement.v.AS(:, 5); 
                    MS_AS_v_contrast = (MS_AS_v - mean(MS_AS_v)) / mean(MS_AS_v); 

                    % Now, we will calculate the phase offset between 
                    % the sensors at this measure 

                    %%%%%%{ Calculate World - MS offset %}%%%%%%
                    [world_AS_phase_diff, world_fit_AS, AS_fit] = calculate_phase_offset(world_t, world_v_contrast, ...
                                                                                         MS_AS_t, MS_AS_v_contrast, ...
                                                                                         frequency...
                                                                                       );

                    %%%%%%{ Calculate World - Pupil offset %}%%%%%%
                    [world_pupil_phase_diff, world_fit_fast, pupil_fit] = calculate_phase_offset(world_t, world_v_contrast, ...
                                                                                                 pupil_t(nPupilFramesToDiscard:end), pupil_v_contrast(nPupilFramesToDiscard:end), ...
                                                                                                 frequency...
                                                                                                );

                    % Calculate the temporal offsets in seconds
                    world_AS_temporal_offset_seconds = (world_AS_phase_diff / (2 * pi)) / frequency;
                    world_pupil_temporal_offset_seconds = (world_pupil_phase_diff / (2 * pi)) / frequency;

                    % Let's save these offsets for this measure 
                    temporal_offset_seconds_by_frequency(ms_temporal_offset_mat_row, mm) = world_AS_temporal_offset_seconds; 
                    temporal_offset_seconds_by_frequency(pupil_temporal_offset_mat_row, mm) = world_pupil_temporal_offset_seconds; 


                    % Let's plot the results for this measure 
                    % TODO: 
                    %}
                                                                                                    

                end % Measure

                % Calculate the mean world - AS temporal offset 
                % and mean world - pupil temporal offset
                mean_temporal_offset_seconds_by_frequency = mean(temporal_offset_seconds_by_frequency, 2); 
                std_temporal_offset_seconds_by_frequency = std(temporal_offset_seconds_by_frequency, [], 2); 

                % Report these values to the terminal 
                fprintf("NDF: %.2f | C: %.2f | F: %.2f | Mean Offset (W-AS): %.3f [ms] | STD Offset: %.3f\n [ms]", NDF, contrast, frequency,...
                                                                                                              mean_temporal_offset_seconds_by_frequency(ms_temporal_offset_mat_row) * 1000,...
                                                                                                              (std_temporal_offset_seconds_by_frequency(ms_temporal_offset_mat_row)* 1000)...
                       );
                fprintf("NDF: %.2f | C: %.2f | F: %.2f | Mean Offset (W-P): %.3f [ms] | STD Offset: %3.f\n [ms]", NDF, contrast, frequency,...
                                                                                                             mean_temporal_offset_seconds_by_frequency(pupil_temporal_offset_mat_row) * 1000,...
                                                                                                             (std_temporal_offset_seconds_by_frequency(pupil_temporal_offset_mat_row)* 1000)...
                       );

            end % Frequency

        end % Contrast 

    end % NDF 



    % Now, we will iterate over the measurements. We will do this to take the average phase
    temporal_offsets_secs = containers.Map({'W-AS', 'W-TS', 'W-P'}, { zeros([n_measures, 1]), zeros([n_measures, 1]), zeros([n_measures, 1]) });

    for nn = 1:n_measures
        % Retrieve the sensor readings for this measurement
        slow_measurement = measurements{1, 1, 1, nn};
        fast_measurement = measurements{1, 1, 2, nn};

        % Adjust the start times to be relative to world camera
        startTime = slow_measurement.W.t(1);
        slow_measurement.W.t = slow_measurement.W.t - startTime;
        slow_measurement.P.t = slow_measurement.P.t - startTime;
        slow_measurement.M.t.AS = slow_measurement.M.t.AS - startTime;
        slow_measurement.S.t = slow_measurement.S.t - startTime;

            startTime = fast_measurement.W.t(1);
        fast_measurement.W.t = fast_measurement.W.t - startTime;
        fast_measurement.P.t = fast_measurement.P.t - startTime;
        fast_measurement.M.t.AS = fast_measurement.M.t.AS - startTime;
        fast_measurement.S.t = fast_measurement.S.t - startTime;
        
        % Calculate the phase difference of the MS chips to the world camera
        % from the slow measurement. Calculate and store the temporal offset
        % between sensors in units of seconds.
        [world_ms_AS_phase_diff, world_fit_slow_AS, AS_fit] = calculate_phase_offset(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), ...
            slow_measurement.M.t.AS, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), ...
            frequencies(1));

        world_ms_AS_temporal_offset_secs = temporal_offsets_secs('W-AS');
        world_ms_AS_temporal_offset_secs(nn) = (world_ms_AS_phase_diff / (2 * pi)) / frequencies(1);
        temporal_offsets_secs('W-AS') = world_ms_AS_temporal_offset_secs;

        [world_ms_TS_phase_diff, world_fit_slow_TS, TS_fit] = calculate_phase_offset(slow_measurement.W.t, convert_to_contrast(slow_measurement.W.v), ...
            slow_measurement.M.t.AS, convert_to_contrast(slow_measurement.M.v.TS(:, 1)), ...
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
        plot(slow_measurement.W.t, world_fit_slow_AS, '-k', 'DisplayName', 'World Fit');
        yyaxis right
        plot(slow_measurement.M.t.AS, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.r','MarkerSize',20, 'DisplayName', 'AS Measured');
        plot(slow_measurement.M.t.AS, AS_fit, '-r', 'DisplayName', 'AS Fit');
        legend show ;

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
        plot(slow_measurement.M.t.AS, convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.', 'DisplayName', 'MS-AS');
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
        plot(slow_measurement.M.t.AS+world_ms_AS_temporal_offset_secs(nn), convert_to_contrast(slow_measurement.M.v.AS(:, 5)), '.', 'DisplayName', 'MS-AS');
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