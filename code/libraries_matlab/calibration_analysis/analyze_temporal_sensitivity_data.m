function analyze_temporal_sensitivity_data(calibration_metadata, measurements)
% Analyze the results of a temporal sensitivity light logger calibration measurement (post-conversion)
%
% Syntax:
%  analyze_temporal_sensitivity_data(calibration_metadata, measurements)
%
% Description:
%   Given the parsed and converted metadata for a temporal sensitivity 
%   calibration measurement, analyze the data and plot the TTF across 
%   light levels. 
%   
% Inputs:
%   calibration_metadata        - Struct. Converted metadata for 
%                                 the temporal sensitivity reading 
%   
%   measurements                - Cell. The parsed + converted 
%                                 temporal sensitivity readings
%                              
%
% Examples:
%{
    path_to_experiment = "/example/path"; 
    converted_light_logger_data = convert_light_logger_calibration_data(path_to_experiment, true, true true, true); 
    analyze_ms_linearity_data(converted_light_logger_data.metdata.temporal_sensitivity, converted_light_logger_data.readings.temporal_sensitivity);
%}
    arguments 
        calibration_metadata; % The parsed + converted metadata for the temporal sensitivity measuremnet 
        measurements; % The parsed + converted metadata for the temporal sensitivity measuremnet 
    end 
    
    % Define the world FPS (Note: if you changed FPS elsewhere, you will need to change it here too)
    world_fps = 120; 

    % Retrieve some relevant metadata about the experiment
    NDFs = calibration_metadata.NDFs; 
    contrast_levels = calibration_metadata.contrast_levels;
    frequencies = calibration_metadata.frequencies;
    n_measures = calibration_metadata.n_measures;

    % Retrieve the contrast attentuation of the CombiLED at higher level frequencies 
    % to combat rolloff. 
    contrast_attenuation_with_frequency = contrastAttenuationByFreq(frequencies); 

    % First, we will iterate over the NDFs 
    for cc = 1:numel(contrast_levels)
        % Retrieve the current contrast level 
        contrast_level = contrast_levels(cc); 

        % Initialize a matrix to save amplitude of response for this contrast level 
        response_amplitude_data = zeros(numel(NDFs), numel(frequencies), n_measures); 

        % Next, we will iterate over the contrast levels 
        % at this NDF level 
        for nn = 1:numel(NDFs)
            % Retrieve the NDF we used for this recording 
            NDF = NDFs(nn); 

            % Iterate over the measurements at each frequency 
            for mm = 1:n_measures
                % Initialize figures for the measured vs fit and AGC settings 
                % for this measurement 
                measured_fit_fig = figure; 
                measured_fit_plot = tiledlayout(measured_fit_fig, numel(frequencies), 1);
                sgtitle(measured_fit_plot, sprintf("Measured vs Fit | C: %.2f | NDF: %.2f | M: %d", contrast_level, NDF, mm)); 

                agc_settings_fig = figure; 
                agc_settings_plot = tiledlayout(agc_settings_fig, numel(frequencies), 1); 
                sgtitle(agc_settings_plot, sprintf("AGC Performance | C: %.2f | NDF: %.2f | M: %d", contrast_level, NDF, mm)); 

                % Iterate over the frequences 
                for ff = 1:numel(frequencies)
                    % Retrieve the current frequency 
                    frequency = frequencies(ff);

                    % Move to the next tile for each of the above plots 
                    measured_fit_ax = nexttile(measured_fit_plot);
                    agc_settings_ax = nexttile(agc_settings_plot); 

                    % Retrieve the given measurement
                    measurement = measurements{nn, cc, ff, mm}; 

                    % Retrieve the relevant world variables
                    world_t = measurement.W.t; % t is time units converted to seconds
                    world_v = measurement.W.v; % v here is the mean pixel of each frame
                    world_settings = measurement.W.settings; 

                    % Now, we will perform the fitting procedure of the observed measurements to
                    % a fit sinusoid. Though first, convert the world to contrast
                    world_v_contrast = (world_v - mean(world_v)) / mean(world_v);
                    [world_r2, world_amplitude, ~, world_fit] = fourierRegression( world_v_contrast, world_t, frequency );

                    % Now, let's save the amplitude 
                    response_amplitude_data(nn, ff, mm) = world_amplitude; 

                    %%%%%%%%%%%%%%%%{ Plot measured vs Fit %}%%%%%%%%%%%%%%%%%
                    local_world_t = world_t - world_t(1);

                    plot(measured_fit_ax, local_world_t, world_v_contrast, 'k.', 'DisplayName', sprintf('WorldMeasured Freq=%2.1f Hz',frequency));
                    hold(measured_fit_ax, 'on'); 
                    plot(measured_fit_ax, local_world_t, world_fit, 'r-', 'DisplayName', sprintf('WorldFit R2=%.3f', world_r2));
                    xlabel(measured_fit_ax, "Time [s]");
                    ylabel(measured_fit_ax, "Contrast");
                    ylim(measured_fit_ax, [-1, 1]); % Set to 8 bit unsigned range

                    % Show the legend for this tile
                    % legend show;
                    title(measured_fit_ax, sprintf("F: %2.2f | R2: %2.2f", frequency, world_r2));
                    legend(measured_fit_ax);


                    %%%%%%%%%%%%%%%%{ Plot AGC performance %}%%%%%%%%%%%%%%%%%
                    axes(agc_settings_ax); 
                    yyaxis left; 
                    plot(world_t, world_settings.Again, '-b', 'DisplayName', 'AnalogueGain'); 
                    hold on; 
                    plot(world_t, world_settings.Dgain, '-r', 'DisplayName', 'DigitalGain');
                    ylabel("Gain Value"); 
        
                    yyaxis right; 
                    plot(world_t, world_settings.exposure, '-g', 'DisplayName', 'Exposure [s]'); 
        
                    % Label the plot 
                    title(sprintf("F: %2.2f | R2: %2.2f", frequency, world_r2));
                    xlabel("Time [s]"); 
                    ylabel("Exposure Time [s]"); 
                    legend show; 

                end % Frequency

            end % Measure 

        end % NDF  

        % Contrast the TTF per measure per contrast and NDF 
        plot_TTF(NDFs, contrast_level, frequencies, response_amplitude_data, world_fps, contrast_attenuation_with_frequency); 

        % Construct the TTF of the average of the measurements
        plot_mean_TTF(NDFs, contrast_level, frequencies, response_amplitude_data, world_fps, contrast_attenuation_with_frequency); 


    end % Contrast

end

%% LOCAL FUNCTIONS
function [filterFreqHz,filterAmp] = approxFreqFilter(frequencies, amplitudes)

    % Trim off the 100 Hz case; this is dominated by signed noise in the
    % amplitude estimation
    frequencies = frequencies(frequencies~=100);
    amplitudes = amplitudes(frequencies~=100);

    myObj = @(p) norm(amplitudes - p(1)*idealDiscreteSampleFilter(frequencies,1/p(2)));
    x0 = [0.8,100];
    lb = [0, 75];
    ub = [1, 200];
    f = fmincon(myObj,x0,[],[],[],[],lb,ub);
    filterFreqHz = f(2);
    filterAmp = f(1);


end

% Local function to plot the TTF of all measures of a certain contrast level 
function plot_TTF(NDFs, contrast_level, frequencies, response_amplitude_data, world_fps, contrast_attenuation_with_frequency)
        % Construct the Mean TTF for this contrast level 
        figure ;    

        % Plot each NDF line 
        for nn = 1:numel(NDFs)
            % Retrieve the mean amplitude for each frequency for this NDF 
            response_amplitude_per_frequency = squeeze(response_amplitude_data(nn, :, :));

            % Plot the different measures of each frequency 
            for mm = 1:size(response_amplitude_per_frequency, 2)
                % Retrieve the response amplitude for the given measure 
                measure_response_amplitude = response_amplitude_per_frequency(:, mm)'; 

                % Adjust for the contrast of the stimulus and for the roll-off in
                % modulation depth with temporal frequency. The splicing here is just for debugging if you 
                % ran with just the first few frequencies 
                measure_amp_data = measure_response_amplitude .* (1./contrast_attenuation_with_frequency(1:numel(frequencies))) ./ contrast_level;

                disp(size(measure_amp_data))

                plot(log10(frequencies), measure_amp_data,...
                    '-x',...
                    'MarkerSize',15,...
                    'LineWidth',2,...
                    'DisplayName', sprintf("NDF %.2f | M: %d", NDFs(nn), mm)...
                    );
                hold on; 

            end % Measure 

        end % NDF 

        % Generate the ideal device curve for the world camera
        xfine = logspace(log10(frequencies(1)),log10(frequencies(end)),100);
        ideal_device = idealDiscreteSampleFilter(xfine, 1/world_fps);
        plot(log10(xfine), ideal_device, "--", "Color",[0.5 0.5 0.5], "DisplayName", "Ideal Device");

        % Add a title
        title(sprintf("World TTF All Measures | C: %.2f", contrast_level));

        % Label and clean up the plot
        xlabel("Frequency [hz]");
        ylabel("Relative Contrast Response");
        legend show;
        a = gca();
        a.XTick = log10(frequencies);
        a.XTickLabels = arrayfun(@(x) num2str(x),frequencies,'UniformOutput',false);
end 


% Local function to plot the mean TTF 
function plot_mean_TTF(NDFs, contrast_level, frequencies, response_amplitude_data, world_fps, contrast_attenuation_with_frequency)
        % Construct the Mean TTF for this contrast level 
        figure ;   

        % We will now take the mean of all measures per frequency 
        mean_response_amplitude_data = mean(response_amplitude_data, 3); 

        % Plot each NDF line 
        for nn = 1:numel(NDFs)
            % Retrieve the mean amplitude for each frequency for this NDF 
            mean_response_amplitude_per_frequency = mean_response_amplitude_data(nn, :, :);

            % Obtain and plot the mean amplitude of response for each frequency across measures.
            % Adjust for the contrast of the stimulus and for the roll-off in
            % modulation depth with temporal frequency. The splicing here is just for debugging if you 
            % ran with just the first few frequencies 
            meanAmpData = mean_response_amplitude_per_frequency .* (1./contrast_attenuation_with_frequency(1:numel(frequencies))) ./ contrast_level;
            plot(log10(frequencies), meanAmpData,...
                '-x',...
                'MarkerSize',15,...
                'LineWidth',2,...
                'DisplayName', sprintf("NDF %.2f", NDFs(nn))...
                );
            hold on; 

        end 

        % Generate the ideal device curve for the world camera
        xfine = logspace(log10(frequencies(1)),log10(frequencies(end)),100);
        ideal_device = idealDiscreteSampleFilter(xfine, 1/world_fps);
        plot(log10(xfine), ideal_device, "--", "Color",[0.5 0.5 0.5], "DisplayName", "Ideal Device");

        % Get the approximate filter frequency and amplitude for the data
        [filterFreqHz,filterAmp] = approxFreqFilter(frequencies, meanAmpData);
        fitVals = filterAmp*idealDiscreteSampleFilter(xfine, 1/filterFreqHz);
        plot(log10(xfine), fitVals, "-", "Color",[1 0 0], "DisplayName", "Data fit");

        % Add a title
        title(sprintf("World TTF Avg Amplitude | C: %.2f | Approx freq %2.2f Hz", contrast_level, filterFreqHz));

        % Label and clean up the plot
        xlabel("Frequency [hz]");
        ylabel("Mean Relative Contrast Response");
        legend show;
        a = gca();
        a.XTick = log10(frequencies);
        a.XTickLabels = arrayfun(@(x) num2str(x),frequencies,'UniformOutput',false);


end 