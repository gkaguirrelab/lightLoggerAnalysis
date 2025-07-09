function ttfFigHandle = analyze_temporal_sensitivity_data(calibration_metadata, measurements)

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


        % Construct the TTF for this contrast level 
        figure ;   

        % We will now take the mean of all measures per frequency 
        mean_response_amplitude_data = mean(response_amplitude_data, 3); 

        % Plot each NDF line 
        for nn = 1:numel(NDFs)
            % Retrieve the mean amplitude for each frequency for this NDF 
            mean_response_amplitude_per_frequency = mean_response_amplitude_data(nn, :, :);

            % Retrieve the amplitudes that we manually measured with the klein for certain frequencies
            test_mod_depth_amplitudes = getpref("lightLoggerAnalysis", "test_mod_depth_amplitudes");

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
        ideal_device = idealDiscreteSampleFilter(xfine, 1/getpref("lightLoggerAnalysis", "world_fps"));
        plot(log10(xfine), ideal_device, "--", "Color",[0.5 0.5 0.5], "DisplayName", "Ideal Device");

        % Get the approximate filter frequency and amplitude for the data
        [filterFreqHz,filterAmp] = approxFreqFilter(frequencies, meanAmpData);
        fitVals = filterAmp*idealDiscreteSampleFilter(xfine, 1/filterFreqHz);
        plot(log10(xfine), fitVals, "-", "Color",[1 0 0], "DisplayName", "Data fit");

        % Add a title
        title(sprintf("World TTF | C: %.2f | Approx freq %2.2f Hz", contrast_level, filterFreqHz));

        % Label and clean up the plot
        xlabel("Frequency [hz]");
        ylabel("Relative Contrast Response");
        legend show;
        a = gca();
        a.XTick = log10(frequencies);
        a.XTickLabels = arrayfun(@(x) num2str(x),frequencies,'UniformOutput',false);

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