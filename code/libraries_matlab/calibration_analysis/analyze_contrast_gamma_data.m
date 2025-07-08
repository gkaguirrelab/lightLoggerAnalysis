function analyze_contrast_gamma_data(calibration_metadata, measurements)
    % Retrieve the NDFs, frequencies and contrast levels used to make the measurements
    % as well as the number of measurements made at each level
    NDFs = calibration_metadata.NDFs; 
    contrast_levels = calibration_metadata.contrast_levels;
    frequencies = calibration_metadata.frequencies;
    n_measures = calibration_metadata.n_measures;


    % First, let's iterate over NDF level 
    for nn = 1:numel(NDFs)
        % Retrieve the current NDF 
        NDF = NDFs(nn); 

        % Initialize a tiled layout plot to display the results 
        [rows, cols] = find_min_figsize(numel(frequencies)); 
        figure; 
        NDF_plot = tiledlayout(rows, cols); 

        % Next, iterate over the frequencies 
        for ff = 1:numel(frequencies)
            % Move to the next tile of the plot 
            frequency_contrast_gamma_plot = nexttile(NDF_plot); 

            % Retrieve the current frequency 
            frequency = frequencies(ff); 

            % Initialize a vector to save the amplitudes per contrast level 
            amplitudes_by_contrasts = zeros(numel(contrast_levels), n_measures); 

            % Next, iterate over the contrast levels 
            % at this frequency 
            for cc = 1:numel(contrast_levels)
                % Iterate over the measures over this frequnecy + contrast + NDF 
                for mm = 1:n_measures
                    % Retrieve the measurement
                    measurement = measurements{nn, cc, ff, mm};

                    % Extract the camera temporal support and v vectors
                    world_t = measurement.W.t;
                    world_v = measurement.W.v;

                    % Fit the modulation to calculate the amplitude
                    % First, convert world to contrat
                    world_v_contrast = (world_v - mean(world_v)) / mean(world_v);
                    [world_r2, world_amplitude, ~, world_fit] = fourierRegression( world_v_contrast, world_t, frequency );

                    % Save the amplitude for this measurement of this contrast + frequency + NDF 
                    amplitudes_by_contrasts(cc, mm) = world_amplitude; 

                end % Measure 

            end % Contrast

            % Calculate the average amplitudes for each contrast level at this frequency 
            mean_amplitudes_by_contrast = mean(amplitudes_by_contrasts, 2); 

            % Plot the mean amplitudes by contrast 
            plot(contrast_levels, mean_amplitudes_by_contrast, '-x', 'DisplayName', 'Data'); 
            hold on; 
            title(sprintf("Contrast Gamma | NDF: %.2f | F: %.2f", NDF, frequency)); 
            xlabel("Contrast"); 
            ylabel("Mean Amplitude"); 
            xlim([-0.05 1.05]);
            ylim([-0.05 1.05]);
            axis square

            % Add a linear fit
            p = polyfit(contrast_levels, mean_amplitudes_by_contrast, 1);
            plot(contrast_levels, polyval(p, contrast_levels), '-r', 'DisplayName', 'Fit');
            plot([0 1], [0 1], ':k', 'DisplayName', 'Reference Line');

            % Display the legend 
            legend show; 

        end % Frequency 

    end % NDF 

end