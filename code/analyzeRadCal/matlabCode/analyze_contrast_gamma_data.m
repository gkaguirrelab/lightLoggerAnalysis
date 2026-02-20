function analyze_contrast_gamma_data(calibration_metadata, measurements)
% Analyze the results of a contrast gamma light logger calibration measurement (post-conversion)
%
% Syntax:
%  analyze_contrast_gamma_data(calibration_metadata, measurements)
%
% Description:
%   Given the parsed and converted metadata for a contrast gamma 
%   calibration measurement, analyze the data and plot the contrast 
%   gamma plot. 
%   
% Inputs:
%   calibration_metadata        - Struct. Converted metadata for 
%                                 the contrast gamma reading 
%   
%   measurements                - Cell. The parsed + converted 
%                                 contrast gamma readings
%                              
%
% Examples:
%{
    path_to_experiment = "/example/path"; 
    converted_light_logger_data = convert_light_logger_calibration_data(path_to_experiment, true, true true, true); 
    analyze_contrast_gamma_data(converted_light_logger_data.metdata.contrast_gamma, converted_light_logger_data.readings.contrast_gamma);
%}  
    arguments 
        calibration_metadata; % The parsed + converted metadata for the contrast gamma measuremnet 
        measurements; % The parsed + conveerted contrast gamma measurements 
    end
    
    % Retrieve the NDFs, frequencies and contrast levels used to make the measurements
    % as well as the number of measurements made at each level
    NDFs = calibration_metadata.NDFs; 
    contrast_levels = calibration_metadata.contrast_levels;
    frequencies = calibration_metadata.frequencies;
    n_measures = calibration_metadata.n_measures;


    ND1Orange = [0.8500, 0.3250, 0.0980];  % Orange we use for ND1 in other cal plots
    %if we ever rerun this with multiple NDFs, we should change the colors.

    % First, let's iterate over NDF level
    fit_slopes_per_NDF = zeros(numel(NDFs), 1);
    for nn = 1:numel(NDFs)
        % Retrieve the current NDF 
        NDF = NDFs(nn); 

        % Next, iterate over the frequencies 
        for ff = 1:numel(frequencies)
            % Retrieve the current frequency 
            frequency = frequencies(ff); 

            % Initialize a vector to save the amplitudes per contrast level 
            amplitudes_by_contrasts = zeros(numel(contrast_levels), n_measures); 

            % Initialize a vector to save measurements an their fits 
            measurements_and_fits = cell(numel(contrast_levels), n_measures); 

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

                    % Save the measurement and the fit 
                    measurements_and_fits{cc, mm} = {{world_t, world_v_contrast}, {world_t, world_fit}}; 

                end % Measure 

            end % Contrast 

            % Initialize a tiled layout plot to display the results 
            [rows, cols] = find_min_figsize(numel(frequencies)); 
            NDF_fig = figure; 
            NDF_plot = tiledlayout(rows, cols, 'Parent', NDF_fig); 


            % Calculate the average amplitudes for each contrast level at this frequency 
            mean_amplitudes_by_contrast = mean(amplitudes_by_contrasts, 2); 

            % Retrieve the NDF plot axis 
            nexttile(NDF_plot); 

            % Plot the mean amplitudes by contrast 
            plot(contrast_levels, sort(mean_amplitudes_by_contrast), 'o',...
                'MarkerEdgeColor', ND1Orange, 'MarkerFaceColor',ND1Orange, 'DisplayName', 'Data'); 
            hold on; 
            title(sprintf("Contrast Gamma | NDF: %.2f | F: %.2f", NDF, frequency)); 
            xlabel("Contrast"); 
            ylabel("Mean Amplitude"); 
            xlim([-0.05 1.05]);
            ylim([-0.05 1.05]);
            axis square

            % Add a linear fit
            p = polyfit(contrast_levels, mean_amplitudes_by_contrast, 1);
            fit_slopes_per_NDF(nn) = p(1); 
            plot(contrast_levels, polyval(p, contrast_levels), '-r', 'DisplayName', 'Fit');
            plot([0 1], [0 1], ':k', 'DisplayName', 'Reference Line');

            % Display the legend 
            legend show; 

        end % Frequency 

    end % NDF 


    % Plot the slopes of the fits per NDF 
    figure; 
    title("Fit Slope By NDF Level"); 
    hold on; 


    plot(NDFs, fit_slopes_per_NDF, '-o', ...
    'DisplayName', 'Slopes', ...
    'MarkerSize', 10, ...
    'LineWidth', 2.5, ...
    'Color', [0 1 1], ...          % cyan line
    'MarkerFaceColor', [0 1 1], ... % filled cyan marker
    'MarkerEdgeColor', [1 1 1]);    % white edge for contrast
    xlabel("NDF"); 
    xticks(NDFs);
    ylabel("Slope"); 
    legend show; 


end