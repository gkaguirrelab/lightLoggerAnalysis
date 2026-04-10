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


    colorList = [
    0.6350, 0.0780, 0.1840;   % Red
    0.8500, 0.3250, 0.0980;   % Orange
    0.9290, 0.6940, 0.1250;   % Yellow
    0.4660, 0.6740, 0.1880;   % Green
    0.3010, 0.7450, 0.9330;   % Light Blue
    0,      0.4470, 0.7410;   % Blue
    0.4940, 0.1840, 0.5560;   % Purple
    ];

   fit_slopes_per_NDF = zeros(numel(NDFs), numel(frequencies));

    for ff = 1:numel(frequencies)
        frequency = frequencies(ff);

        figure;
        hold on;
        title(sprintf("Contrast Gamma | F: %.2f", frequency));
        xlabel("Contrast");
        ylabel("Mean Amplitude");
        xlim([-0.05 1.05]);
        ylim([-0.05 1.05]);
        axis square

        for nn = 1:numel(NDFs)
            NDF = NDFs(nn);

            amplitudes_by_contrasts = zeros(numel(contrast_levels), n_measures);
            measurements_and_fits = cell(numel(contrast_levels), n_measures);

            for cc = 1:numel(contrast_levels)
                for mm = 1:n_measures
                    measurement = measurements{nn, cc, ff, mm};

                    world_t = measurement.W.t;
                    world_v = measurement.W.v;

                    world_v_contrast = (world_v - mean(world_v)) / mean(world_v);
                    [~, world_amplitude, ~, world_fit] = fourierRegression(world_v_contrast, world_t, frequency);

                    amplitudes_by_contrasts(cc, mm) = world_amplitude;
                    measurements_and_fits{cc, mm} = {{world_t, world_v_contrast}, {world_t, world_fit}};
                end
            end

            mean_amplitudes_by_contrast = mean(amplitudes_by_contrasts, 2);

            plot(contrast_levels, mean_amplitudes_by_contrast, '-o', ...
                'Color', colorList(nn,:), ...
                'MarkerFaceColor', colorList(nn,:), ...
                'MarkerEdgeColor', colorList(nn,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('NDF%.0f', NDF));

            p = polyfit(contrast_levels, mean_amplitudes_by_contrast, 1);
            fit_slopes_per_NDF(nn, ff) = p(1);
        end

        plot([0 1], [0 1], ':k', 'DisplayName', 'Reference Line');
        legend('Location', 'best');
        hold off;
    end

end