function analyze_contrast_gamma_data(CalibrationData, sorted_measurements, NDF, Pi_util)
% Retrieve the frequencies and contrast levels used to make the measurements
% as well as the number of measurements made at each level (this should be 1 here)
contrast_levels = CalibrationData.contrast_levels;
frequencies = CalibrationData.frequencies;
n_measures = CalibrationData.n_measures;

assert(numel(frequencies) == 1); % Assert this was all made at a single frequency

% Sorted data is based into here as a 3D cell array in the shape
% (contrastIdx, frequencyIdx, measurementIdx)

% First, retrieve the single frequency exposed under all conditions
frequency = frequencies(1);

% Initialize an array of avg amplitudes that we will then
% fill in for each contrast level
avg_amplitudes = zeros(size(contrast_levels));

% Then, let's iterate over the contrast levels
for cc = 1:numel(contrast_levels)
    % Retrieve the contrast at this level
    contrast = contrast_levels(cc);

    % Then, iterate over the measures at this contrast level
    local_amplitudes = zeros([n_measures, 1]);
    for nn = 1:n_measures
        % Retrieve the measurement
        measurement = sorted_measurements{cc, 1, nn};

        % Extract the camera temporal support and v vectors
        world_t = measurement.W.t;
        world_v = measurement.W.v;

        % Fit the modulation to calculate the amplitude
        % First, convert world to contrat
        world_v_contrast = (world_v - mean(world_v)) / mean(world_v);
        [world_r2, world_amplitude, ~, world_fit] = fourierRegression( world_v_contrast, world_t, frequency );

        % Check that the fits are good
        if contrast > 0 && world_r2<0.95
            warning('poor R2 fit')
        end

        % Illustrate this fit
        %{
            figure ; 
            title(sprintf("World Contrast Gamma Fit C: %.3f | N: %d", contrast, nn)); 
            hold on ; 
            plot(world_t, world_v_contrast, '-x', 'DisplayName', 'Measured'); 
            plot(world_t, world_fit, '-x', 'DisplayName', 'Fit'); 
            
            xlabel("Time [t]");
            ylabel("Contrast"); 

            legend show; 
        %}

        % Add to the list of local amplitudes
        local_amplitudes(nn) = world_amplitude;
    end

    % Save the mean amplitude for this contrast level
    avg_amplitudes(cc) = mean(local_amplitudes);
end

% Plot the contrast gamma function
figure ;
title(sprintf("World Contrast Gamma Function at %d Hz, ND%d",frequency,NDF));
hold on;
plot(contrast_levels, avg_amplitudes, '.k','MarkerSize',15, 'DisplayName', 'Observed');
hold on

% Add a linear fit
p = polyfit(contrast_levels, avg_amplitudes,1);
plot(contrast_levels,polyval(p,contrast_levels),'-r');
plot([0 1],[0 1],':k');

xlabel("Contrast Level");
ylabel("Avg Amplitude");
xlim([-0.05 1.05]);
ylim([-0.05 1.05]);
axis square

end