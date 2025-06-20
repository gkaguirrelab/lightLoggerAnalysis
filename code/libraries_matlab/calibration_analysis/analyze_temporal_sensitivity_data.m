function ttfFigHandle = analyze_temporal_sensitivity_data(CalibrationData, sorted_measurements, NDF, ttfFigHandle)

% We measured the modulation depth of the stimulus source as a function of
% contrast using the Klein photometer and the "ChromaSurf" sofwate. The
% values are given here, and used for correction of the observed response
% amplitude below
testModDepthFrequencies = [0.1, 0.25,0.5,1,3,6,12,25,50,100];
testModDepthAmplitudes = [1,1,1,1,1,0.983,0.975,0.939,0.90,0.66];

% Retrieve the contrast levels and frequencies shown to the device
contrast_levels = CalibrationData.contrast_levels;
frequencies = CalibrationData.frequencies;
n_measures = CalibrationData.n_measures;

% Assert that we have readings for every frequency and contrast level
intended_size = [numel(contrast_levels), numel(frequencies), n_measures];
actual_size = size(sorted_measurements);
assert(isequal(intended_size, actual_size));

% Initialize a map between the contrast levels
% and the mean amplituded of the repeated measures
% at each frequency
world_mean_measure_amplitude_by_contrast = containers.Map(contrast_levels, zeros([1, numel(frequencies)]));

world_amplitude_data = zeros(n_measures,numel(frequencies));

% Sorted data is based into here as a 3D cell array in the shape
% (contrastIdx, frequencyIdx, measurementIdx)

% Create a figre for the TTF across contrast levels
if isempty(ttfFigHandle)
    ttfFigHandle = figure;
else
    figure(ttfFigHandle);
end
hold on ;

% Now, we will iterate over the contrast levels and frequencies
for cc = 1:numel(contrast_levels)
    % Retrieve the contrast level for this recording
    contrast = contrast_levels(cc);

    % Iterate over the measurements at this contrast and frequency
    % Generate a matrix to save the amplitudes per measure per sensor
    local_amplitudes = zeros([n_measures, 4]);
    for nn = 1:n_measures
        % Generate a figure for the sinusoidal fits
        fitsByMeasureHandle = figure('Position',[0,0,500,1500]);
        t = tiledlayout(numel(frequencies),1);
        sgtitle(sprintf("Sensor fits by measure. C: %.3f | M: %d", contrast, nn));

        % Iterate over the frequency levels at this contrast
        for ff = 1:numel(frequencies)
            % Get this frequency
            frequency = frequencies(ff);

            % Retrieve the measurement for this combination
            % of contrast and frequency and measurement number
            measurement = sorted_measurements{cc, ff, nn};

            % Initialize variables to hold the flattened measurements
            world_t = measurement.W.t; % t is time units converted to seconds
            world_v = measurement.W.v; % v here is the mean pixel of each frame

            % Now, we will perform the fitting procedure of the observed measurements to
            % a fit sinusoid. Though first, convert the world to contrast
            world_v_contrast = (world_v - mean(world_v)) / mean(world_v);
            [world_r2, world_amplitude, ~, world_fit] = fourierRegression( world_v_contrast, world_t, frequency );

            % Show the world fit for this combination of contrast measure and frequency
            nexttile(t) ;

            % Plot the world camera and its fit sinusoid
            local_world_t = world_t - world_t(1);

            plot(local_world_t, world_v_contrast, 'k.', 'DisplayName', sprintf('WorldMeasured Freq=%2.1f Hz',frequency));
            hold on
            plot(local_world_t, world_fit, 'r-', 'DisplayName', sprintf('WorldFit R2=%.3f', world_r2));
            xlabel("Time [s]");
            ylabel("Contrast");
            ylim([-1, 1]); % Set to 8 bit unsigned range

            % Show the legend for this tile
            % legend show;
            title(sprintf("F: %2.2f | R2: %2.2f", frequency, world_r2));

            % Save the amplitude for this measure
            local_amplitudes(nn, 1) = world_amplitude;
            world_amplitude_data(nn,ff) = world_amplitude;

        end

        % Generate a figure to illustrate AGC settings performance
        AGCByMeasureHandle = figure('Position',[0,0,500,1500]);
        t2 = tiledlayout(numel(frequencies), 1);
        sgtitle(sprintf("AGC performance by measure. C: %.3f | M: %d", contrast, nn));

        for ff = 1:numel(frequencies)
            % Get this frequency
            frequency = frequencies(ff);

            % Retrieve the measurement for this combination
            % of contrast and frequency and measurement number
            measurement = sorted_measurements{cc, ff, nn};

            % Retreive the world settings for this measure 
            world_t = measurement.W.t; 
            world_settings = measurement.W.settings;

            % Show the world AGC performance for this measure 
            nexttile(t2); 
        
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
        end 

    end

    % Update the TTF plot with the data for this contrast level
    figure(ttfFigHandle);

    % Obtain and plot the mean amplitude of response for each frequency.
    % Adjust for the contrast of the stimulus and for the roll-off in
    % modulation depth with temporal frequency.
    % testModDepthFrequencies = [0.25,5,3,6,12,25,50,100];
    meanAmpData = mean(world_amplitude_data,1).*(1./testModDepthAmplitudes)./(contrast_levels(cc));
    plot(log10(frequencies), meanAmpData,...
        '.k',...
        'MarkerSize',15,...
        'LineWidth',2,...
        'DisplayName', sprintf("C: %.3f", contrast_levels(cc)));

    % Generate the ideal device curve for the world camera
    world_fps = 120;
    xfine = logspace(log10(frequencies(1)),log10(frequencies(end)),100);
    ideal_device = idealDiscreteSampleFilter(xfine, 1/world_fps);
    plot(log10(xfine), ideal_device, "--", "Color",[0.5 0.5 0.5], "DisplayName", "Ideal Device");

    % Get the approximate filter frequency and amplitude for the data
    [filterFreqHz,filterAmp] = approxFreqFilter(frequencies,meanAmpData);
    fitVals = filterAmp*idealDiscreteSampleFilter(xfine, 1/filterFreqHz);
    plot(log10(xfine), fitVals, "-", "Color",[1 0 0], "DisplayName", "Data fit");

    % Add a title
    title(sprintf("World TTF | NDF: %.3f | Approx freq %2.2f Hz", NDF, filterFreqHz));

end

% Label and clean up the plot
xlabel("Frequency [hz]");
ylabel("Relative Contrast Response");
legend show;
a = gca();
a.XTick = log10(frequencies);
a.XTickLabels = arrayfun(@(x) num2str(x),frequencies,'UniformOutput',false);

end

%% LOCAL FUNCTIONS

function [filterFreqHz,filterAmp] = approxFreqFilter(frequencies,amplitudes)

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