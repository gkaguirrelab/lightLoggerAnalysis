function delayByLineSecs = generateRollingShutterPhaseMap(CalibrationData, sorted_measurements)

%{
LightLoggerCalibrationData = ndf0_meanedCols;
CalibrationData = LightLoggerCalibrationData.CalibrationData;
parsed_readings = LightLoggerCalibrationData.ParsedReadings;
delayByLineSecs = generateRollingShutterPhaseMap(CalibrationData.temporal_sensitivity,parsed_readings.temporal_sensitivity);
%}

% Retrieve the contrast levels and frequencies shown to the device
contrast_levels = CalibrationData.contrast_levels;
frequencies = CalibrationData.frequencies;
n_measures = CalibrationData.n_measures;

% Assert that we have readings for every frequency and contrast level
intended_size = [numel(contrast_levels), numel(frequencies), n_measures];
actual_size = size(sorted_measurements);
assert(isequal(intended_size, actual_size));

% Initialize a map between the contrast levels and the mean amplitude of
% the repeated measures at each frequency
world_mean_measure_amplitude_by_contrast = containers.Map(contrast_levels, zeros([1, numel(frequencies)]));

world_amplitude_data = zeros(n_measures,numel(frequencies));

% Sorted data is based into here as a 3D cell array in the shape
% (contrastIdx, frequencyIdx, measurementIdx)


% Figure out how many vertical lines are in each image
tmp = sorted_measurements{1, 1, 1};
nLines = size(tmp.W.v,2);

% Now, we will iterate over the contrast levels and frequencies
for cc = 1:numel(contrast_levels)

    % Retrieve the contrast level for this recording
    contrast = contrast_levels(cc);

    % Define a data matrix that will hold the temporal delays calculated across
    % frequencies and measures
    temporalDelayByLineSecs = zeros(n_measures, numel(frequencies),nLines);
    temporalDelayByFrameSecs = zeros(n_measures, numel(frequencies));

    % Iterate over the measurements at this contrast and frequency
    for nn = 1:n_measures

        % Iterate over the frequency levels at this contrast
        for ff = 1:numel(frequencies)

            % Get this frequency
            frequency = frequencies(ff);

            % Retrieve the measurement for this combination
            % of contrast and frequency and measurement number
            measurement = sorted_measurements{cc, ff, nn};

            % Initialize variables to hold the flattened measurements
            world_t = measurement.W.t; % t is time units converted to seconds
            world_v = measurement.W.v; % v here is the mean pixel of each frame for each line

            % Obtain and store the mean offset for the frame
            world_v_frame = mean(world_v,2);
            world_v_contrast = (world_v_frame - mean(world_v_frame)) / mean(world_v_frame);
            [r2, ~, meanPhase, ] = fourierRegression( world_v_contrast, world_t, frequency );
            if r2<0.9
                warning('bad fit')
            end
            % Convert phase to temporal delay and store the value
            temporalDelayByFrameSecs(nn,ff) = (meanPhase/(2*pi))/frequency;


            % Loop through the lines and obtain the phase
            for ll = 1:nLines
                world_v_line = world_v(:,ll);
                world_v_contrast = (world_v_line - mean(world_v_line)) / mean(world_v_line);
                [r2, ~, phase, ] = fourierRegression( world_v_contrast, world_t, frequency );

                % Get the phase difference
                phaseDiff = phase- meanPhase;

                % Wrap phase to ±pi/2
                if phaseDiff > pi/2
                    phaseDiff = phaseDiff - pi;
                end

                % Convert phase to temporal delay relative to the frame mean and store the value
                temporalDelayByLineSecs(nn,ff,ll) = ((phaseDiff/(2*pi))/frequency);
            end

        end

    end

end

% Calculate and return the mean delay by line. Just use the first 8
% frequencies for this, as the higher frequencies produce weird phase
% wraps.
tmp = temporalDelayByLineSecs(:,1:8,:);
tmp = mean(tmp,1); tmp = mean(tmp,2);
delayByLineSecs = squeeze(tmp);

% Make a figure
figure
plot(delayByLineSecs*1e3,'-r');
hold on
plot(delayByLineSecs*1e3,'.k');
xlim([0 480]);
axis square
ylabel('Temporal offset [ms]');
xlabel('Image line index');
box off
a = gca();
a.TickDir = 'out';
a.XTick = 0:80:480;

figure
plot([230 250],[0,0],':k');
hold on
plot([240 240],[-0.3 0.3],':k');
plot(delayByLineSecs*1e3,'-r');
plot(delayByLineSecs*1e3,'.k','MarkerSize',15);
xlim([230 250]);
axis square
ylabel('Temporal offset [ms]');
xlabel('Image line index');
box off
a = gca();
a.TickDir = 'out';
a.XTick = 230:5:250;

% Show the source sinusoid at the hyper resolved time
ff = 3; nn = 1;
measurement = sorted_measurements{cc, ff, nn};
world_t = measurement.W.t; % t is time units converted to seconds
world_v = measurement.W.v; % v here is the mean pixel of each frame for each line

% Get the mean difference in amplitude between the odd and even lines. This
% happens because the even lines contain the blue sensors, while the odd
% lines have just red and green
oddIdx = 1:2:nLines;
evenIdx = 2:2:nLines;
oddLineMean = mean(mean(world_v(:,oddIdx)));
evenLineMean = mean(mean(world_v(:,evenIdx)));
lineAmpAdjust = [...
    (oddLineMean / mean([oddLineMean evenLineMean])) ...
    (evenLineMean / mean([oddLineMean evenLineMean])) ...
    ];

delayByLineSecs = linspace(-(1/400),(1/400),480);

lineWeight = squeeze(mean(world_v,1));

% Loop through the frames and lines and generate the hyper resolved version
t_hyper = [];
v_hyper = [];
lineIdx = [];
for tt = 1:length(world_t)
    for ll = 1:nLines
        v_hyper(end+1) = world_v(tt,ll) / lineWeight(ll);
        t_hyper(end+1) = world_t(tt) + delayByLineSecs(ll);
        lineIdx(end+1) = mod(ll,2);
    end
end
lineIdx = logical(lineIdx);
figure
plot(t_hyper(lineIdx),v_hyper(lineIdx),'.r');
hold on
plot(t_hyper(~lineIdx),v_hyper(~lineIdx),'.b');
%xlim([t_hyper(1),t_hyper(20000)]);



end