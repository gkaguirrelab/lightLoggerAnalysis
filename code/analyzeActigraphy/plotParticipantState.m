function plotParticipantState(IMUdata, eyeStateData, blinkData, winSizeSec)
    % plotParticipantState: Overlays ENMO Activity and Eye Openness
    % ENMO (Euclidean Norm Minus One) is standard for processing raw accelerometer data. 
    % Mathematically simple and correlated with energy expenditure.
    
    %% Calculate Motion Magnitude
    %Convert to seconds and start at 0
    timeSecIMU = (IMUdata.timestamp_ns_ - IMUdata.timestamp_ns_(1)) / 1e9;
    timeMinIMU = timeSecIMU / 60;

    % Calculate sampling frequency (Fs) to convert winSizeSec to number of samples
    dt = mean(diff(timeSecIMU));
    winSizeSamples = round(winSizeSec / dt);
    
    %Calculate the vector magnitude (L2 Norm)
    % L2 Norm is euclidean magnitude (sqrt of each direction squared and
    % % summed)
    magAccel = sqrt(sum(IMUdata{:, {'accelerationX_g_', 'accelerationY_g_', 'accelerationZ_g_'}}.^2, 2));

    % Apply the ENMO
    % "Minus 1" is subtracting 1g becuase of gravity. 
    %  Rectify to get rid of noise when still and we care about energy
    %  against gravity most.
    enmo = max(0, magAccel - 1);

    % Calculate activity counts (mean ENMO per window)
    activityIndex = movmean(enmo, winSizeSamples); %movmean returns the local k-point mean values, where each
    % mean is calculated over a sliding window of length k across neighboring
    % elements of A.

    %% Process eye openness data
    % align eye timestamps to the same T0 as the IMU
    timeSecEye = (double(eyeStateData.timestamp_ns_) - double(IMUdata.timestamp_ns_(1))) / 1e9;
    timeMinEye = timeSecEye / 60;
    dtEye = mean(diff(timeSecEye));
    winSizeSamplesEye = round(winSizeSec / dtEye);

    % Mean left and right eyes
    avgAperture = mean(eyeStateData{:, {'eyelidApertureLeft_mm_', 'eyelidApertureRight_mm_'}}, 2, 'omitnan');

    % Nan blink intervals
    isBlink = false(size(avgAperture));

    % Loop through each recorded blink and flag those time points
    for i = 1:height(blinkData)
        % Find indices in eyeStateData that fall within this blink's start/end
        blinkMask = (eyeStateData.timestamp_ns_ >= blinkData.startTimestamp_ns_(i)) & ...
                    (eyeStateData.timestamp_ns_ <= blinkData.endTimestamp_ns_(i));
        
        isBlink = isBlink | blinkMask;
    end

    % Set aperture to NaN during those identified intervals
    avgAperture(isBlink) = NaN;
    %smooth the data with a rolling mean becuase it is noisy
    smoothAperture = movmean(avgAperture, winSizeSamplesEye, 'omitnan');

    % Plot 
    figure('Color', 'w');
    hold on;
    
    % Left Axis: Physical Movement (Smoothed)
    yyaxis left
    pMove = plot(timeMinIMU, activityIndex, '-', 'LineWidth', 1.0, 'Color', [0 0.45 0.74]);
    ylabel(['Physical Activity (ENMO, g)']);
    ylim([0, max(activityIndex)*1.2]);

    % Right Axis: Eye Aperture (Raw)
    yyaxis right
    % Plot a small line at the bottom for every blink start
    timeMinBlinks = (double(blinkData.startTimestamp_ns_) - double(IMUdata.timestamp_ns_(1))) / 1e9 / 60;
    yTickPos = zeros(size(timeMinBlinks)); 
    hBlink = stem(timeMinBlinks, yTickPos + 0.2, 'Color', [0.7 0.7 0.7], ...
             'Marker', 'none', 'LineWidth', 0.5);

    pEyeSmooth = plot(timeMinEye, smoothAperture, 'Color', [0.85 0.33 0.1], 'LineWidth', 0.5, 'LineStyle', '-'); 
    ylabel('Mean Eyelid Aperture (mm)');
    ylim([0, max(smoothAperture, [], 'omitnan')*1.1]);

    % Match the axis color to the line for clarity
    ax = gca;
    ax.YAxis(2).Color = [0.85 0.33 0.1];
    ax.YAxis(1).Color = [0 0.45 0.74];

    xlabel('Time (min)');
    title(['Participant State  - ' num2str(winSizeSec) 's Sliding Window']);
    legend([pMove, pEyeSmooth, hBlink], {'Movement (ENMO)', 'Smoothed Eye Aperture', 'Blinks'}, ...
        'Location', 'best');
end