function plotParticipantState(IMUdata, eyeStateData, blinkData, gazeData, winSizeSec, activityName, options)
    arguments 
        IMUdata, 
        eyeStateData,
        blinkData, 
        gazeData, 
        winSizeSec, 
        activityName,
        options.input_dir; 
        options.output_dir; 
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {}; 
        options.force_recalc = false; 
    end

    % plotParticipantState: Stacked subplots with synchronized color coding

    % Load in the Python utility libraries if not loaded before     
    persistent world_util ms_util;
    if isempty(world_util) || isempty(ms_util) || options.force_recalc
        world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));
        ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path")); 
    end

    % Load in the world and ms timestamps in neon timespace
    path_to_raw_recording = "/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026/FLIC_2001/chat/GKA"
    [world_t, ms_t, ms_v] = load_ms_world_timestamps_and_data(path_to_raw_recording, world_util, ms_util); 

    disp(size(world_t))

    return 

    %% Calculate Time Offsets (T0 based on IMU start)
    t0 = IMUdata.timestamp_ns_(1);
    
    % IMU Timing
    timeMinIMU = (double(IMUdata.timestamp_ns_) - double(t0)) / 1e9 / 60;
    dt = mean(diff(timeMinIMU * 60)); % sampling interval in seconds
    winSizeSamples = round(winSizeSec / dt);
    
    % Gaze/Eye Timing
    timeMinGaze = (double(gazeData.timestamp_ns_) - double(t0)) / 1e9 / 60;
    timeMinEye = (double(eyeStateData.timestamp_ns_) - double(t0)) / 1e9 / 60;
    dtEye = mean(diff(timeMinEye * 60));
    winSizeSamplesEye = round(winSizeSec / dtEye);
    
    %% Calculate Motion Magnitude (ENMO)
    magAccel = sqrt(sum(IMUdata{:, {'accelerationX_g_', 'accelerationY_g_', 'accelerationZ_g_'}}.^2, 2));
    enmo = max(0, magAccel - 1);
    activityIndex = movmean(enmo, winSizeSamples); 
    
    %% Calculate Rotational Velocities
    unwrappedYaw = unwrap(deg2rad(IMUdata.yaw_deg_)) * (180/pi);
    vRoll  = gradient(IMUdata.roll_deg_, dt);
    vPitch = gradient(IMUdata.pitch_deg_, dt);
    vYaw   = gradient(unwrappedYaw, dt);
    
    %% Process Eye State & Blink Masking
    avgAperture = mean(eyeStateData{:, {'eyelidApertureLeft_mm_', 'eyelidApertureRight_mm_'}}, 2, 'omitnan');
    avgPupil    = mean(eyeStateData{:, {'pupilDiameterLeft_mm_', 'pupilDiameterRight_mm_'}}, 2, 'omitnan');
    
    % Prepare Gaze Data for Masking
    gazeElev = gazeData.elevation_deg_;
    gazeAzim = gazeData.azimuth_deg_;

    % Create Blink Masks for Eye Data and Gaze Data
    isBlinkEye = false(size(avgAperture));
    for i = 1:height(blinkData)
        mask = (eyeStateData.timestamp_ns_ >= blinkData.startTimestamp_ns_(i)) & ...
               (eyeStateData.timestamp_ns_ <= blinkData.endTimestamp_ns_(i));
        isBlinkEye = isBlinkEye | mask;
    end
    
    isBlinkGaze = false(size(gazeElev));
    for i = 1:height(blinkData)
        mask = (gazeData.timestamp_ns_ >= blinkData.startTimestamp_ns_(i)) & ...
               (gazeData.timestamp_ns_ <= blinkData.endTimestamp_ns_(i));
        isBlinkGaze = isBlinkGaze | mask;
    end
    
    % Apply Masking (Set blinks to NaN)
    avgAperture(isBlinkEye) = NaN;
    avgPupil(isBlinkEye)    = NaN;
    gazeElev(isBlinkGaze)   = NaN;
    gazeAzim(isBlinkGaze)   = NaN;
    
    % Apply Windowing to Eye Data
    smoothAperture = movmean(avgAperture, winSizeSamplesEye, 'omitnan');
    smoothPupil    = movmean(avgPupil, winSizeSamplesEye, 'omitnan');

    %% Define Color Scheme
    cRoll  = [0 0.4470 0.7410];      % Blue
    cPitch = [0.8500 0.3250 0.0980]; % Orange
    cYaw   = [0.9290 0.6940 0.1250]; % Yellow
    cPupil = [0 0.5 0.5];            % Teal
    cApert = [0.6350 0.0780 0.1840]; % Maroon
    
    %% Figure 1: Session Summary
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.1 0.3 0.8], 'Name', 'Summary');
    tlo1 = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    ax1 = nexttile(tlo1);
    plot(ax1, timeMinIMU, activityIndex, 'k', 'LineWidth', 1.2, 'DisplayName', 'ENMO');
    ylabel(ax1, 'Activity (g)'); title(ax1, [activityName ' - ' num2str(winSizeSec) 's Window']);
    legend(ax1, 'Location', 'northeast', 'FontSize', 8); grid(ax1, 'on');

    ax2 = nexttile(tlo1); hold(ax2, 'on');
    plot(ax2, timeMinIMU, vRoll, 'Color', cRoll, 'DisplayName', 'Roll');
    plot(ax2, timeMinIMU, vPitch, 'Color', cPitch, 'DisplayName', 'Pitch');
    plot(ax2, timeMinIMU, vYaw, 'Color', cYaw, 'DisplayName', 'Yaw');
    ylabel(ax2, 'Rot. Vel (deg/s)'); legend(ax2, 'Location', 'northeast', 'FontSize', 8); grid(ax2, 'on');
    
    ax3 = nexttile(tlo1); hold(ax3, 'on');
    plot(ax3, timeMinGaze, gazeElev, 'Color', cPitch, 'DisplayName', 'Elev');
    plot(ax3, timeMinGaze, gazeAzim, 'Color', cYaw, 'DisplayName', 'Azim');
    ylabel(ax3, 'Gaze (deg)'); legend(ax3, 'Location', 'northeast', 'FontSize', 8); grid(ax3, 'on');
    
    ax4 = nexttile(tlo1); hold(ax4, 'on');
    yyaxis left
    timeMinBlinks = (double(blinkData.startTimestamp_ns_) - double(t0)) / 1e9 / 60;
    hBlink = stem(ax4, timeMinBlinks, zeros(size(timeMinBlinks)) + 2.0, 'Color', [0.6 0.6 0.6], 'Marker', 'none', 'DisplayName', 'Blinks');
    pAperture = plot(ax4, timeMinEye, smoothAperture, '-', 'Color', cApert, 'LineWidth', 1.2, 'DisplayName', 'Aperture'); 
    ylabel('Aperture (mm)'); ax4.YAxis(1).Color = cApert; 
    yyaxis right
    pPupil = plot(ax4, timeMinEye, smoothPupil, '-', 'Color', cPupil, 'LineWidth', 1.2, 'DisplayName', 'Pupil');
    ylabel('Pupil (mm)'); ax4.YAxis(2).Color = cPupil; 
    xlabel(ax4, 'Time (min)'); legend(ax4, [pAperture, pPupil, hBlink], 'Location', 'northeast', 'FontSize', 8); grid(ax4, 'on');
    
    linkaxes([ax1, ax2, ax3, ax4], 'x');
    drawnow;

    %% Figure 2: Gaze-Head Correlation
    targetTime = double(gazeData.timestamp_ns_);
    imuTime = double(IMUdata.timestamp_ns_);
    interpPitch = interp1(imuTime, IMUdata.pitch_deg_, targetTime, 'linear', 'extrap');
    interpYaw   = interp1(imuTime, unwrappedYaw, targetTime, 'linear', 'extrap');
    validP = find(~isnan(interpPitch) & ~isnan(gazeElev));
    validY = find(~isnan(interpYaw) & ~isnan(gazeAzim));
    
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.36 0.5 0.3 0.4], 'Name', 'Gaze-Head Corr');
    tlo2 = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    axA = nexttile(tlo2); hold(axA, 'on');
    scatter(axA, interpPitch(validP), gazeElev(validP), 10, cPitch, 'filled', 'MarkerFaceAlpha', 0.1);
    p1 = polyfit(interpPitch(validP), gazeElev(validP), 1);
    plot(axA, interpPitch(validP), polyval(p1, interpPitch(validP)), 'k', 'LineWidth', 1.5);
    title(sprintf('Vertical (R=%.2f)', corr(interpPitch(validP), gazeElev(validP))));
    grid on; axis square; xlabel('Pitch'); ylabel('Elev');

    axB = nexttile(tlo2); hold(axB, 'on');
    scatter(axB, interpYaw(validY), gazeAzim(validY), 10, cYaw, 'filled', 'MarkerFaceAlpha', 0.1);
    p2 = polyfit(interpYaw(validY), gazeAzim(validY), 1);
    plot(axB, interpYaw(validY), polyval(p2, interpYaw(validY)), 'k', 'LineWidth', 1.5);
    title(sprintf('Horizontal (R=%.2f)', corr(interpYaw(validY), gazeAzim(validY))));
    grid on; axis square; xlabel('Yaw'); ylabel('Azim');
    drawnow;

    %% Figure 3: Activity vs. Eye State Correlation
    eyeTime = double(eyeStateData.timestamp_ns_);
    interpAperture = interp1(eyeTime, smoothAperture, imuTime, 'linear', 'extrap');
    interpPupil    = interp1(eyeTime, smoothPupil, imuTime, 'linear', 'extrap');
    
    validA = find(~isnan(activityIndex) & ~isnan(interpAperture));
    validPup = find(~isnan(activityIndex) & ~isnan(interpPupil));

    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.36 0.05 0.3 0.4], 'Name', 'Activity-Eye Corr');
    tlo3 = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    
    axC = nexttile(tlo3); hold(axC, 'on');
    if ~isempty(validA)
        scatter(axC, activityIndex(validA), interpAperture(validA), 10, cApert, 'filled', 'MarkerFaceAlpha', 0.1);
        pf = polyfit(activityIndex(validA), interpAperture(validA), 1);
        plot(axC, activityIndex(validA), polyval(pf, activityIndex(validA)), 'k', 'LineWidth', 1.5);
        title(sprintf('Aperture (R=%.2f)', corr(activityIndex(validA), interpAperture(validA))));
    end
    xlabel('Activity (ENMO)'); ylabel('Aperture (mm)'); grid on; axis square;

    axD = nexttile(tlo3); hold(axD, 'on');
    if ~isempty(validPup)
        scatter(axD, activityIndex(validPup), interpPupil(validPup), 10, cPupil, 'filled', 'MarkerFaceAlpha', 0.1);
        pf = polyfit(activityIndex(validPup), interpPupil(validPup), 1);
        plot(axD, activityIndex(validPup), polyval(pf, activityIndex(validPup)), 'k', 'LineWidth', 1.5);
        title(sprintf('Pupil (R=%.2f)', corr(activityIndex(validPup), interpPupil(validPup))));
    end
    xlabel('Activity (ENMO)'); ylabel('Pupil (mm)'); grid on; axis square;
    
    title(tlo3, ['Activity vs Eye State: ' activityName], 'FontWeight', 'bold');
    drawnow;
end


% Load in the world timestamps 
function world_t = load_world_gka_neon_timestamps(path)
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');

    % Force timestamp column to int64 (preserves ns precision)
    opts = setvartype(opts, 'timestamp [ns]', 'uint64');

    world_timestamps_table = readtable(path, opts); 
    world_t = uint64(world_timestamps_table.('timestamp [ns]')); 
    
    return 

end 


function [world_t, ms_t, ms_v] = load_ms_world_timestamps_and_data(path_to_raw_recording, world_util, ms_util)
    % Load in the world timesstamps and the MS data/ timestamps in nanoseconds from the RAW recording
    %
    % Timestamps are w.r.t the start of the light logger device NOT the recording, 
    % so you will need to find a common start point to match
    
    % Load in the timestamps in GKA time 
    gka_world_timestamps_gka_time = int64( world_util.world_timestamps_from_chunks(path_to_raw_recording, py.False) )'; % nanoseconds  
    ms_data_and_timestamps = cell(ms_util.ms_data_from_chunks(path_to_raw_recording));
    ms_data = double(ms_data_and_timestamps{1}); 
    ms_timestamps_gka_time = int64(ms_data_and_timestamps{2} * (10 ^ 9))'; % nanoseconds 

    path_to_gka_world_timestamps_neon_time = "/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026/FLIC_2001/chat/Neon/egocentric_mapper_results/alternative_camera_timestamps.csv"; 
    gka_world_timestamps_neon_time = int64(load_world_gka_neon_timestamps(path_to_gka_world_timestamps_neon_time)); 

    % Assert these two have the same length
    num_gka_world_timestamps_gka_time = numel(gka_world_timestamps_gka_time);
    num_gka_world_timestamps_neon_time = numel(gka_world_timestamps_neon_time);
    if(num_gka_world_timestamps_gka_time ~= num_gka_world_timestamps_neon_time)
        error("Num timestamps GKA time: %d ~+ num timestamps Neon time: %d\n", num_gka_world_timestamps_gka_time, num_gka_world_timestamps_neon_time);
    end   
    
    % Now, let's do matching between the world timestamps and the MS 
    
    % First, let's find the nearest neighbor timestamp 
    % because it is very unlikely that timestamps will be 
    % exactly equal on the nanosecond level 
    double_ms_timestamps_gka_time = double(ms_timestamps_gka_time);
    double_world_timestamps_gka_time = double(gka_world_timestamps_gka_time); 
    anchorpoint_diff_magnitude = inf; 
    world_anchorpoint_idx = inf; 
    ms_anchorpoint_idx = inf; 
    for ii = 1:numel(gka_world_timestamps_gka_time)
        gka_world_current_timestamp = double_world_timestamps_gka_time(ii);
        [nearest_neighbor_diff, nearest_neighbor_idx] = min(abs(gka_world_current_timestamp - double_ms_timestamps_gka_time));
        if(nearest_neighbor_diff < anchorpoint_diff_magnitude)
            world_anchorpoint_idx = ii; 
            ms_anchorpoint_idx = nearest_neighbor_idx; 
            anchorpoint_diff_magnitude = nearest_neighbor_diff;
        end 
    end
    

    % Now we have the anchorpoint. We will set all other timestamps to be with respect to this 
    ms_t = zeros(numel(ms_timestamps_gka_time), 1, 'int64'); 
    anchorpoint_diff_with_directionality = ms_timestamps_gka_time(ms_anchorpoint_idx) - gka_world_timestamps_gka_time(world_anchorpoint_idx);
    ms_t(ms_anchorpoint_idx) = gka_world_timestamps_neon_time(world_anchorpoint_idx) + anchorpoint_diff_with_directionality; 

    % Now we do two loops. We do one from the anchorpoint to the 
    % start and one from the anchorpoint until the end to fill in these timestamps 

    % First, let's go backwards towards the start of the array
    for ii = ms_anchorpoint_idx - 1 : -1 : 1
        current_timestamp_gka_time = ms_timestamps_gka_time(ii);
        next_timestamp_gka_time = ms_timestamps_gka_time(ii + 1);
        next_timestamp_neon_time = ms_t(ii + 1);

        % Find the time delta between these two poinst 
        time_delta = next_timestamp_gka_time - current_timestamp_gka_time; 

        % Save the updated timestamp 
        ms_t(ii) = next_timestamp_neon_time - time_delta; 
    end 


    % Next, we will do a loop from the anchorpoint until the end of the array 
    for ii = ms_anchorpoint_idx + 1 : numel(ms_t)
        current_timestamp_gka_time = ms_timestamps_gka_time(ii);
        previous_timestamp_gka_time = ms_timestamps_gka_time(ii - 1); 
        previous_timestamp_neon_time = ms_t(ii - 1); 

        % Find the time delta between these points 
        time_delta = current_timestamp_gka_time - previous_timestamp_gka_time; 

        % Save the updated timestamp 
        ms_t(ii) = previous_timestamp_neon_time + time_delta; 
    end 

    % Return the converted values 
    world_t = gka_world_timestamps_neon_time;
    ms_v = ms_data; 

    return ; 

end 