function plotParticipantState(raw_dir, processing_dir, output_dir, subject_id, activity, figureTitle, IMUdata, eyeStateData, blinkData, gazeData, options)
% Generate multi-panel summary figures of participant actigraphy, gaze, and light data
%
% Syntax:
%   plotParticipantState(raw_dir, processing_dir, output_dir, subject_id, activity, figureTitle, IMUdata, eyeStateData, blinkData, gazeData)
%   plotParticipantState(..., 'save_figures', true)
%
% Description:
%   Produces four diagnostic figures for a single participant and activity
%   session: (1) a five-panel time-series summary showing ENMO activity,
%   rotational velocity, gaze angles, eye state (aperture, pupil, blinks),
%   and illuminance; (2) gaze-vs-head correlation scatter plots; (3)
%   activity-vs-eye-state correlations; and (4) illuminance-vs-eye-state
%   correlations. Sensor data are temporally aligned using Neon timestamps,
%   world camera timestamps, and a configurable world-to-AS temporal offset.
%
% Inputs:
%   raw_dir               - String. Path to the raw recording directory.
%   processing_dir        - String. Path to the processing directory.
%   output_dir            - String. Path to the output directory for saved
%                           figures.
%   subject_id            - String. Subject identifier (e.g., 'FLIC_2001').
%   activity              - String. Activity name (e.g., 'walkOutdoor').
%   figureTitle           - String. Title string used in correlation figure
%                           titles.
%   IMUdata               - Table or char/string. IMU sensor data table or
%                           path to CSV file.
%   eyeStateData          - Table or char/string. Eye state data table or
%                           path to CSV file.
%   blinkData             - Table or char/string. Blink event data table or
%                           path to CSV file.
%   gazeData              - Table or char/string. Gaze angle data table or
%                           path to CSV file.
%
% Optional key/value pairs:
%   'save_figures'        - Logical (default: false). Export figures as PDF
%                           and EPS to output_dir.
%   'winSizeSec'          - Scalar double (default: 5). Sliding window
%                           size in seconds for smoothing signals.
%   'force_recalc'        - Logical (default: false). Force reload of
%                           Python utility libraries.
%   'w_as_offset'         - Scalar double (default: -139.354). World-to-AS
%                           temporal offset in milliseconds.
%   'active_threshold'    - Scalar double (default: 0.025). ENMO threshold
%                           for active/inactive classification.
%   'outdoor_threshold'   - Scalar double (default: 10^2.5). Illuminance
%                           threshold in lux for outdoor classification.
%
% Outputs:
%   none
%
% Examples:
%{
    plotParticipantState(raw_dir, processing_dir, output_dir, ...
        'FLIC_2001', 'walkOutdoor', 'Walk Outdoor', ...
        IMUdata, eyeStateData, blinkData, gazeData, ...
        'save_figures', true);
%}
    arguments
        raw_dir;
        processing_dir;    
        output_dir; 
        subject_id; 
        activity; 
        figureTitle; 

        IMUdata, 
        eyeStateData,
        blinkData, 
        gazeData, 
        
        options.save_figures = false; 
        options.winSizeSec = 5;  
        options.force_recalc = false; 
        options.w_as_offset = -139.354; % milliseconds
        options.active_threshold = 0.025;
        options.outdoor_threshold = 10 ^ 2.5;

    end
    winSizeSec = options.winSizeSec; 

    %% Load Utility Libraries & Data
    persistent world_util ms_util;
    if isempty(world_util) || isempty(ms_util) || options.force_recalc
        world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));
        ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path")); 
    end
    
    % Load in the actigraphy data if needed 
    [IMUdata, eyeStateData, blinkData, gazeData] = load_actigraphy_data(IMUdata, eyeStateData, blinkData, gazeData);
    
    % First, we need to find and load in the timestamp data so we can do temporal alignment
    path_to_raw_recording = fullfile(raw_dir, subject_id, activity, "GKA"); 
    path_to_activity_start_end = fullfile(processing_dir, subject_id, activity, "tag_task_start_end.mat"); 
    path_to_gka_world_timestamps_neon_time = fullfile(processing_dir, subject_id, activity, "Neon", "egocentric_mapper_results", "alternative_camera_timestamps.csv"); 

    assert(exist(path_to_raw_recording, 'dir') ~= 0, ...
    sprintf('Path %s does not exist', path_to_raw_recording));

    assert(exist(path_to_gka_world_timestamps_neon_time, 'file') ~= 0, ...
    sprintf('Path %s does not exist', path_to_gka_world_timestamps_neon_time));

    % At this point, everything is in Neon time
    [world_t, ms_t, ms_v] = load_ms_world_timestamps_and_data(path_to_raw_recording, path_to_gka_world_timestamps_neon_time, world_util, ms_util); 

    % With this information, we can also get the start/ending time of the activity in neon time 
    tag_task_start_end_frames = load(path_to_activity_start_end).tag_task_start_end; 

    % Clip the tag start end to the world t length 
    % This is due to rare rounding error in the neon alignment that clips off 
    % some small number of frames sometimes 
    tag_start_end_frames = tag_task_start_end_frames.tag; 
    task_start_end_frames = tag_task_start_end_frames.task; 
    tag_start_end_frames  = max(min(tag_start_end_frames,  numel(world_t)), 1);
    task_start_end_frames = max(min(task_start_end_frames, numel(world_t)), 1);
    tag_start_end_neon_time = world_t(tag_start_end_frames); 
    task_start_end_neon_time = world_t(task_start_end_frames);


    % We know the temporal offset (W-AS) is options.w_as_offset milliseconds. 
    % Therefore, we apply this correction now to the Nanosecond time 
    % Given that it is negative, it implies that the world is BEHIND 
    % the MS, so subtract this offset from the MS to bring them onto the same 
    % time
    offset_in_nanoseconds = options.w_as_offset * 1e6; 
    ms_t = ms_t + offset_in_nanoseconds; 

    % Convert MS counts to absolute illuminance
    MS2illum_lux = msCounts2Illuminance(ms_v);
    meanLux = mean(MS2illum_lux, 'omitnan');
    fprintf('%.6f\n', meanLux);
    
    % Classify outdoor/indoor based on this 
    outdoor_threshold = options.outdoor_threshold; 
    outdoor_indoor_classifications = ~classifyIndoorOutdoorPeriods(path_to_raw_recording, ...
                                                                  "window_size_seconds", options.winSizeSec,...
                                                                  "outdoor_threshold", outdoor_threshold...
                                                                  ); 

    %% Calculate Time Offsets (T0 based on IMU start)
    t0 = IMUdata.timestamp_ns_(1);
    
    % Align MS Timing to T0 (Minutes)
    timeMinMS = (double(ms_t) - double(t0)) / 1e9 / 60;
    timeMinTask = (double(task_start_end_neon_time) - double(t0)) / 1e9 / 60;
    
    % IMU Timing
    timeMinIMU = (double(IMUdata.timestamp_ns_) - double(t0)) / 1e9 / 60;
    dt = mean(diff(timeMinIMU * 60)); 
    winSizeSamples = round(winSizeSec / dt);
    
    % Gaze/Eye Timing
    timeMinGaze = (double(gazeData.timestamp_ns_) - double(t0)) / 1e9 / 60;
    timeMinEye = (double(eyeStateData.timestamp_ns_) - double(t0)) / 1e9 / 60;
    dtEye = mean(diff(timeMinEye * 60));
    winSizeSamplesEye = round(winSizeSec / dtEye);
    
    %% Signal Processing
    % ENMO
    magAccel = sqrt(sum(IMUdata{:, {'accelerationX_g_', 'accelerationY_g_', 'accelerationZ_g_'}}.^2, 2));
    enmo = max(0, magAccel - 1); 
    activityIndex = movmean(enmo, winSizeSamples, 'omitnan'); 
    active_threshold = options.active_threshold; 
    active_inactive_classifications = classifyActiveInactivePeriods(IMUdata, ...
                                                                    "window_size_seconds", options.winSizeSec,...
                                                                    "active_threshold", active_threshold...
                                                                   );  % Classify active/inactive based on this 

    % Rotational Velocities
    unwrappedYaw = unwrap(deg2rad(IMUdata.yaw_deg_)) * (180/pi);
    vRoll  = gradient(IMUdata.roll_deg_, dt);
    vPitch = gradient(IMUdata.pitch_deg_, dt);
    vYaw   = gradient(unwrappedYaw, dt);
    
    % Eye State & Blink Masking
    avgAperture = mean(eyeStateData{:, {'eyelidApertureLeft_mm_', 'eyelidApertureRight_mm_'}}, 2, 'omitnan');
    avgPupil    = mean(eyeStateData{:, {'pupilDiameterLeft_mm_', 'pupilDiameterRight_mm_'}}, 2, 'omitnan');
    gazeElev = gazeData.elevation_deg_;
    gazeAzim = gazeData.azimuth_deg_;

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
    
    avgAperture(isBlinkEye) = NaN; avgPupil(isBlinkEye) = NaN;
    gazeElev(isBlinkGaze) = NaN; gazeAzim(isBlinkGaze) = NaN;
    
    smoothAperture = movmean(avgAperture, winSizeSamplesEye, 'omitnan');
    smoothPupil    = movmean(avgPupil, winSizeSamplesEye, 'omitnan');

    %% Colors
    cRoll  = [0 0.4470 0.7410]; cPitch = [0.8500 0.3250 0.0980]; cYaw = [0.9290 0.6940 0.1250];
    cPupil = [0 0.5 0.5]; cApert = [0.6350 0.0780 0.1840]; cLux = [0.4660 0.6740 0.1880]; % Green
    cOutdoor = [1 1 0];   % yellow
    cActive  = [0 1 0];   % green
    cTaskBounds = [0 0.65 1];
    summaryAxisFontSize = 11;
    summaryLabelFontSize = 13;
    summaryTitleFontSize = 24;
    summaryLegendFontSize = 9;
    
    %----------------------------------------------
    %% Figure 1: Session Summary (Updated to 5 rows)
    figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 14.5 8.32], 'Name', 'Summary');
    tlo1 = tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    tlo1.OuterPosition = [0.05 0.06 0.66 0.90];
    
    % Subplot 1: ENMO
    ax1 = nexttile(tlo1); 
    hold(ax1, 'on');

    formattedTaskTitle = format_activity_title(string(activity));
    pENMO = plot(ax1, timeMinIMU, activityIndex, 'k', 'LineWidth', 1.2, 'DisplayName', 'ENMO');
    ylabel(ax1, 'Activity (g)', 'FontSize', summaryLabelFontSize); 
    title(ax1, formattedTaskTitle, 'FontSize', summaryTitleFontSize, 'FontWeight', 'bold');
    ylim(ax1, [0, 0.2]);
    grid(ax1, 'on');
    ax1.FontSize = summaryAxisFontSize;
    ax1.XTickLabel = [];


    % Subplot 2: Rotational Vel
    ax2 = nexttile(tlo1); hold(ax2, 'on');
    pRoll = plot(ax2, timeMinIMU, vRoll,  'Color', [cRoll  0.3], 'LineWidth', 1.2, 'DisplayName', 'Roll');
    pPitch = plot(ax2, timeMinIMU, vPitch, 'Color', [cPitch 0.3], 'LineWidth', 1.2, 'DisplayName', 'Pitch');
    pYaw = plot(ax2, timeMinIMU, vYaw,   'Color', [cYaw   0.3], 'LineWidth', 1.2, 'DisplayName', 'Yaw');
    ylabel(ax2, 'Rot. Vel (deg/s)', 'FontSize', summaryLabelFontSize); grid on;
    ylim([-360, 360]);
    ax2.FontSize = summaryAxisFontSize;
    ax2.XTickLabel = [];
    lgd2 = legend(ax2, [pRoll, pPitch, pYaw], {'Roll', 'Pitch', 'Yaw'}, 'FontSize', summaryLegendFontSize, 'Location', 'northeastoutside');
    lgd2.Box = 'off';
    

    % Subplot 3: Gaze
    ax3 = nexttile(tlo1); hold(ax3, 'on');
    pElev = plot(ax3, timeMinGaze, gazeElev, 'Color', cPitch, 'LineWidth', 1.2, 'DisplayName', 'Elev');
    pAzim = plot(ax3, timeMinGaze, gazeAzim, 'Color', cYaw, 'LineWidth', 1.2, 'DisplayName', 'Azim');
    ylabel(ax3, 'Gaze (deg)', 'FontSize', summaryLabelFontSize); grid on;
    ylim([-60, 60]);
    ax3.FontSize = summaryAxisFontSize;
    ax3.XTickLabel = [];
    lgd3 = legend(ax3, [pElev, pAzim], {'Elev', 'Azim'}, 'FontSize', summaryLegendFontSize, 'Location', 'northeastoutside');
    lgd3.Box = 'off';
    
    % Subplot 4: Eye State
    ax4 = nexttile(tlo1); hold(ax4, 'on');
    yyaxis left
    timeMinBlinks = (double(blinkData.startTimestamp_ns_) - double(t0)) / 1e9 / 60;
    pBlink = stem(ax4, timeMinBlinks, zeros(size(timeMinBlinks)) + 2.0, 'Color', [0.6 0.6 0.6], 'Marker', 'none', 'DisplayName', 'Blinks');
    pAperture = plot(ax4, timeMinEye, smoothAperture, '-', 'Color', cApert, 'LineWidth', 1.2, 'DisplayName', 'Aperture'); 
    ylabel('Eyelid Aperture (mm)', 'FontSize', summaryLabelFontSize); ax4.YAxis(1).Color = cApert; 
    ylim(ax4, [0, 15]);
    yyaxis right
    pPupil = plot(ax4, timeMinEye, smoothPupil, '-', 'Color', cPupil, 'LineWidth', 1.2, 'DisplayName', 'Pupil');
    ylim(ax4, [2, 6]);
    ylabel('Pupil (mm)', 'FontSize', summaryLabelFontSize); ax4.YAxis(2).Color = cPupil; 
    grid on;
    ax4.FontSize = summaryAxisFontSize;
    ax4.XTickLabel = [];
    lgd4 = legend(ax4, [pAperture, pPupil, pBlink], {'Eyelid openness', 'Pupil size', 'Blinks'}, 'FontSize', summaryLegendFontSize, 'Location', 'northeastoutside');
    lgd4.Box = 'off';

    % Subplot 5: Absolute Illuminance (Lux)
    ax5 = nexttile(tlo1);
    hold(ax5, 'on');

    pLux = plot(ax5, timeMinMS, MS2illum_lux, '.', 'Color', cLux, 'LineWidth', 1.1, 'DisplayName', 'Illum');

    ax5.YScale = 'log';
    ax5.YMinorGrid = 'off'; 
    ax5.YGrid = 'on';      
    ylabel(ax5, 'Illum (Lux)', 'FontSize', summaryLabelFontSize); 
    xlabel(ax5, 'Time (min)', 'FontSize', summaryLabelFontSize);
    ylim(ax5, [1, 10e4]); 
    grid(ax5, 'on');
    ax5.FontSize = summaryAxisFontSize;

    linkaxes([ax1, ax2, ax3, ax4, ax5], 'x');
    % Find the maximum time across your data to set a tight limit
    maxTime = max([timeMinIMU; timeMinGaze; timeMinEye; timeMinMS]);
    xlim(ax1, [0, maxTime]);
    add_task_boundary_lines([ax1, ax2, ax3, ax4, ax5], timeMinTask, cTaskBounds);
    drawnow;

    xl1 = xlim(ax1);
    yl1 = ylim(ax1);

    hActive = patch(ax1, ...
        [xl1(1) xl1(2) xl1(2) xl1(1)], ...
        [active_threshold active_threshold yl1(2) yl1(2)], ...
        cActive, ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none', ...
        'DisplayName', 'Active');

    uistack(hActive, 'bottom');

    plot(ax1, xl1, [active_threshold active_threshold], ...
        'Color', cActive, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    pTaskPeriod1 = plot(ax1, nan, nan, '-', 'Color', cTaskBounds, 'LineWidth', 1.8, 'DisplayName', 'Task period');
    lgd1 = legend(ax1, [pENMO, hActive, pTaskPeriod1], {'ENMO', 'Active threshold', 'Task period'}, 'FontSize', summaryLegendFontSize, 'Location', 'northeastoutside');
    lgd1.Box = 'off';


    xl5 = xlim(ax5);
    yl5 = ylim(ax5);

    hOutdoor = patch(ax5, ...
        [xl5(1) xl5(2) xl5(2) xl5(1)], ...
        [outdoor_threshold outdoor_threshold yl5(2) yl5(2)], ...
        cOutdoor, ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none', ...
        'DisplayName', 'Outdoor');

    uistack(hOutdoor, 'bottom');

    plot(ax5, xl5, [outdoor_threshold outdoor_threshold], ...
        'Color', cOutdoor, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    lgd5 = legend(ax5, [pLux, hOutdoor], {'Illum', 'Outdoor threshold'}, 'FontSize', summaryLegendFontSize, 'Location', 'northeastoutside');
    lgd5.Box = 'off';
    
    % Keep the legends at MATLAB's native northeastoutside placement.
    lgd1.Location = 'northeastoutside';
    lgd2.Location = 'northeastoutside';
    lgd3.Location = 'northeastoutside';
    lgd4.Location = 'northeastoutside';
    lgd5.Location = 'northeastoutside';
    drawnow;

    % Move all summary legends except the eyelid/pupil legend further to
    % the right by a fixed amount in inches.
    legend_shift_in = 0.4;
    lgd1.Units = 'inches'; lgd1.Position(1) = lgd1.Position(1) + legend_shift_in;
    lgd2.Units = 'inches'; lgd2.Position(1) = lgd2.Position(1) + legend_shift_in;
    lgd3.Units = 'inches'; lgd3.Position(1) = lgd3.Position(1) + legend_shift_in;
    lgd5.Units = 'inches'; lgd5.Position(1) = lgd5.Position(1) + legend_shift_in;

    % Save the figure if desired 
    if(options.save_figures)
        export_figure_dual(gcf, output_dir, "actigraphy_summary");
        close(gcf); 
    end 

    %----------------------------------------------
    %% Figure 2: Gaze-Head Correlation
    targetTime = double(gazeData.timestamp_ns_);
    imuTime = double(IMUdata.timestamp_ns_);

    interpPitch = interp1(imuTime, IMUdata.pitch_deg_, targetTime, 'nearest', NaN);
    interpYaw   = interp1(imuTime, unwrappedYaw, targetTime, 'nearest', NaN);
    validP = find(~isnan(interpPitch) & ~isnan(gazeElev));
    validY = find(~isnan(interpYaw) & ~isnan(gazeAzim));

    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.52 0.52 0.34], 'Name', 'Gaze-Head Corr');
    tlo2 = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    axA = nexttile(tlo2); hold(axA, 'on');
    scatter(axA, interpPitch(validP), gazeElev(validP), 10, cPitch, 'filled',...
        'Marker', 'o', ...
        'MarkerEdgeAlpha', 0.05, ...
        'MarkerFaceAlpha', 0.05);
    p1 = polyfit(interpPitch(validP), gazeElev(validP), 1);
    plot(axA, interpPitch(validP), polyval(p1, interpPitch(validP)), 'k', 'LineWidth', 1.5);
    title(sprintf('Vertical (R=%.2f)', corr(interpPitch(validP), gazeElev(validP))));
    ylim([-60, 60]);
    grid on; axis square; xlabel('Pitch'); ylabel('Elev');

    axB = nexttile(tlo2); hold(axB, 'on');
    scatter(axB, interpYaw(validY), gazeAzim(validY), 10, cYaw, 'filled',...
        'Marker', 'o', ...
        'MarkerEdgeAlpha', 0.05, ...
        'MarkerFaceAlpha', 0.05);
    p2 = polyfit(interpYaw(validY), gazeAzim(validY), 1);
    plot(axB, interpYaw(validY), polyval(p2, interpYaw(validY)), 'k', 'LineWidth', 1.5);
    title(sprintf('Horizontal (R=%.2f)', corr(interpYaw(validY), gazeAzim(validY))));
    grid on; axis square; xlabel('Yaw'); ylabel('Azim');
    title(tlo2, ['Gaze Pos vs Head Pos: ' figureTitle], 'FontWeight', 'bold');
    ylim([-60, 60]);
    drawnow;
    % Save the figure if desired 
    if(options.save_figures)
        export_figure_dual(gcf, output_dir, "gaze_head_correlation");
        close(gcf); 
    end 


    %----------------------------------------------
    %% Figure 3: Activity vs. Eye State Correlation
    eyeTime = double(eyeStateData.timestamp_ns_);
    interpAperture = interp1(eyeTime, smoothAperture, imuTime, 'nearest', NaN);
    interpPupil    = interp1(eyeTime, smoothPupil, imuTime, 'nearest', NaN);
    
    validA = find(~isnan(activityIndex) & ~isnan(interpAperture));
    validPup = find(~isnan(activityIndex) & ~isnan(interpPupil));

    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.10 0.52 0.34], 'Name', 'Activity-Eye Corr');
    tlo3 = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    
    axC = nexttile(tlo3); hold(axC, 'on');
    if ~isempty(validA)
        scatter(axC, activityIndex(validA), interpAperture(validA), 10, cApert,...
            'Marker', 'o', ...
            'MarkerEdgeAlpha', 0.05, ...
            'MarkerFaceAlpha', 0.05);
        pf = polyfit(activityIndex(validA), interpAperture(validA), 1);
        plot(axC, activityIndex(validA), polyval(pf, activityIndex(validA)), 'k', 'LineWidth', 1.5);
        title(sprintf('Aperture (R=%.2f)', corr(activityIndex(validA), interpAperture(validA))));
    end
    xlabel('Activity (ENMO)'); ylabel('Aperture (mm)'); grid on; axis square;

    axD = nexttile(tlo3); hold(axD, 'on');
    if ~isempty(validPup)
        scatter(axD, activityIndex(validPup), interpPupil(validPup), 10, cPupil, 'filled',...
            'Marker', 'o', ...
            'MarkerEdgeAlpha', 0.05, ...
            'MarkerFaceAlpha', 0.05);
        pf = polyfit(activityIndex(validPup), interpPupil(validPup), 1);
        plot(axD, activityIndex(validPup), polyval(pf, activityIndex(validPup)), 'k', 'LineWidth', 1.5);
        title(sprintf('Pupil (R=%.2f)', corr(activityIndex(validPup), interpPupil(validPup))));
    end
    xlabel('Activity (ENMO)'); ylabel('Pupil (mm)'); grid on; axis square;
    
    title(tlo3, ['Activity vs Eye State: ' figureTitle], 'FontWeight', 'bold');
    drawnow;
    % Save the figure if desired 
    if(options.save_figures)
        export_figure_dual(gcf, output_dir, "activity_vs_eyestate");
        close(gcf); 
    end 

    %----------------------------------------------
    %% Figure 4: Illuminance vs. Eye State Correlation
    % Use log10 because the pupillary light reflex is logarithmic
    logLux = log10(MS2illum_lux + eps); 
    interpLogLux_Eye = interp1(double(ms_t), logLux, eyeTime, 'nearest', NaN);

    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.10 0.58 0.36], 'Name', 'Illuminance Correlations');
    tlo4 = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
    
    % Panel 1: Log(Lux) vs. Eyelid Aperture
    axE = nexttile(tlo4); hold(axE, 'on');
    validLuxA = ~isnan(interpLogLux_Eye) & ~isnan(smoothAperture);
    if any(validLuxA)
        scatter(axE, interpLogLux_Eye(validLuxA), smoothAperture(validLuxA), 10, cApert, 'filled',...
            'Marker', 'o', ...
            'MarkerEdgeAlpha', 0.05, ...
            'MarkerFaceAlpha', 0.05);
        pf3 = polyfit(interpLogLux_Eye(validLuxA), smoothAperture(validLuxA), 1);
        plot(axE, interpLogLux_Eye(validLuxA), polyval(pf3, interpLogLux_Eye(validLuxA)), 'k', 'LineWidth', 1.5);
        title(axE, sprintf('Lux vs Eye Openness (R=%.2f)', corr(interpLogLux_Eye(validLuxA), smoothAperture(validLuxA))));
    end
    xlabel(axE, 'Log10(Lux)'); ylabel(axE, 'Eye Openness (mm)'); grid on; axis square;
    ylim([0, 15]); 

    % Panel 2: Log(Lux) vs. Pupil Diameter
    axF = nexttile(tlo4); hold(axF, 'on');
    validLuxP = ~isnan(interpLogLux_Eye) & ~isnan(smoothPupil);
    if any(validLuxP)
        scatter(axF, interpLogLux_Eye(validLuxP), smoothPupil(validLuxP), 10, cPupil, 'filled',...
            'Marker', 'o', ...
            'MarkerEdgeAlpha', 0.05, ...
            'MarkerFaceAlpha', 0.05);
        pf4 = polyfit(interpLogLux_Eye(validLuxP), smoothPupil(validLuxP), 1);
        plot(axF, interpLogLux_Eye(validLuxP), polyval(pf4, interpLogLux_Eye(validLuxP)), 'k', 'LineWidth', 1.5);
        title(axF, sprintf('Lux vs Pupil (R=%.2f)', corr(interpLogLux_Eye(validLuxP), smoothPupil(validLuxP))));
    end
    xlabel(axF, 'Log10(Lux)'); ylabel(axF, 'Pupil (mm)'); grid on; axis square;
    ylim([2, 6]);
    
    title(tlo4, ['Environmental Light vs Eye State: ' figureTitle], 'FontWeight', 'bold');
    drawnow;
        % Save the figure if desired 
    if(options.save_figures)
        export_figure_dual(gcf, output_dir, "illuminance_vs_eyestate");
        close(gcf); 
    end 

end

% Load in the world timestamps 
function world_t = load_world_gka_neon_timestamps(path)
% Internal helper to load world gka neon timestamps.
%
% Syntax:
%   world_t = load_world_gka_neon_timestamps(path)
%
% Description:
%   This local helper function internal helper to load world gka neon timestamps within its parent workflow.
% Inputs:
%   path                     - Path-like input used by the function.
%
% Outputs:
%   world_t                  - Output produced by the function.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');

    % Force timestamp column to int64 (preserves ns precision)
    opts = setvartype(opts, 'timestamp [ns]', 'uint64');

    world_timestamps_table = readtable(path, opts); 
    world_t = uint64(world_timestamps_table.('timestamp [ns]')); 
    
    return 

end 


function [world_t, ms_t, ms_v] = load_ms_world_timestamps_and_data(path_to_raw_recording,  path_to_gka_world_timestamps_neon_time, world_util, ms_util)
% Internal helper to load ms world timestamps and data.
%
% Syntax:
%   world_t, ms_t, ms_v = load_ms_world_timestamps_and_data(path_to_raw_recording, path_to_gka_world_timestamps_neon_time, world_util, ms_util)
%
% Description:
%   This local helper function internal helper to load ms world timestamps and data within its parent workflow.
% Inputs:
%   path_to_raw_recording    - Path-like input used by the function.
%   path_to_gka_world_timestamps_neon_time - Path-like input used by the function.
%   world_util               - Input used by the function.
%   ms_util                  - Input used by the function.
%
% Outputs:
%   world_t                  - Output produced by the function.
%   ms_t                     - Output produced by the function.
%   ms_v                     - Output produced by the function.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    gka_world_timestamps_gka_time = int64( world_util.world_timestamps_from_chunks(path_to_raw_recording, py.False) )'; % nanoseconds  
    ms_data_and_timestamps = cell(ms_util.ms_data_from_chunks(path_to_raw_recording));
    ms_data = double(ms_data_and_timestamps{1}); 
    ms_timestamps_gka_time = int64(ms_data_and_timestamps{2} * (10 ^ 9))'; % nanoseconds 

    %load egocentric video mapper timestamps for converting between gka and
    %neon
    gka_world_timestamps_neon_time = int64(load_world_gka_neon_timestamps(path_to_gka_world_timestamps_neon_time)); 

    % Assert these two have the same length
    num_gka_world_timestamps_gka_time = numel(gka_world_timestamps_gka_time);
    num_gka_world_timestamps_neon_time = numel(gka_world_timestamps_neon_time);
    num_frames_diff = abs(num_gka_world_timestamps_gka_time - num_gka_world_timestamps_neon_time); 
    if(num_gka_world_timestamps_gka_time ~= num_gka_world_timestamps_neon_time)
        warning("Num timestamps GKA time: %d ~= num timestamps Neon time: %d\n", num_gka_world_timestamps_gka_time, num_gka_world_timestamps_neon_time);

        % If the difference is a rounding error, just take the minimum number existing
        if(num_frames_diff > 10)
            error("Difference is %d frames. Likely programmatic error somewhere", num_frames_diff); 
        end
        
        warning("Difference is =< 10 frames. Cropping to minimum size");

        if(num_gka_world_timestamps_gka_time < num_gka_world_timestamps_neon_time)
            gka_world_timestamps_neon_time = gka_world_timestamps_neon_time(1:num_gka_world_timestamps_gka_time);
        else    
            num_gka_world_timestamps_gka_time = gka_world_timestamps_gka_time(1:num_gka_world_timestamps_neon_time);
        end 



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


% Local function to do loading in of the needed .csv files 
function [IMUdata, eyeStateData, blinkData, gazeData] = load_actigraphy_data(IMUdata, eyeStateData, blinkData, gazeData)
% Internal helper to load actigraphy data.
%
% Syntax:
%   IMUdata, eyeStateData, blinkData, gazeData = load_actigraphy_data(IMUdata, eyeStateData, blinkData, gazeData)
%
% Description:
%   This local helper function internal helper to load actigraphy data within its parent workflow.
% Inputs:
%   IMUdata                  - Input used by the function.
%   eyeStateData             - Input used by the function.
%   blinkData                - Input used by the function.
%   gazeData                 - Input used by the function.
%
% Outputs:
%   IMUdata                  - Output produced by the function.
%   eyeStateData             - Output produced by the function.
%   blinkData                - Output produced by the function.
%   gazeData                 - Output produced by the function.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    if(isstring(IMUdata) || ischar(IMUdata))
        IMUdata = readtable(IMUdata);
    end 

    if(isstring(eyeStateData) || ischar(eyeStateData))
        eyeStateData = readtable(eyeStateData);
    end 

    if(isstring(blinkData) || ischar(blinkData))
        blinkData = readtable(blinkData);
    end 

    if(isstring(gazeData) || ischar(gazeData))
        gazeData = readtable(gazeData);
    end 

    return ; 
end 

function add_task_boundary_lines(axesHandles, timeMinTask, lineColor)
% Internal helper to add task boundary lines.
%
% Syntax:
%   add_task_boundary_lines(axesHandles, timeMinTask, lineColor)
%
% Description:
%   This local helper function internal helper to add task boundary lines within its parent workflow.
% Inputs:
%   axesHandles              - Input used by the function.
%   timeMinTask              - Input used by the function.
%   lineColor                - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    if(isempty(timeMinTask) || numel(timeMinTask) ~= 2 || any(~isfinite(timeMinTask)))
        return;
    end

    for ax = axesHandles
        hold(ax, 'on');
        xline(ax, timeMinTask(1), '-', 'Color', lineColor, 'LineWidth', 1.8, 'HandleVisibility', 'off');
        xline(ax, timeMinTask(2), '-', 'Color', lineColor, 'LineWidth', 1.8, 'HandleVisibility', 'off');
    end
end

function export_figure_dual(figHandle, output_dir, basename)
% Internal helper to export figure dual.
%
% Syntax:
%   export_figure_dual(figHandle, output_dir, basename)
%
% Description:
%   This local helper function internal helper to export figure dual within its parent workflow.
% Inputs:
%   figHandle                - Input used by the function.
%   output_dir               - Path-like input used by the function.
%   basename                 - Input used by the function.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    pdf_path = fullfile(output_dir, basename + ".pdf");
    eps_path = fullfile(output_dir, basename + ".eps");

    exportgraphics(figHandle, pdf_path, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white', ...
        'Resolution', 300);

    exportgraphics(figHandle, eps_path, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white', ...
        'Resolution', 300);
end

function formatted_title = format_activity_title(activity_name)
% Internal helper to format activity title.
%
% Syntax:
%   formatted_title = format_activity_title(activity_name)
%
% Description:
%   This local helper function internal helper to format activity title within its parent workflow.
% Inputs:
%   activity_name            - Input used by the function.
%
% Outputs:
%   formatted_title          - Output produced by the function.
%
% Examples:
%{
    % See plotParticipantState.m for usage context.
%}

    formatted_title = regexprep(char(activity_name), '([a-z])([A-Z])', '$1 $2');
    formatted_title(1) = upper(formatted_title(1));
end
