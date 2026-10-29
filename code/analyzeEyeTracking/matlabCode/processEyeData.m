% =========================================================================
% Pupil Labs Data Cleaning, Imputation, and Kinematic Smoothing Routine
% =========================================================================

% 1. Load the data files while preserving the original CSV column names
opts = detectImportOptions('3d_eye_states.csv');
opts.VariableNamingRule = 'preserve';
eye_data = readtable('3d_eye_states.csv', opts);

opts_gaze = detectImportOptions('alternative_camera_gaze.csv');
opts_gaze.VariableNamingRule = 'preserve';
gaze_data = readtable('alternative_camera_gaze.csv', opts_gaze);

opts_blinks = detectImportOptions('blinks.csv');
opts_blinks.VariableNamingRule = 'preserve';
blinks_data = readtable('blinks.csv', opts_blinks);

% 2. Define the Master Timeline (using gaze as the lower-resolution base)
gaze_t_ns = gaze_data.("timestamp [ns]");
eye_t_ns = eye_data.("timestamp [ns]");

% Restrict the timeline to the overlapping domain of both sensors to prevent 
% interpolation artifacts (wild deflections) at the start or end of the data.
min_valid_t_ns = max(gaze_t_ns(1), eye_t_ns(1));
max_valid_t_ns = min(gaze_t_ns(end), eye_t_ns(end));
valid_time_idx = (gaze_t_ns >= min_valid_t_ns) & (gaze_t_ns <= max_valid_t_ns);

t_ns = gaze_t_ns(valid_time_idx);
t_sec = (t_ns - t_ns(1)) / 1e9; % Convert to seconds for reference/plotting

% Extract and interpolate bilateral vectors from 3d_eye_states onto the valid gaze timebase
pupil_L = interp1(eye_t_ns, eye_data.("pupil diameter left [mm]"), t_ns, 'linear');
pupil_R = interp1(eye_t_ns, eye_data.("pupil diameter right [mm]"), t_ns, 'linear');
aperture_L = interp1(eye_t_ns, eye_data.("eyelid aperture left [mm]"), t_ns, 'linear');
aperture_R = interp1(eye_t_ns, eye_data.("eyelid aperture right [mm]"), t_ns, 'linear');

% 3. Extract Gaze Data (Native to this timebase, bound to the overlapping domain)
gaze_x = gaze_data.("gaze x [px]")(valid_time_idx);
gaze_y = gaze_data.("gaze y [px]")(valid_time_idx);

% First-stage spatial bounds filtering: NaN any coordinates outside the sensor frame
out_of_bounds_gaze = (gaze_x < 0) | (gaze_x > 640) | (gaze_y < 0) | (gaze_y > 480);
gaze_x(out_of_bounds_gaze) = NaN;
gaze_y(out_of_bounds_gaze) = NaN;

% 4. Asymmetric Blink Masking 
is_blink = false(size(t_ns));
is_blink_raw = false(size(t_ns)); % Store the unpadded blink vector
pad_ns_pre  = 100 * 1e6; % 100 ms pre-blink padding to catch the downstroke
pad_ns_post = 150 * 1e6; % 150 ms post-blink padding to catch the slower upstroke

for i = 1:height(blinks_data)
    raw_start = blinks_data.("start timestamp [ns]")(i);
    raw_end   = blinks_data.("end timestamp [ns]")(i);
    
    % Record the raw, unpadded blinks for visualization
    is_blink_raw = is_blink_raw | (t_ns >= raw_start & t_ns <= raw_end);
    
    % Apply padding for the actual removal mask
    start_t = raw_start - pad_ns_pre;
    end_t   = raw_end + pad_ns_post;
    is_blink = is_blink | (t_ns >= start_t & t_ns <= end_t);
end

% 5. Validity Bounds (Remove physically impossible sizes independently)
invalid_L = (pupil_L < 1.5) | (pupil_L > 8.5);
invalid_R = (pupil_R < 1.5) | (pupil_R > 8.5);
pupil_L(invalid_L) = NaN;
pupil_R(invalid_R) = NaN;

% Combine Eyes for Pupil and generate temporary combined aperture for spike detection
pupil_combined = mean([pupil_L, pupil_R], 2, 'omitnan');
aperture_combined_temp = mean([aperture_L, aperture_R], 2, 'omitnan');

% --- SAVE RAW VECTORS FOR PLOTTING BEFORE MASKING ---
pupil_raw = pupil_combined;
aperture_L_raw = aperture_L;
aperture_R_raw = aperture_R;
aperture_mean_raw = mean([aperture_L_raw, aperture_R_raw], 2, 'omitnan');
aperture_diff_raw = aperture_L_raw - aperture_R_raw;
gaze_x_raw = gaze_x;
gaze_y_raw = gaze_y;

% 6. Kinematic Spike Detection and Independent Mask Dilation
% Calculate frame-to-frame absolute derivatives for physiological states
diff_pupil = [0; abs(diff(pupil_combined))];
diff_aperture = [0; abs(diff(aperture_combined_temp))];

% Calculate frame-to-frame absolute derivatives (velocity) for spatial gaze
diff_gaze_x = [0; abs(diff(gaze_x))];
diff_gaze_y = [0; abs(diff(gaze_y))];
diff_gaze_mag = sqrt(diff_gaze_x.^2 + diff_gaze_y.^2);

% Thresholds adjusted for typical lower-res gaze sample rates
is_physio_spike = (diff_pupil > 0.25) | (diff_aperture > 0.3);
is_gaze_spike = diff_gaze_mag > 40; % Threshold for anomalous gaze jumps (> 30 px/frame)

% Dilate the respective masks using movmax
dilation_window_samples = 5; 
is_blink_dilated = movmax(is_blink, dilation_window_samples) > 0;
is_physio_spike_dilated = movmax(is_physio_spike, dilation_window_samples) > 0;
is_gaze_spike_dilated = movmax(is_gaze_spike, dilation_window_samples) > 0;

% Combine for the master physiological removal mask
physio_mask = is_blink_dilated | is_physio_spike_dilated;

% Targeted Gaze Scrutiny: Only censor gaze if it exhibits a high derivative spike 
% *during* a padded blink region. This prevents over-filtering stable tracking.
gaze_mask = is_blink_dilated & is_gaze_spike_dilated;

% Censor variables using their respective domain masks
pupil_combined(physio_mask) = NaN;
aperture_L(physio_mask) = NaN;
aperture_R(physio_mask) = NaN;
gaze_x(gaze_mask) = NaN;
gaze_y(gaze_mask) = NaN;

% =========================================================================
% Prevent end-of-vector pchip extrapolation spikes
% Find the last valid observation and carry it forward to overwrite trailing NaNs
last_valid_p = find(~isnan(pupil_combined), 1, 'last');
if ~isempty(last_valid_p) && last_valid_p < length(pupil_combined)
    pupil_combined(last_valid_p+1:end) = pupil_combined(last_valid_p);
end

last_valid_aL = find(~isnan(aperture_L), 1, 'last');
if ~isempty(last_valid_aL) && last_valid_aL < length(aperture_L)
    aperture_L(last_valid_aL+1:end) = aperture_L(last_valid_aL);
end

last_valid_aR = find(~isnan(aperture_R), 1, 'last');
if ~isempty(last_valid_aR) && last_valid_aR < length(aperture_R)
    aperture_R(last_valid_aR+1:end) = aperture_R(last_valid_aR);
end

% Prevent start-of-vector pchip extrapolation spikes
% Find the first valid observation and carry it backward to overwrite leading NaNs
first_valid_p = find(~isnan(pupil_combined), 1, 'first');
if ~isempty(first_valid_p) && first_valid_p > 1
    pupil_combined(1:first_valid_p-1) = pupil_combined(first_valid_p);
end

first_valid_aL = find(~isnan(aperture_L), 1, 'first');
if ~isempty(first_valid_aL) && first_valid_aL > 1
    aperture_L(1:first_valid_aL-1) = aperture_L(first_valid_aL);
end

first_valid_aR = find(~isnan(aperture_R), 1, 'first');
if ~isempty(first_valid_aR) && first_valid_aR > 1
    aperture_R(1:first_valid_aR-1) = aperture_R(first_valid_aR);
end
% =========================================================================

% 7. Targeted Imputation (With Physiological Safety Limits)
% Sample gaps adjusted assuming lower sampling rate (e.g. 60-120Hz)
max_general_gap_samples = 25;  
max_blink_gap_samples = 200;   

% Pass 1: Interpolate for blinks, but CAP at the physiological blink maximum (2 seconds).
pupil_full = fillmissing(pupil_combined, 'pchip', 'MaxGap', max_blink_gap_samples);
aperture_L_full = fillmissing(aperture_L, 'pchip', 'MaxGap', max_blink_gap_samples);
aperture_R_full = fillmissing(aperture_R, 'pchip', 'MaxGap', max_blink_gap_samples);
gaze_x_full = fillmissing(gaze_x, 'nearest', 'MaxGap', max_blink_gap_samples);
gaze_y_full = fillmissing(gaze_y, 'nearest', 'MaxGap', max_blink_gap_samples);

% Pass 2: Apply the strict interpolation for general noise or tracking loss (250 ms).
pupil_clean = fillmissing(pupil_combined, 'pchip', 'MaxGap', max_general_gap_samples);
aperture_L_clean = fillmissing(aperture_L, 'pchip', 'MaxGap', max_general_gap_samples);
aperture_R_clean = fillmissing(aperture_R, 'pchip', 'MaxGap', max_general_gap_samples);
gaze_x_clean = fillmissing(gaze_x, 'nearest', 'MaxGap', max_general_gap_samples);
gaze_y_clean = fillmissing(gaze_y, 'nearest', 'MaxGap', max_general_gap_samples);

% Pass 3: Explicitly overwrite the targeted mask periods with the fully interpolated data.
% Physiological variables overwrite entire blink regions
pupil_clean(is_blink_dilated) = pupil_full(is_blink_dilated);
aperture_L_clean(is_blink_dilated) = aperture_L_full(is_blink_dilated);
aperture_R_clean(is_blink_dilated) = aperture_R_full(is_blink_dilated);

% Gaze variables overwrite only the isolated artifact spikes detected during blinks
gaze_x_clean(gaze_mask) = gaze_x_full(gaze_mask);
gaze_y_clean(gaze_mask) = gaze_y_full(gaze_mask);

% 8. Final Kinematic and State Smoothing
window_size = 3; % Adjusted window for lower frequency data

% Median filter for physiological states (preserves tonic baseline shifts)
pupil_final = medfilt1(pupil_clean, window_size, 'omitnan', 'truncate');
aperture_L_final = medfilt1(aperture_L_clean, window_size, 'omitnan', 'truncate');
aperture_R_final = medfilt1(aperture_R_clean, window_size, 'omitnan', 'truncate');

% Calculate Mean and Difference vectors from the independent cleaned signals
aperture_mean_final = mean([aperture_L_final, aperture_R_final], 2, 'omitnan');
aperture_diff_final = aperture_L_final - aperture_R_final;

% Savitzky-Golay filter for spatial gaze kinematics (preserves saccadic peak velocity)
gaze_x_final = smoothdata(gaze_x_clean, 'sgolay', window_size);
gaze_y_final = smoothdata(gaze_y_clean, 'sgolay', window_size);

% Output matrix combining the final cleaned vectors
cleaned_data_matrix = [t_sec, pupil_final, aperture_mean_final, aperture_diff_final, gaze_x_final, gaze_y_final];

% =========================================================================
% 9. Visualization: Raw vs. Cleaned Data Vectors
% =========================================================================

figure('Name', 'Eye Data Processing: Raw vs Cleaned', 'Position', [100, 100, 800, 1400]);

% Panel 1: Pupil Diameter
subplot(6,1,1); hold on;
plot(t_sec, pupil_raw, 'k.');
plot(t_sec, pupil_final, 'r-', 'LineWidth', 1);
title('Pupil Diameter');
ylabel('Diameter [mm]');
legend('Raw', 'Cleaned', 'Location', 'best');
grid on; hold off;

% Panel 2: Eyelid Aperture (Mean)
subplot(6,1,2); hold on;
plot(t_sec, aperture_mean_raw, 'k.');
plot(t_sec, aperture_mean_final, 'r-', 'LineWidth', 1);
title('Mean Eyelid Aperture (L/R Combined)');
ylabel('Aperture [mm]');
grid on; hold off;

% Panel 3: Eyelid Aperture (Difference)
subplot(6,1,3); hold on;
plot(t_sec, aperture_diff_raw, 'k.');
plot(t_sec, aperture_diff_final, 'r-', 'LineWidth', 1);
title('Eyelid Aperture Difference (Left - Right)');
ylabel('Difference [mm]');
grid on; hold off;

% Panel 4: Gaze X
subplot(6,1,4); hold on;
plot(t_sec, gaze_x_raw, 'k.');
plot(t_sec, gaze_x_final, 'r-', 'LineWidth', 1);
title('Gaze X');
ylabel('Position [px]');
grid on; hold off;

% Panel 5: Gaze Y
subplot(6,1,5); hold on;
plot(t_sec, gaze_y_raw, 'k.');
plot(t_sec, gaze_y_final, 'r-', 'LineWidth', 1);
title('Gaze Y');
ylabel('Position [px]');
grid on; hold off;

% Panel 6: Raw Blinks
subplot(6,1,6); hold on;
plot(t_sec, is_blink_raw, 'k-', 'LineWidth', 1.5);
title('Raw Blink Events');
xlabel('Time [s]');
ylabel('Blink Active');
ylim([-0.1 1.1]);
grid on; hold off;