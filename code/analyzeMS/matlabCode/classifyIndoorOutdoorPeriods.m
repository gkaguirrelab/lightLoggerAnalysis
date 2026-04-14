function classifications = classifyIndoorOutdoorPeriods()
    arguments
        raw_recording_dir;
        options.force_recalc = false; 
        options.window_size_seconds = 5; 
        options.outdoor_threshold = 10 ^ 2.5; % This lux value and above is outdoor  
    end 

    %% Load Utility Libraries & Data
    persistent ms_util;
    if isempty(ms_util) || options.force_recalc
        ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path")); 
    end

    % Load in the timestamps and values of the MS 
    % for the given recording 
    [ms_t, ms_v_raw] = load_ms_data(raw_recording_dir, ms_util);
    num_timestamps = numel(ms_t); 

    % Convert ms sensor values to illuminance 
    ms_v_illum = msCounts2Illuminance(ms_v_raw);

    % Preallocate the return variable that classifies each timepoint 
    % as indoor or outdoor 
    classifications = zeros(num_timestamps, 1);

    % Given our timestamps are returned in nanoseconds, 
    % convert the window size seconds to nanoseconds too 
    window_size_nanoseconds = options.window_size_seconds * (10 ^ 9); 

    % Because the FPS is non-constant, we CANNOT use movemean. 
    % We will implement our own sliding window technique 
    counts_movemean = msTimestampSlidingWindow(ms_t, ms_v_illum, window_size_nanoseconds, "mean"); 

    % Now, let's find the portions that are below the outdoor threshold 
    indoor_frames = counts_movemean < options.outdoor_threshold; 
    classifications(indoor_frames) = 1; 

    return ; 
end 


% Local funciton to load in the MS timestamps and data 
function [ms_t, ms_v] = load_ms_data(raw_recording_dir, ms_util) 
    ms_data_and_timestamps = cell(ms_util.ms_data_from_chunks(raw_recording_dir));
    ms_data = double(ms_data_and_timestamps{1}); 
    ms_timestamps_gka_time = int64(ms_data_and_timestamps{2} * (10 ^ 9))'; % nanoseconds 

    % Save MS values (v) and timestamps (t)
    ms_v = ms_data; 
    ms_t = ms_timestamps_gka_time; 

end 