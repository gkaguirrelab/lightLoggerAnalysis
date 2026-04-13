function classifications = classifyIndoorOutdoorPeriods()
    arguments
        raw_recording_dir;
        options.force_recalc = false; 
        options.window_size_seconds = 5; 
    end 

    %% Load Utility Libraries & Data
    persistent ms_util;
    if isempty(ms_util) || options.force_recalc
        ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path")); 
    end

    % Load in the timestamps and values of the MS 
    % for the given recording 
    [ms_t, ms_v] = load_ms_data(raw_recording_dir, ms_util);
    num_timestamps = numel(ms_t); 

    % Preallocate the return variable that classifies each timepoint 
    % as indoor or outdoor 
    classifications = nan(num_tiemstamps);

    % Given our timestamps are returned in nanoseconds, 
    % convert the window size seconds to nanoseconds too 
    winows_size_nanoseconds = options.window_size_nanseconds * (10 ^ 9); 



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