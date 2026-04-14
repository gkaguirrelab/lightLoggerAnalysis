function classifications = classifyActiveInactivePeriods()
    arguments 
        IMUdata 
        options.window_size_seconds = 5; 
        options.active_threshold = 0.05; % Anything above this value is considered active. 
    end 

    % We will accept either a path to the data or a table of the data 
    % if it's a path, read it in now 
    if(isstring(IMUdata) || ischar(IMUdata))
        IMUdata = readtable(IMUdata); 
    end 

    % Get the table size to kow how many readings we have 
    num_readings = size(IMUdata, 1); 

    % Take the magnitude of acceleration in all dimensions (vector norm)
    magAccel = sqrt(sum(IMUdata{:, {'accelerationX_g_', 'accelerationY_g_', 'accelerationZ_g_'}}.^2, 2));
    enmo = max(0, magAccel - 1); 
    activityIndex = movmean(enmo, winSizeSamples, 'omitnan');  


    % Preallocate the return variable that classifies each timepoint 
    % as active or inactive 
    classifications = zeros(num_readings, 1);
    outdoor_readings = activityIndex > options.active_threshold; 
    classifications(outdoor_readings) = 1; 

    return ; 
end 
