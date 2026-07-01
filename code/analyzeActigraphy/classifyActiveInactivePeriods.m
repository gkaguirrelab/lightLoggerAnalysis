function classifications = classifyActiveInactivePeriods(IMUdata, options)
% Classify each IMU timepoint as active or inactive based on accelerometer data
%
% Syntax:
%   classifications = classifyActiveInactivePeriods(IMUdata)
%   classifications = classifyActiveInactivePeriods(IMUdata, 'window_size_seconds', 10)
%
% Description:
%   Given IMU accelerometer data (either as a table or a path to a CSV
%   file), compute the Euclidean Norm Minus One (ENMO) activity index
%   using a sliding window average, then classify each timepoint as
%   active (1) or inactive (0) based on a threshold.
%
% Inputs:
%   IMUdata               - Table or char/string. Either a table of IMU
%                           readings with columns 'timestamp_ns_',
%                           'accelerationX_g_', 'accelerationY_g_',
%                           'accelerationZ_g_', or a path to a CSV file
%                           containing such a table.
%
% Optional key/value pairs:
%   'window_size_seconds' - Scalar double (default: 5). Size of the
%                           sliding window in seconds for smoothing the
%                           ENMO signal.
%   'active_threshold'    - Scalar double (default: 0.05). ENMO values
%                           above this threshold are classified as active.
%
% Outputs:
%   classifications       - Numeric vector. Binary vector (0 or 1) of the
%                           same length as the number of IMU readings,
%                           where 1 indicates active and 0 indicates
%                           inactive.
%
% Examples:
%{
    IMUdata = readtable('imu.csv');
    classifications = classifyActiveInactivePeriods(IMUdata, 'active_threshold', 0.03);
%}
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

    % Calculate the window size
    t0 = IMUdata.timestamp_ns_(1);
    timeMinIMU = (double(IMUdata.timestamp_ns_) - double(t0)) / 1e9 / 60;
    dt = mean(diff(timeMinIMU * 60)); 
    winSizeSamples = round(options.window_size_seconds / dt);

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
