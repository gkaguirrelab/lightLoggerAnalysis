function [minutes, seconds, ms] = frame2time(frame_num, fps)
% FRAME2TIME Converts a frame number (or vector of frame numbers) into a timestamp 
% (minutes, seconds, ms).
%
% Usage:
%   [m, s, ms] = frame2time(frame_number, fps)
%
% Inputs:
%   frame_num - The sequential number(s) of the frame. This can be a 
%                  scalar or a column vector.
%                  (Numeric vector/scalar)
%   fps          - Frames Per Second rate. Defaults to 120 if not provided.
%                  (Numeric scalar, must be > 0)
%
% Outputs:
%   minutes      - The whole number of minutes (Column vector).
%   seconds      - The whole number of seconds remaining after minutes (Column vector).
%   ms           - The whole number of milliseconds remaining after seconds
%                  (rounded to the nearest millisecond, Column vector).
%

    % --- Input Handling and Default FPS ---
    % Set the default FPS value if the second argument is not provided.
    if nargin < 2
        fps = 120;
    end

    % Input validation (optional, but good practice)
    if fps <= 0
        error('FPS must be a positive value.');
    end
    if any(frame_num < 0)
        error('Frame number cannot be negative.');
    end
    
    % Ensure input is treated as a column vector for consistent output structure
    if isrow(frame_num)
        frame_num = frame_num';
    end

    % --- Time Calculation (Vectorized Operations) ---

    % 1. Calculate the total time in seconds
    total_time_s = frame_num / fps;

    % 2. Calculate the whole minutes part
    minutes = floor(total_time_s / 60);

    % 3. Calculate the remaining time in seconds (including fractional part)
    remaining_s = total_time_s - (minutes * 60);

    % 4. Calculate the whole seconds part
    seconds = floor(remaining_s);

    % 5. Calculate the remaining time as milliseconds and round it
    % Subtract the whole seconds part to get the fractional second
    fractional_s = remaining_s - seconds;
    ms = round(fractional_s * 1000);

    % --- Vectorized Roll-over Logic ---

    % 6. Handle ms rolling over from 999 to 1000 (which becomes 0ms and +1s)
    rollover_ms = (ms >= 1000);
    ms(rollover_ms) = 0; % Reset ms to 0 for those frames
    seconds(rollover_ms) = seconds(rollover_ms) + 1; % Increment seconds

    % 7. Handle seconds rolling over from 59 to 60 (which becomes 0s and +1m)
    rollover_s = (seconds >= 60);
    seconds(rollover_s) = seconds(rollover_s) - 60; % Reset seconds by subtracting 60
    minutes(rollover_s) = minutes(rollover_s) + 1; % Increment minutes

end
