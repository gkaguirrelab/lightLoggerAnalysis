function [output] = frame2time(frame_num, fps)
% Convert frame numbers to timestamps in [minutes, seconds, milliseconds]
%
% Syntax:
%   output = frame2time(frame_num)
%   output = frame2time(frame_num, fps)
%
% Description:
%   Converts one or more frame numbers into human-readable timestamps
%   expressed as [minutes, seconds, milliseconds]. Handles rollover of
%   milliseconds to seconds and seconds to minutes. Defaults to 120 fps
%   if no frame rate is specified.
%
% Inputs:
%   frame_num             - Numeric scalar or vector. The sequential frame
%                           number(s) to convert. Must be non-negative.
%   fps                   - Scalar double (default: 120). Frames per second
%                           rate. Must be positive.
%
% Outputs:
%   output                - Nx3 double. Matrix with columns [minutes,
%                           seconds, milliseconds] where N is the number
%                           of input frame numbers.
%
% Examples:
%{
    output = frame2time(14400);       % 2 minutes at 120 fps
    output = frame2time([120; 240], 120);
%}

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

    output = [minutes seconds ms]

end
