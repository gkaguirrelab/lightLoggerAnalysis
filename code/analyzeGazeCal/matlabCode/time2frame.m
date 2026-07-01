function frame = time2frame(time, fps)
% Convert a video timestamp to a frame number
%
% Syntax:
%   frame = time2frame(time)
%   frame = time2frame(time, fps)
%
% Description:
%   Given a timestamp in [minutes, seconds, milliseconds] format and a
%   frame rate, computes the corresponding frame number. Supports
%   vectorized input where each row of time is a separate timestamp.
%   Defaults to 120 fps if no frame rate is specified.
%
% Inputs:
%   time                  - Nx3 double. Timestamp(s) as [minutes, seconds,
%                           milliseconds]. Each row is one timestamp.
%   fps                   - Scalar double (default: 120). Frame rate in Hz.
%
% Outputs:
%   frame                 - Nx1 double. Corresponding frame number(s).
%
% Examples:
%{
    fps = 120;
    time = [1, 23, 683];
    frame = time2frame(time, fps);
%}
if nargin < 2
    fps = 120;
end
timeSec = time(:,1)*60 + time(:,2) +time(:,3)/1000; % convert to seconds
frame = timeSec.*fps;
end