function frame = estimateFrameFromTime(fps, time)
% given frame rate (fps) and a timepoint in the video, figure out frame
% number
%Inputs
%   fps (Hz)
%   time (minute second millisecond) The human observer should determine
%   this using IINA
%Example 
%{
    fps = 120;
    time = [1, 23, 683]
    frame = estimateFrameFromTime(fps, time);

%}
timeSec = time(:,1)*60 + time(:,2) +time(:,3)/1000; % convert to seconds
frame = timeSec.*fps;
end