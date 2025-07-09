function runGazeCalibrationStimulus()
% Displays 13-dot gaze calibration stimulus at fixed visual angles,
% with a brief beep signaling each dot onset.

% Hard-coded parameters
viewingDistCm = 50;       % viewing distance (cm)
dotRadiusDeg = 0.2;       % dot radius (degrees)
dotTime      = 1;         % display time per dot (s)
repetitions  = 1;         % cycles through all dots
bgColor      = [0 0 0];
fgColor      = [255 255 255];

AssertOpenGL;
screenNum = max(Screen('Screens'));

% Initialize PsychPortAudio for beep
InitializePsychSound(1);
sampleRate    = 44100;    % audio sample rate (Hz)
beepFreq      = 1000;     % beep frequency (Hz)
beepLengthSec = 0.1;      % beep duration (s)
nrchannels    = 1;        % mono beep
% Generate beep
beep = MakeBeep(beepFreq, beepLengthSec, sampleRate);
% Open audio device
pahandle = PsychPortAudio('Open', [], 1, 1, sampleRate, nrchannels);
% Fill buffer with beep waveform
PsychPortAudio('FillBuffer', pahandle, beep);

% Query physical display size (mm) → convert to cm
[widthMm, heightMm] = Screen('DisplaySize', screenNum);
widthCm  = widthMm  / 10;
heightCm = heightMm / 10;

% Open window and get its pixel resolution
[win, winRect] = Screen('OpenWindow', screenNum, bgColor);
[xCen, yCen]  = RectCenter(winRect);
winWidthPx     = winRect(3) - winRect(1);
winHeightPx    = winRect(4) - winRect(2);

% Compute pixels-per-cm from actual window
pxPerCmX = winWidthPx  / widthCm;
pxPerCmY = winHeightPx / heightCm;

% Degree→pixel conversion lambdas
deg2pxX     = @(deg) viewingDistCm * tand(deg) * pxPerCmX;
deg2pxY     = @(deg) viewingDistCm * tand(deg) * pxPerCmY;
dotRadiusPx = viewingDistCm * tand(dotRadiusDeg) * pxPerCmX;

% Define 13 calibration positions in degrees [xDeg, yDeg]
degPositions = [ ...
      0,   0;   -10,  10;   -10, -10;    10,  10;    10, -10; ...
      0,  10;    0, -10;   -10,   0;     10,   0;   -5,   5; ...
      5,   5;   -5,  -5;     5,  -5;     0,   0];
nDots = size(degPositions, 1);

% CHECK CM DISTANCE
xCm = viewingDistCm * tand(degPositions(:,1));
yCm = viewingDistCm * tand(degPositions(:,2));
rCm = hypot(xCm, yCm);

fprintf(' Dot |  x°   |  y°   |   x_cm   |   y_cm   |  r_cm\n');
fprintf('-----+-------+-------+----------+----------+--------\n');
for i = 1:nDots
    fprintf('%4d | %5.1f° | %5.1f° | %7.2f  | %7.2f  | %6.2f\n', ...
        i, degPositions(i,1), degPositions(i,2), xCm(i), yCm(i), rCm(i));
end
fprintf('\n');

% Convert to pixel coordinates
positions = zeros(nDots, 2);
for i = 1:nDots
    xOff = deg2pxX(degPositions(i,1));
    yOff = deg2pxY(degPositions(i,2));
    positions(i,:) = [xCen + xOff, yCen - yOff];
end

% Draw loop with beep on each dot onset
HideCursor;
Priority(MaxPriority(win));
for rep = 1:repetitions
    for i = 1:nDots
        % Play beep (single repetition)
        PsychPortAudio('Start', pahandle, 1, 0, 0);

        % Draw dot
        pos  = positions(i,:);
        rect = [pos(1)-dotRadiusPx, pos(2)-dotRadiusPx, ...
                pos(1)+dotRadiusPx, pos(2)+dotRadiusPx];
        Screen('FillOval', win, fgColor, rect);
        Screen('Flip', win);
        WaitSecs(dotTime);
    end
end
Priority(0);
ShowCursor;
Screen('CloseAll');
% Close audio device
tPsychPortAudio('Close', pahandle);
end