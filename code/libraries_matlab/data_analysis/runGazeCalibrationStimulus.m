function runGazeCalibrationStimulus
% Displays 13-dot gaze calibration stimulus at fixed visual angles,
% with a brief beep signaling each dot onset.

% Hard-coded parameters
viewingDistCm = 40;
dotRadiusDeg = 0.6;
dotTime = 1; 
repetitions = 1;
bgColor = [0 0 0];
fgColor = [255 255 255];
redColor  = [255   0   0];
innerFrac = 0.2;

AssertOpenGL;
screenNum = max(Screen('Screens'));

% Initialize PsychPortAudio for beep
InitializePsychSound(1);
sampleRate = 44100;
beepFreq = 1000;
beepLengthSec = 0.1;
nrchannels = 1; 
% Generate beep
beep = MakeBeep(beepFreq, beepLengthSec, sampleRate);
% Open audio device
pahandle = PsychPortAudio('Open', [], 1, 1, sampleRate, nrchannels);
% Fill buffer with beep waveform
PsychPortAudio('FillBuffer', pahandle, beep);

% Query physical display size (mm) and convert to cm
[widthMm, heightMm] = Screen('DisplaySize', screenNum);
widthCm  = widthMm / 10;
heightCm = heightMm / 10;

% Open window and get its pixel resolution
[win, winRect] = Screen('OpenWindow', screenNum, bgColor);
[xCen, yCen] = RectCenter(winRect);
winWidthPx = winRect(3) - winRect(1);
winHeightPx = winRect(4) - winRect(2);

% Print physical size and pixels to confirm display
fprintf('→ Screen %d: %.1f×%.1f cm physical → %d×%d px window\n', ...
        screenNum, widthCm, heightCm, winWidthPx, winHeightPx);

% Compute pixels-per-cm from actual window
pxPerCmX = winWidthPx  / widthCm;
pxPerCmY = winHeightPx / heightCm;

% Degree→pixel conversion lambdas
deg2pxX = @(deg) viewingDistCm * tand(deg) * pxPerCmX;
deg2pxY = @(deg) viewingDistCm * tand(deg) * pxPerCmY;
dotRadiusPx = viewingDistCm * tand(dotRadiusDeg) * pxPerCmX;

% Define 13 calibration positions in degrees [xDeg, yDeg]
degPositions = [ ...
      0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
      0, 20;   0, -20;   -20,   0;   20, 0; ...
     -15, 15;  15, 15;   -15, -15;   15, -15; ...

      -10,  10;   -10, -10;    10,  10;    10, -10; ...
      0,  10;    0, -10;   -10,   0;     10,   0; ...
      -5,   5; 5,   5;   -5,  -5;     5,  -5;     0,   0];
nDots = size(degPositions, 1);

% Check CM distance
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

% Instruction text
text = [
    'You will see 13 calibration dots appear one by one on the screen,' newline ...
    'each signaled by a brief beep. Please fix your gaze on each dot when it appears,' newline ...
    'and do not try to anticipate the location of the next dot.' newline newline ...
    'Press any key to begin.'
];
Screen('TextSize', win, 24);
DrawFormattedText(win, text, 'center', 'center', fgColor);

% Unify key names across platforms
KbName('UnifyKeyNames');
KbReleaseWait;
FlushEvents('keyDown');
Screen('Flip', win);

% Wait for key press
KbWait(-1);

% Draw loop with beep on each dot onset
HideCursor;
Priority(MaxPriority(win));
for rep = 1:repetitions
    for i = 1:nDots
        % Play beep (single repetition)
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        
        % Draw dot
        pos  = positions(i,:);
        % Outer (white) circle
        outerRect = [pos(1)-dotRadiusPx, pos(2)-dotRadiusPx, ...
                     pos(1)+dotRadiusPx, pos(2)+dotRadiusPx];
        Screen('FillOval', win, fgColor, outerRect);
        
        % Inner (red) circle
        innerRadiusPx = dotRadiusPx * innerFrac;
        innerRect = [pos(1)-innerRadiusPx, pos(2)-innerRadiusPx, ...
                     pos(1)+innerRadiusPx, pos(2)+innerRadiusPx];
        Screen('FillOval', win, redColor, innerRect);
        WaitSecs(dotTime);

        Screen('Flip', win);
        WaitSecs(dotTime);
    end
    Screen('Flip', win);
    WaitSecs(dotTime);
end
Priority(0);
ShowCursor;
Screen('CloseAll');
end