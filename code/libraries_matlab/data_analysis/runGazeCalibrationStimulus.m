function runGazeCalibrationStimulus(monitorDiag, aspectRatio)
% Displays 13-dot gaze calibration stimulus at fixed visual angles.
% 
% Inputs:
%   monitorDiag         - Scalar. Monitor diagonal in inches.
%   aspectRatio         - Scalar. width/height ratio.

% Harcoded parameters. Subject to change.
viewingDistCm = 50;        % viewing distance from participant (cm)
dotSpacingDeg = 10;       % degrees belibtween dots from center
dotRadiusDeg = 0.1;       % radius of dot in degrees
dotTime = 1;              % time each dot is displayed (seconds)
repetitions = 1;          % number of times to cycle through all dots

backgroundColor = [0 0 0];
foregroundColor = [255 255 255];

% Get screen specs
AssertOpenGL;
screenNumber = max(Screen('Screens'));
res = Screen('Resolution', screenNumber);
monitorRes = res.width;

% Compute monitor width from diagonal and aspect ratio
monitorDiagCm = monitorDiag * 2.54;
monitorWidthCm = monitorDiagCm * (aspectRatio / sqrt(1 + aspectRatio^2));
pixelsPerCm = monitorRes / monitorWidthCm;

% Compute dot spacing and radius in pixels
dotSpacingPx = viewingDistCm * tand(dotSpacingDeg) * pixelsPerCm;
dotRadiusPx = 2 * viewingDistCm * tand(dotRadiusDeg / 2) * pixelsPerCm;

% Open window
[win, winRect] = Screen('OpenWindow', screenNumber, backgroundColor);
[xCen, yCen] = RectCenter(winRect);

% Generate 13 target positions (9 in 3x3 + 4 diagonals)
positions = [];
offsets = [-dotSpacingPx, 0, dotSpacingPx];
for x = offsets
    for y = offsets
        positions = [positions; [xCen + x, yCen + y]];
    end
end
% Generate 4 diagonal positions (halfway to corners)
diagOffset = dotSpacingPx / sqrt(2);
diagCoords = [ diagOffset,  diagOffset;
               diagOffset, -diagOffset;
              -diagOffset,  diagOffset;
              -diagOffset, -diagOffset];
for i = 1:4
    positions = [positions; [xCen + diagCoords(i,1), yCen + diagCoords(i,2)]];
end

% Play video
HideCursor;
Priority(MaxPriority(win));

for rep = 1:repetitions
    for i = 1:size(positions,1)
        % Draw dot
        pos = positions(i,:);
        Screen('FillOval', win, foregroundColor, [pos(1)-dotRadiusPx, pos(2)-dotRadiusPx, pos(1)+dotRadiusPx, pos(2)+dotRadiusPx]);
        Screen('Flip', win);
        WaitSecs(dotTime);
    end
end

Priority(0);
ShowCursor;
Screen('CloseAll');

end
