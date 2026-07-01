function degPositions = displayGazeCalGrid(heightCm, widthCm, viewingDistCm)
% Display the full gaze-calibration target grid on screen.
%
% Syntax:
%   degPositions = displayGazeCalGrid(heightCm, widthCm, viewingDistCm)
%
% Description:
%   This function opens a Psychtoolbox window, converts a fixed set of
%   target coordinates from degrees of visual angle into screen pixels,
%   and draws all calibration targets simultaneously. The routine is meant
%   as a quick visualization or alignment aid rather than the timed
%   calibration task itself. It waits for a key press before closing the
%   display and returns the degree-space target locations that were shown.
%
% Inputs:
%   heightCm                 - Scalar. Physical height of the display in
%                              centimeters.
%   widthCm                  - Scalar. Physical width of the display in
%                              centimeters.
%   viewingDistCm            - Scalar. Viewing distance from the observer
%                              to the display in centimeters.
%
% Outputs:
%   degPositions             - Numeric matrix. One [xDeg, yDeg] pair per
%                              calibration target, expressed in degrees of
%                              visual angle relative to screen center.
%
% Examples:
%{
    degPositions = displayGazeCalGrid(106.7, 192.4, 100);
%}

arguments
    heightCm = 106.7; % 2nd floor LGTV
    widthCm = 192.4; % 2nd floor LGTV
    viewingDistCm = 100; % 1 meter standard
end

AssertOpenGL;
KbName('UnifyKeyNames');

%% ---------------- Parameters ----------------
bgColor         = [0 0 0];
outerDotColor   = [255 255 255];
innerDotColor   = [255 0 0];

outerDotRadiusDeg = 1.5;
innerDotRadiusDeg = 0.25;

%% ---------------- Dot layout (same as your code) ----------------
degPositions = [ ...
     0,    0;
   -15,   15;
   -15,  -15;
    15,   15;
    15,  -15;
     0,   15;
     0,  -15;
   -15,    0;
    15,    0;
   -7.5,  7.5;
   -7.5, -7.5;
    7.5,  7.5;
    7.5, -7.5;
     0,   7.5;
     0,  -7.5;
   -7.5,  0;
    7.5,  0];

nDots = size(degPositions,1);

%% ---------------- Open screen ----------------
screenNum = max(Screen('Screens'));
[win, winRect] = Screen('OpenWindow', screenNum, bgColor);
[xCen, yCen] = RectCenter(winRect);

winWidthPx  = winRect(3);
winHeightPx = winRect(4);

pxPerCmX = winWidthPx  / widthCm;
pxPerCmY = winHeightPx / heightCm;

deg2pxX = @(deg) viewingDistCm * tand(deg) * pxPerCmX;
deg2pxY = @(deg) viewingDistCm * tand(deg) * pxPerCmY;

outerRadiusPx = viewingDistCm * tand(outerDotRadiusDeg) * pxPerCmX;
innerRadiusPx = viewingDistCm * tand(innerDotRadiusDeg) * pxPerCmX;

%% ---------------- Compute dot rectangles ----------------
outerRects = zeros(4, nDots);
innerRects = zeros(4, nDots);

for i = 1:nDots
    xPx = xCen + deg2pxX(degPositions(i,1));
    yPx = yCen - deg2pxY(degPositions(i,2));

    outerRects(:,i) = [ ...
        xPx - outerRadiusPx;
        yPx - outerRadiusPx;
        xPx + outerRadiusPx;
        yPx + outerRadiusPx ];

    innerRects(:,i) = [ ...
        xPx - innerRadiusPx;
        yPx - innerRadiusPx;
        xPx + innerRadiusPx;
        yPx + innerRadiusPx ];
end

%% ---------------- Draw all dots ----------------
Screen('FillRect', win, bgColor);
Screen('FillOval', win, outerDotColor, outerRects);
Screen('FillOval', win, innerDotColor, innerRects);

Screen('Flip', win);

%% ---------------- Wait for input ----------------
disp('Press any key to quit');
KbReleaseWait;  % ensure no key is held down
KbWait(-1);

%% ---------------- Cleanup ----------------
sca;
ShowCursor;

end
