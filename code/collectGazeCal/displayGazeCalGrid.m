function degPositions = displayGazeCalGrid(heightCm, widthCm, viewingDistCm)
% DISPLAYGAZECALGRID
% Displays all gaze calibration dots simultaneously until a button is pressed.
%
% Inputs:
%   heightCm       - Screen height in cm
%   widthCm        - Screen width in cm
%   viewingDistCm  - Viewing distance in cm
%
% Output:
%   degPositions   - Nx2 array of dot positions in degrees visual angle
% Example
%{
    displayGazeCalGrid(106.7, 192.4, 100)
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