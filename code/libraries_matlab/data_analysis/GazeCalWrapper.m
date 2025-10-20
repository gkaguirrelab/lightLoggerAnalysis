%GazeCalibrationWrapper

% STEP 1: make a perimeter file from raw data
perimeterFile = []; % path to perimeter file

% STEP 2: find the start frame from the playable pupil camera video using
%   IINA. Also calculate the duration of the dots from the playable world
%   camera video. sometimes the first dot is shorter than the rest.

%% STEP 3: find gaze frames to use in scene geometry estimation
startTime = [1, 23, 683]; % [minutes, seconds, milliseconds]
targetDurSec = 3.267;
gazeTargetsDeg = [ ...
    0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
    0, 15; 0, -15; -15, 0; 15, 0;...
    -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
    0, 10; 0, -7.5; -7.5, 0; 7.5, 0;...
    0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
    0, 15; 0, -15; -15, 0; 15, 0;...
    -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
    0, 10; 0, -7.5; -7.5, 0; 7.5, 0];

fullFrameSet = findGazeFrames(startTime, gazeTargetsDeg, targetDurSec);

%% STEP 4: estimate scene geometry
% first, use only 5 points and use the p output for your second round of
% searching

% Define the input variables for this particular gaze cal video
gazeSubsetIdx = [1,6,7,8,9];
frameSet = fullFrameSet(gazeSubsetIdx);
gazeTargets = gazeTargetsDeg(:,gazeTargetsDeg);

% Define some properties of the eye and of the scene that will be fixed
% for the scene search
sceneArgs = {
    'cameraGlintSourceRelative',[3.6;-3.2;6.9], ...
    'intrinsicCameraMatrix',[561.471804, 0.0, 200; 0.0, 562.494105, 200; 0.0, 0.0, 1.0],...
    'sensorResolution',[400 400],...
    'radialDistortionVector',[0 0]};

% Add the args for this particular observer
observerArgs = {'sphericalAmetropia',2.25,'spectacleLens',[2.25,2,80]};

% Combine the two argument sets
setupArgs = [sceneArgs observerArgs];


% This is the x0, in case we want to pass that
% THIS NEEDS TO BE UPDATED 
x0 = [-26.8458  -14.9894   48.2704   25.2714   -1.9474    9.5655    0.9503    1.0497   13.9882   43.9961   44.9963];

% Run the routine
[sceneGeometry,p] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs);

% Search again starting from the prior search result
[sceneGeometry,p] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p);

% now again with all the gaze targets
frameSet = fullFrameSet;
gazeTargets = gazeTargetsDeg;
[sceneGeometry,p] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p);
