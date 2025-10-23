%GazeCalibrationShepherd
subjectID = 'FLIC_2002';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% STEP 1: make a perimeter file from raw data
perimeterFile = [dropboxBasedir, '/FLIC_data/lightLogger/Processing/', subjectID '_gazeCalibration_session1_perimeter.mat']; % path to perimeter file

% STEP 2: find the start frame from the playable pupil camera video using
%   IINA. Also calculate the duration of the dots from the playable world
%   camera video. sometimes the first dot is shorter than the rest.

%% STEP 3: find gaze frames to use in scene geometry estimation
% load run file for this participant
folders = ['/FLIC_data/lightLogger/GazeCalRunFileData/', subjectID];
searchPattern = [dropboxBasedir, folders, '/', subjectID, '_GazeCalibration_session1*'];
fileList = dir(searchPattern);
if isempty(fileList)
    error('No files found matching the pattern: %s', searchPattern);
end
fileName = fileList(1).name;
runDataPath = fullfile(fileList(1).folder, fileName);

runData = load(runDataPath, 'taskData');
runData = runData.taskData;
% pull out gaze target positions from this file
gazeTargetsDeg = runData.gaze_target_positions_deg;

%check timing inputs (human!)
startTime = [1, 20, 933]; % [minutes, seconds, milliseconds]
targetDurSec = 3.267;
%determine frame numbers to analyze
fullFrameSet = findGazeFrames(startTime, gazeTargetsDeg, perimeterFile, targetDurSec);

%% STEP 4: estimate scene geometry
% first, use only 5 points and use the p output for your second round of
% searching

% Define the input variables for this particular gaze cal video
gazeSubsetIdx = [1,6,7,8,9];
frameSet = fullFrameSet(gazeSubsetIdx);
gazeTargets = (gazeTargetsDeg(gazeSubsetIdx,:)).*[-1,1];

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
x0 = [-29.9355  -10.3699   52.3664   24.2541    2.1374   15.5410    0.9895    1.0024   17.6172   43.0417   41.0665];
% Run the routine
[sceneGeometry,p] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs);

% Search again starting from the prior search result
[sceneGeometry,p5] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.

%% now again with the first half of the gaze targets
frameSet = fullFrameSet(1:17);
gazeTargets = gazeTargetsDeg(1:17,:).*[-1,1];
[sceneGeometry,p17] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.
%% now again with the second half of the gaze targets
frameSet = fullFrameSet(18:end);
gazeTargets = gazeTargetsDeg(18:end,:).*[-1,1];
[sceneGeometry,p18] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.
%% now run with all gaze targets
% use mean x0 from 1st and 2nd half as starting point
pMean = mean([p17(:), p18(:)],2);

frameSet = fullFrameSet;
gazeTargets = gazeTargetsDeg.*[-1,1];
[sceneGeometry,p34] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', pMean);

%% Save things so we could regenerate the scene geometry file if needed!
% if everything looks good, save the scene geometry, p34, gaze offset, x0
% and frames used for this participant in a file!

