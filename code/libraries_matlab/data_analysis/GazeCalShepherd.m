function GazeCalShepherd
%GazeCalibrationShepherd
subjectID = 'FLIC_2003';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% STEP 1: make a perimeter file from raw data
perimeterFile = [dropboxBasedir, '/FLIC_data/lightLogger/scriptedIndoorOutdoor/', subjectID, '/', subjectID, '_gazeCalibration_session1_perimeter.mat']; % path to perimeter file
perimeter = load(perimeterFile, 'perimeter');
perimeter = perimeter.perimeter;
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
gazeTargetsDeg = vertcat(runData.gaze_target_positions_deg,runData.gaze_target_positions_deg);

%check timing inputs (human!)
startTime = [1, 20, 108]; % [minutes, seconds, milliseconds]
targetDurSec = 3.5;
onset_delay_s = 0.8; % again, human should calculate this based on the difference betwen start frame and first eye movement. What is that duration compared to the intended?
%determine frame numbers to analyze
fullFrameSet = findGazeFrames(startTime, gazeTargetsDeg, perimeterFile, targetDurSec, onset_delay_s);
goodIdx = find(~isnan(fullFrameSet));
fullFrameSet = fullFrameSet(goodIdx);
gazeTargetsDeg = gazeTargetsDeg(goodIdx,:);
%% STEP 4: estimate scene geometry
% first, use only 5 points and use the p output for your second round of
% searching

% Define the input variables for this particular gaze cal video
gazeSubsetIdx = [1,6,7,8,9]; % NEEDS TO BE ADJUSTED IF THERE ARE NANs
gazeSubsetIdx = [1,6:9]; % NEEDS TO BE ADJUSTED IF THERE ARE NANs

frameSet = fullFrameSet(gazeSubsetIdx);
gazeTargets = (gazeTargetsDeg(gazeSubsetIdx,:)).*[-1,1];

% use this to make sure the points look like they make a cross
figure; for ii=1:length(frameSet); Xp = perimeter.data{frameSet(ii)}.Xp; Yp = perimeter.data{frameSet(ii)}.Yp; plot(Xp,Yp,'x'); hold on; pause; end

% Define some properties of the eye and of the scene that will be fixed
% for the scene search
sceneArgs = {
    'cameraGlintSourceRelative',[3.6;-3.2;6.9], ...
    'intrinsicCameraMatrix',[561.471804, 0.0, 200; 0.0, 562.494105, 200; 0.0, 0.0, 1.0],...
    'sensorResolution',[400 400],...
    'radialDistortionVector',[0 0]};

% Add the args for this particular observer
%observerArgs = {'sphericalAmetropia',-1.25,'spectacleLens',[-1.25,0,0]};
observerArgs = {'sphericalAmetropia',-1.25};


% Combine the two argument sets
setupArgs = [sceneArgs observerArgs];


% This is the x0, in case we want to pass that
x0 = [-28.6484   -7.3094   51.0564   24.3158    0.5042   12.1706    0.9918 0.9927   18.8754   49.3395   40.5355];

[sceneGeometry,p5] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', x0);

[sceneGeometry,p5] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.

%% now again with the first half of the gaze targets
fullFrameSet = sort([fullFrameSet; 13776; 21995]);
fullFrameSet = fullFrameSet(1:end-2);
frameSet = fullFrameSet(1:16);
gazeTargets = gazeTargetsDeg(1:16,:).*[-1,1];
[sceneGeometry,p17] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.
%% now again with the second half of the gaze targets
subset = [17:33];
frameSet = fullFrameSet(subset);
gazeTargets = gazeTargetsDeg(subset,:).*[-1,1];

figure; for ii=14:16 %length(frameSet); 
Xp = perimeter.data{frameSet(ii)}.Xp; Yp = perimeter.data{frameSet(ii)}.Yp; plot(Xp,Yp,'x'); hold on; pause; end

plotCenters(fullFrameSet, perimeter, [17:33]);

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
[sceneGeometry,p34] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', pMean');

%% Save things so we could regenerate the scene geometry file if needed!
% if everything looks good, save the scene geometry, p34, gaze offset, x0
% and frames used for this participant in a file!

% what is the gaze offset?? HUMAN
gazeOffset = []; % [azi, ele]

saveFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/', subjectID, 'SceneGeometry.mat'];
saveFileMeta = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/', subjectID, 'SceneGeometryMetadata.mat'];
save(saveFile, 'sceneGeometry')
save(saveFileMeta, "p34", "gazeOffset", "fullFrameSet", "gazeTargets")
end
