function GazeCalShepherd
%GazeCalibrationShepherd
%attempt save
subjectID = 'FLIC_200';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% STEP 1: make a perimeter file from raw data
saveFolders = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/'];
perimeterFile = [saveFolders, subjectID, '_gazeCal_perimeter.mat']; % path to perimeter file
perimeter = load(perimeterFile, 'perimeter');
perimeter = perimeter.perimeter;
% STEP 2: find the start frame from the playable pupil camera video using
%   IINA. Also calculate the duration of the dots from the playable world
%   camera video. sometimes the first dot is shorter than the rest.
startTime = [1, 25, 400]; % [minutes, seconds, milliseconds]
firstDotEnd = [1 26 608];
secondDotEnd = [1 30 000];
thirdDotEnd =[1 39 983];

fps = 120;
% calculate the duration of the first dot becuase it is often shorter
% becuase of eye opening.
firstDotDurFrames = time2frame(firstDotEnd) - time2frame(startTime);
firstDotDurS = firstDotDurFrames/fps;
% calculate the duration of the second and third dots and average them
secondDotDurFrames = time2frame(secondDotEnd) - time2frame(firstDotEnd);
secondDotDurS = secondDotDurFrames/fps;

thirdDotDurFrames = time2frame(thirdDotEnd) - time2frame(secondDotEnd);
thirdDotDurS = thirdDotDurFrames/fps;

targetDurSec = mean([thirdDotDurS, secondDotDurS]);
targetDurSec = 3.34;
%% STEP 3: find gaze frames to use in scene geometry estimation
% load run file for this participant
folders = ['/FLIC_data/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/'];
searchPattern = [dropboxBasedir, folders, '/', subjectID, '_gazeCal_runData.mat*'];
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

%determine frame numbers to analyze
confidenceCutoff = 0.8; % FLIC_2004

fullFrameSet = findGazeFrames(startTime, gazeTargetsDeg, perimeterFile, targetDurSec, firstDotDurS, confidenceCutoff);
goodIdx = find(~isnan(fullFrameSet));
fullFrameSet = fullFrameSet(goodIdx);
gazeTargetsDeg = gazeTargetsDeg(goodIdx,:);
%% STEP 4: estimate scene geometry
% first, use only 5 points and use the p output for your second round of
% searching

% Define the input variables for this particular gaze cal video
gazeSubsetIdx = [1,5:9]; % NEEDS TO BE ADJUSTED IF THERE ARE NANs

frameSet = fullFrameSet(gazeSubsetIdx);
gazeTargets = (gazeTargetsDeg(gazeSubsetIdx,:)).*[-1,1];

% use this to make sure the points look like they make a cross
plotPupilCenters(fullFrameSet, perimeter, [1:33]);

% Define some properties of the eye and of the scene that will be fixed
% for the scene search
sceneArgs = {
    'cameraGlintSourceRelative',[3.6;-3.2;6.9], ...
    'intrinsicCameraMatrix',[561.471804, 0.0, 200; 0.0, 562.494105, 200; 0.0, 0.0, 1.0],...
    'sensorResolution',[400 400],...
    'radialDistortionVector',[0 0]};

% Add the args for this particular observer
correctionType = 'spectacleLens';
%correctionType = 'contactLens'; 
% note contact lens only takes 2 key value pairs
%observerArgs = {'sphericalAmetropia',-1.25,'spectacleLens',[-1.25,0,0]};
%observerArgs = {'sphericalAmetropia',-1.25};
observerArgs = {'sphericalAmetropia',-1.00,correctionType,[-1.00,-0.25,0]};


% Combine the two argument sets
setupArgs = [sceneArgs observerArgs];
confidenceThreshold = confidenceCutoff;


% This is the x0, in case we want to pass that
% pick u here
x0 = [-28.6484   -7.3094   51.0564   24.3158    0.5042   12.1706    0.9918 0.9927   18.8754   49.3395   40.5355];

[sceneGeometry,p5] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', x0, 'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);

[sceneGeometry,p5] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5, 'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.

%% now again with the first half of the gaze targets
%fullFrameSet = sort([fullFrameSet; 13776; 21995]); % use if adding frames
%manually
%fullFrameSet = fullFrameSet(1:end-2);
frameSet = fullFrameSet(1:16);
gazeTargets = gazeTargetsDeg(1:16,:).*[-1,1];
[sceneGeometry,p17] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5, 'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);

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

plotPupilCenters(fullFrameSet, perimeter, [17:33]);

[sceneGeometry,p18] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p5, 'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);

% CHECK the graphs. Do the xs and os overlap well? Is the f value below 4?
%If no, investigate the playable video of the pupil camera and see if any
%of the points look poorly outlined. They may need to be omitted from the
%procedure.
%% now run with all gaze targets
close all
% use mean x0 from 1st and 2nd half as starting point
pMean = mean([p17(:), p18(:)],2);

frameSet = fullFrameSet;
gazeTargets = gazeTargetsDeg.*[-1,1];
[sceneGeometry,p34] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets,...
    'setupArgs', setupArgs, 'x0', pMean',...
    'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);
[sceneGeometry,p34] = estimateSceneGeometry(perimeterFile, fullFrameSet, gazeTargets,...
    'setupArgs', setupArgs, 'x0', p34,...
    'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);
% to redo quickly:
%[sceneGeometry,p34] = estimateSceneGeometry(perimeterFile, fullFrameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p34, 'confidenceThreshold', confidenceThreshold, 'nWorkers', 6);

%% Save things so we could regenerate the scene geometry file if needed!
% if everything looks good, save the scene geometry, p34, gaze offset, x0
% and frames used for this participant in a file!

% what is the gaze offset?? HUMAN
gazeOffset = [-3.6, -0.7]; % [azi, ele]

sceneGeometryFile = [saveFolders, subjectID, '_gazeCal_session-1_SceneGeometry.mat'];
saveFileMeta = [saveFolders, subjectID, '_gazeCal_session-1_SceneGeometryMetadata.mat'];
save(sceneGeometryFile, 'sceneGeometry')
save(saveFileMeta, "p34", "gazeOffset", "fullFrameSet", "gazeTargets", "startTime", "observerArgs", "confidenceThreshold");

% Save the figure 1 as a MATLAB figure file
fig1_handle = figure(35);
saveas(fig1_handle, [saveFolders, subjectID, '_gazeCal_SceneGeometryTargetsPlot'], 'fig');
%% How to turn pupil perimeters into gaze angles now that you have scene geometry
% Define variables for the path to the sceneGeometry file, perimeter file, and a _pupilData.mat file (which is to be created).
% Issue this command: fitPupilPerimeter(perimeterFileName, pupilFileName,'sceneGeometryFileName',sceneGeometryFileName,'useParallel',true,'verbose',true);
pupilFileName = [saveFolders, subjectID, '_gazeCal_pupilData.mat'];
fitPupilPerimeter(perimeterFile,pupilFileName, 'sceneGeometryFileName',sceneGeometryFile,'useParallel',true,'verbose',true, 'nWorkers', 6);
%% Clean up: remove low confidence points and something.
load(pupilFileName);
RMSECutoff = 2;

% Create indices for bad data with bad RMSE or fitAtBound
RMSEVals = pupilData.sceneConstrained.ellipses.RMSE;
badRMSEIdx = RMSEVals > RMSECutoff;
badFitAtBoundIdx = pupilData.sceneConstrained.eyePoses.fitAtBound > 0;
badIdx = badRMSEIdx | badFitAtBoundIdx; 
pupilData.sceneConstrained.eyePoses.values(badIdx, [1 2]) = NaN;
%% Smooth the pupil perimeters
aziLowerBound = -25 + gazeOffset(1,1);
aziUpperBound = 25 + gazeOffset(1,1);
eleLowerBound = -25 + gazeOffset(1,2);
eleUpperBound = 25 + gazeOffset(1,2);

[pupilData] = smoothPupilRadius(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'useParallel', true, 'nWorkers', 6,...
    'eyePoseLB', [aziLowerBound, eleLowerBound, 0, 0.5], 'eyePoseUB', [aziUpperBound, eleUpperBound, 0, 0.5]);

load([saveFolders, subjectID, '_gazeCal_pupilData.mat'])
figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,2), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,2), '.-')
ElevationFigHandle = figure(70);
saveas(ElevationFigHandle, [saveFolders, subjectID, '_gazeCal_eyePosEle_smoothed'], 'fig');

figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,1), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,1), '.-')
ElevationFigHandle = figure(71);
saveas(ElevationFigHandle, [saveFolders, subjectID, '_gazeCal_eyePosAzi_smoothed'], 'fig');
end
