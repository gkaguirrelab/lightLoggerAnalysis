subjectID = 'FLIC_2001';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% STEP 1: make a perimeter file from raw data
saveFolders = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/walkIndoor/temporalFrequency/'];
sceneGeometryFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometry.mat'];
perimeterFile = [saveFolders, subjectID, '_walkIndoor_tf_perimeter_contrast-1x5.mat']; % path to perimeter file
gazeCalMetaFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometryMetadata.mat'];
load(gazeCalMetaFile);

gazeOffset = [0.9512, 8.0238];
pupilFileName = [saveFolders, subjectID, '_walkIndoor_pupilData_contrast-1x5.mat'];
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
aziLowerBound = -40 + gazeOffset(1,1);
aziUpperBound = 40 + gazeOffset(1,1);
eleLowerBound = -60 + gazeOffset(1,2);
eleUpperBound = 40 + gazeOffset(1,2);

[pupilData] = smoothPupilRadius(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'useParallel', true, 'nWorkers', 6,...
    'eyePoseLB', [aziLowerBound, eleLowerBound, 0, 0.5], 'eyePoseUB', [aziUpperBound, eleUpperBound, 0, 0.5], 'exponentialTauParam', 20);

load([saveFolders, subjectID, '_gazeCal_pupilData.mat'])
figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,2), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,2), '.-')
ElevationFigHandle = figure(1);
saveas(ElevationFigHandle, [saveFolders, subjectID, '_gazeCal_eyePosEle_smoothed'], 'pdf');

figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,1), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,1), '.-')
ElevationFigHandle = figure(2);
saveas(ElevationFigHandle, [saveFolders, subjectID, '_gazeCal_eyePosAzi_smoothed'], 'pdf');