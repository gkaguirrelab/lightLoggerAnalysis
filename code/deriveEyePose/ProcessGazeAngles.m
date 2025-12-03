subjectID = 'FLIC_2005';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% load files
saveFolders = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/lunch/temporalFrequency/'];
sceneGeometryFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometry.mat'];
perimeterFile = [saveFolders, subjectID, '_lunch_tf_perimeter_contrast1x25gamma1.mat']; % path to perimeter file
gazeCalMetaFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometryMetadata.mat'];
load(gazeCalMetaFile);

%gazeOffset = [0.9512, 8.0238];
pupilFileName = [saveFolders, subjectID, '_lunch_pupilData_contrast1x25gamma1.mat'];
% FLIC_2002 crashing. try commenting out parfor and switching to for loop
% in pupilFitting code. turn on warnings to see what is going wrong
confidenceThreshold = 0.80;
fitPupilPerimeter(perimeterFile,pupilFileName, 'sceneGeometryFileName',sceneGeometryFile,'useParallel',true,'verbose',true, 'nWorkers', 5, 'confidenceThreshold', confidenceThreshold);
%% Clean up: remove low confidence points and something.
load(pupilFileName);
RMSECutoff = 2;

% Create indices for bad data with bad RMSE or fitAtBound
RMSEVals = pupilData.sceneConstrained.ellipses.RMSE;
badRMSEIdx = RMSEVals > RMSECutoff;
badFitAtBoundIdx = pupilData.sceneConstrained.eyePoses.fitAtBound > 0;
badIdx = badRMSEIdx | badFitAtBoundIdx;
percentBad = (sum(badIdx)/length(RMSEVals))*100
pupilData.sceneConstrained.eyePoses.values(badIdx, [1 2]) = NaN;
%% Smooth the pupil perimeters
figure;
plot(pupilData.sceneConstrained.eyePoses.values(:,2), '.-')
figure;
plot(pupilData.sceneConstrained.eyePoses.values(:,1), '.-')


aziLowerBound = -40 + gazeOffset(1,1);
aziUpperBound = 40 + gazeOffset(1,1);
eleLowerBound = -60 + gazeOffset(1,2);
eleUpperBound = 40 + gazeOffset(1,2);

[pupilData] = smoothPupilRadius(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'useParallel', true, 'nWorkers', 6,...
    'eyePoseLB', [aziLowerBound, eleLowerBound, 0, 0.5], 'eyePoseUB', [aziUpperBound, eleUpperBound, 0, 0.5], 'exponentialTauParam', 20);

figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,2), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,2), '.-')

figure; hold on
plot(pupilData.sceneConstrained.eyePoses.values(:,1), '.-')
plot(pupilData.radiusSmoothed.eyePoses.values(:,1), '.-')


startAprilTag = 8681;
endApriltag = 10375;
figure; hold on
plot(pupilData.radiusSmoothed.eyePoses.values(startAprilTag:endApriltag,1),...
    pupilData.radiusSmoothed.eyePoses.values(startAprilTag:endApriltag,2), 'o-')
ylabel('degrees elevation');
xlabel('degrees azimuth');
