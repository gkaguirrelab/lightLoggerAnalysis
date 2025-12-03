subjectID = 'FLIC_2001';
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

% load files
saveFolders = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/lunch/temporalFrequency/'];
sceneGeometryFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometry.mat'];
perimeterFile = [saveFolders, subjectID, '_lunch_tf_perimeter_contrast1x25gamma1.mat']; % path to perimeter file
gazeCalMetaFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometryMetadata.mat'];
load(gazeCalMetaFile);

pupilFileName = [saveFolders, subjectID, '_lunch_pupilData_contrast1x25gamma1.mat'];
load(pupilFileName);

%% see how things would look after smoothing
[~,priorPupilRadius, finalGoodIdx] = smoothPupilTime(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'verbose', true, 'pupilVarThresh',0.1,...
    'ellipseRMSEThresh',2, 'confidenceThreshold',0.6,'confidenceThresholdDown',0.5,'currField', 'sceneConstrained', ...
    'ellipseRMSEThreshDown', 5, 'minPerimPointsDown', 3 , 'dryRun', true);

% plotting to check
nFrames = length(priorPupilRadius);
timeIndex = 1:nFrames;
goodIdxOrig = ~isnan(pupilData.smoothPupilTime_02.eyePoses.values(:,4));
figure; hold on;
plot(pupilData.sceneConstrained.eyePoses.values(:,4), 'DisplayName', 'Original Prior Radius');
plot(priorPupilRadius,'DisplayName', 'Smoothed Radius Constraint');
plot(timeIndex(finalGoodIdx), priorPupilRadius(finalGoodIdx), '.', 'DisplayName', 'Frames Kept by All Filters');
plot(pupilData.smoothPupilTime_02.eyePoses.values(:,4), '.', 'DisplayName', 'original smoothed');

xlabel('Frame');
ylabel('Pupil Radius (mm)');
title('Pupil Radius Time Series and Filtering');
legend('show', 'Location', 'best');
hold off;