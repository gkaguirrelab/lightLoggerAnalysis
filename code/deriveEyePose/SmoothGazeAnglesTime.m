subjectID = {'FLIC_2001', 'FLIC_2003', 'FLIC_2004', 'FLIC_2005', 'FLIC_2006'};
dropboxBasedir = fullfile(getpref("lightLoggerAnalysis", 'dropboxBaseDir'));

%activity = 'walkIndoor';
activity = 'work'; % 'lunch'
% load files
for subjIdx = 1:length(subjectID)
    saveFolders = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID{subjIdx}, '/', activity, '/temporalFrequency/'];
    sceneGeometryFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID{subjIdx}, '/gazeCalibration/temporalFrequency/', subjectID{subjIdx}, '_gazeCal_SceneGeometry.mat'];
    %searchPattern = [saveFolders, subjectID{subjIdx}, '_', activity, '_tf_perimeter_contrast*'];
    searchPattern = [saveFolders, subjectID{subjIdx}, '_', activity, '_tf_perimeter*'];
    fileList = dir(searchPattern);
    if isempty(fileList)
        error('No files found matching the pattern: %s', searchPattern);
    end
    fileName = fileList(1).name;
    perimeterFile = [saveFolders, fileName];
    % gazeCalMetaFile = [dropboxBasedir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/', subjectID, '_gazeCal_SceneGeometryMetadata.mat'];
    % load(gazeCalMetaFile);

    pupilFileName = [saveFolders, subjectID{subjIdx}, '_', activity, '_pupilData_contrast1x25gamma1.mat'];
    if ~isfile(pupilFileName)
        confidenceThreshold = 0.7;
        fitPupilPerimeter(perimeterFile,pupilFileName, 'sceneGeometryFileName',sceneGeometryFile,'useParallel',true,'verbose',true, 'nWorkers', 5, 'confidenceThreshold', confidenceThreshold);
    end
    update = [subjectID{subjIdx}, ' pupilData complete'];
    disp(update);
end

%start here for smoothing.
load(pupilFileName);

%% see how things would look after smoothing
[~,priorPupilRadius] = smoothPupilTime(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'verbose', true, 'justReturnPrior', true, 'pupilVarThresh',0.8,...
    'ellipseRMSEThresh',2, 'confidenceThreshold',0.7);

nFrames = length(priorPupilRadius);
badIdx = pupilData.sceneConstrained.ellipses.RMSE > 2;
nBad = sum(badIdx);
goodIdx = find(badIdx ==0);
perceptBad = nBad/nFrames*100
% plotting to check

timeIndex = 1:nFrames;
%goodIdxOrig = ~isnan(pupilData.smoothPupilTime_02.eyePoses.values(:,4));
figure; hold on;
plot(pupilData.sceneConstrained.eyePoses.values(:,4), 'DisplayName', 'Original Prior Radius');
plot(priorPupilRadius,'DisplayName', 'Smoothed Radius Constraint');
plot(timeIndex(goodIdx), priorPupilRadius(goodIdx), '.', 'DisplayName', 'Good RMSE');

xlabel('Frame');
ylabel('Pupil Radius (mm)');
title('Pupil Radius Time Series and Filtering');
legend('show', 'Location', 'best');
hold off;

%% If that looks good, do the smoothing!
[pupilData,priorPupilRadius] = smoothPupilTime(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'verbose', true, 'justReturnPrior', false, 'pupilVarThresh',0.8,...
    'ellipseRMSEThresh',2, 'confidenceThreshold',0.7, 'nWorkers', 8);



plot(pupilData.smoothPupilTime_02.eyePoses.values(:,4), '.', 'DisplayName', 'original smoothed');
