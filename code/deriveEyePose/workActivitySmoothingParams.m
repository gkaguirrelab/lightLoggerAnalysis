%Smoothing parameters

%2001 Lunch
[pupilData,priorPupilRadius] = smoothPupilTime(perimeterFile, pupilFileName,...
    sceneGeometryFile, 'verbose', true, 'justReturnPrior', false, 'pupilVarThresh',0.8,...
    'ellipseRMSEThresh',2, 'confidenceThreshold',0.7, 'nWorkers', 8);


% 2003 work
using ellipseRMSEThresh = 1.5, confidenceThreshold = 0.7

% 2004 Work
eyePoseLB = [-45,-45,0,1.5];
eyePoseUB = [45,45,0,3];
smoothPupilTime(perimeterFileName, pupilFileName, sceneGeometryFile,'ellipseRMSEThresh',4,'confidenceThreshold',0.7,'eyePoseUB',eyePoseUB,'eyePoseLB',eyePoseLB);

%2006
[pupilData, priorPupilRadius] = smoothPupilTime(perimeterFileName, pupilFileName, sceneGeometryFile,'ellipseRMSEThresh',4.5,'justReturnPrior',false,'confidenceThreshold',0.7,'nWorkers',10,'verbose',true,'useParallel',true);