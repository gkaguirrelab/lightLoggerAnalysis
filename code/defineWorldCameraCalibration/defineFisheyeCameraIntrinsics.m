% This minimal script serves to document the paths to the materials and
% output of the use of the Matlab camera lens calibration app to
% characterize the intrinsics of the IMX219 camera chip equipped with the
% M12 wide-angle lens. A README file located in the
% data/fisheyeLensCalibration directory provides details regarding the
% measurement.

% Housekeeping
clear
close all

% Path to the matlab camera calibration app session file
sessionFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'fisheyeLensCalibration',...
    'intrinsics_calibration_session.mat');

% Launch the camera calibration app
cameraCalibrator(sessionFileName)

% The resulting intrinsics file is saved as:
%   /derived/arducamB0392cameraInstrinsics.mat
