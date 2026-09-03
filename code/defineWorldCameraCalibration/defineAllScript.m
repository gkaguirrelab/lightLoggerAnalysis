% Run all of the world camera define functions to create the set of derived
% camera parameters. These steps must be run in the specified order.

% Housekeeping
clear
close all
addpath(fullfile(fileparts(mfilename('fullpath')),'utilities'))

defineDarkSignal
defineFullWellCapacityEffect
defineFlatFieldingFunction
defineRadiometricWeights
% defineFisheyeCameraIntrinsics -- This is run in an interactive GUI
defineDeltaSteradians
defineAGCToMeanRadiance

validateDerivedFilesHaveREADME(fullfile( ...
    tbLocateProjectSilent('lightLoggerAnalysis'), 'derived'));
