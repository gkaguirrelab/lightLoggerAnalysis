% Run all of the world camera define functions to create the set of derived
% camera parameters. These steps must be run in the specified order.

% Housekeeping
clear all
close all

defineDarkSignal
defineFullWellCapacityEffect
defineFlatFieldingFunction
defineRadiometricWeights
% defineFisheyeCameraIntrinsics -- This is run in an interactive GUI
defineDeltaSteradians
defineAGCToIntegratedRadiance
defineMinispectRadianceWeights
