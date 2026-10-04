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
defineCameraToVisualAngles
defineAGCToIntegratedRadianceViaMacbeth
defineSensorToConeMapping

% Not strictly part of the IMX219 camera calibration pathway, but
% nonetheless lives here.
defineMinispectRadianceWeights

% We can define the AGC -> integrated radiance conversion using either:
%{
    defineAGCToIntegratedRadianceViaSphere
    defineAGCToIntegratedRadianceViaMacbeth
%}
% The first of these has a stronger theoretical motivation. The definition
% via sphere assumes a fixed -1 (log) slope between integrated radiance and
% camera sensitivity, and then derives the intercept from an analysis of
% calibration measurements made in the the light integrating sphere.
%
% The second routine conducts a search across possible function parameters
% to find the conversion that best aligns the set of measured radiance
% values from the macbeth color checker in different environments to the
% predicted radiance as measured with the PR670.