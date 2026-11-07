% defineSensorToConeMatrix.m
%
% This script derives the 3x3 transformation matrix to convert IMX219
% camera RGB sensor values into LMS cone isomerization rates. It bases the
% cone sensitivities on 2-degree foveal fundamentals for a standard
% observer. It also derives a 1x3 vector to convert RGB to photopic
% luminance (V_lambda). It saves the resulting matrices to the "derived"
% directory.

clear all
close all

% Load the world camera channel spectral sensitivity functions
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName, 'T');

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing)
commonS = [380, 1, 352];

% Extract wavelength support and spline sensitivities to the common domain
sensorWls = T{:, 1};
T_common = SplineRaw(sensorWls, [T.red, T.green, T.blue], SToWls(commonS));

% Scale the camera spectral sensitivity functions so their maximum value is unity
sensorSensitivities = T_common ./ max(T_common, [], 1);

% Get the sensitivities for the foveal cone classes using the common domain
fieldSizeDegrees = 2;
observerAgeInYears = 30;
pupilDiameterMm = 2;
T_receptors = GetHumanPhotoreceptorSS(commonS, ...
    {'LConeTabulatedAbsorbance2Deg', ...
    'MConeTabulatedAbsorbance2Deg', ...
    'SConeTabulatedAbsorbance2Deg'}, ...
    fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], [], [], []);

% Calculate the 3x3 transformation matrix using the pseudoinverse
% T_receptors is 3xN, sensorSensitivities is Nx3 (so its transpose is 3xN)
rgb2lmsMatrix = T_receptors * pinv(sensorSensitivities');

% Load the CIE 1931 XYZ color matching functions to get the luminance
% function (V_lambda)
load T_xyz1931

% T_xyz1931(2,:) corresponds to the y-bar (V_lambda) function
T_vLambda_common = SplineCmf(S_xyz1931, T_xyz1931(2,:), commonS);

% Define peak luminous efficacy (Km) in lumens per watt
Km = 683;

% Calculate the 1x3 transformation vector to convert RGB to absolute
% Luminance (Scaled by Km to yield photometric units rather than a
% normalized response)
rgb2lumVector = Km * (T_vLambda_common * pinv(sensorSensitivities'));

% Save the resulting matrix to the derived directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'sensorToConeMatrix.mat');

readme = ['Created by defineSensorToConeMatrix.\n'...
    'rgb2lmsMatrix -- 3x3 matrix to project IMX219 RGB radiance to LMS cone fundamentals.\n' ...
    'rgb2lumVector -- 1x3 vector to project IMX219 RGB radiance to photopic luminance (V_lambda).\n'];

save(saveFileName, 'readme', 'rgb2lmsMatrix', 'rgb2lumVector');

fprintf('\nTransformation matrices derived and saved to:\n%s\n\n', saveFileName);