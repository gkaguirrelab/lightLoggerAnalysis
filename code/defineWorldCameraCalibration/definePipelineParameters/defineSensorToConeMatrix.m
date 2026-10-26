% defineSensorToConeMatrix.m
%
% This script derives the 3x3 transformation matrix to convert IMX219 
% camera RGB sensor values into LMS cone isomerization rates. It bases the 
% cone sensitivities on 2-degree foveal fundamentals for a standard observer 
% and saves the resulting matrix to the "derived" directory.

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

% Save the resulting matrix to the derived directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'sensorToConeMatrix.mat');

readme = ['Created by defineSensorToConeMatrix.\n'...
    'rgb2lmsMatrix -- 3x3 matrix to project IMX219 RGB radiance to LMS cone fundamentals.\n'];

save(saveFileName, 'readme', 'rgb2lmsMatrix');

fprintf('\nTransformation matrix derived and saved to:\n%s\n\n', saveFileName);