function coneVec = cameraToCones(rgbVec, camera)
% Convert an n x 3 matrix of camera rgb measurements to an n x 3 matrix of
% relative lms cone weights. This is all quite crude at the moment.

% Table that contains the average spectral sensitivity functions of the
% camera sensors in commercial mobile phones, reported in:
%   Tominaga S, Nishi S, Ohtera R. Measurement and estimation of spectral
%   sensitivity functions for mobile phone cameras. Sensors. 2021 Jul
%   22;21(15):4985.

% Default
if nargin < 2
    camera = 'standard';
end

% Load camera spectral sensitivity functions
switch lower(camera)
    case 'standard'
        filename = fullfile(fileparts(mfilename('fullpath')), 'averageCameraSPDs.xlsx');
        Table_camera = readtable(filename);
        wls = Table_camera.wavelength;
        T_camera = table2array(Table_camera(:,2:4))';

    case 'imx219'
        filename = fullfile(fileparts(mfilename('fullpath')), 'IMX219_spectralSensitivity.mat');
        data = load(filename);
        Table_camera = data.T;
        wls = Table_camera.wls;
        T_camera = table2array(Table_camera(:,{'red', 'green', 'blue'}))';

    otherwise
        error('Unknown camera type. Use ''standard'' or ''imx219''.');
end

% Convert wavelengths to sampling format
S = WlsToS(wls);

% Get the sensitivities for the foveal cone classes
fieldSizeDegrees = 30;
observerAgeInYears = 30;
pupilDiameterMm = 2;
T_receptors = GetHumanPhotoreceptorSS(S, ...
    {'LConeTabulatedAbsorbance2Deg', ...
    'MConeTabulatedAbsorbance2Deg', ...
    'SConeTabulatedAbsorbance2Deg'}, ...
    fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], [], [], []);

% Create the spectrum implied by the rgb camera weights, and then project
% that on the receptors
coneVec = (T_receptors*(rgbVec*T_camera)')';

end