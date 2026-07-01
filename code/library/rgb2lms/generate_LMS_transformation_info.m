function [T_receptors, T_camera] = generate_LMS_transformation_info(camera)
% Load camera and cone spectral sensitivities for RGB-to-LMS conversion.
%
% Syntax:
%   [T_receptors, T_camera] = generate_LMS_transformation_info(camera)
%
% Description:
%   This function prepares the two spectral bases needed by `rgb2lms`: one
%   for the camera sensors and one for human L, M, and S cones. Depending
%   on the requested camera type, it reads either tabulated average mobile
%   camera sensitivities or a saved IMX219 sensitivity measurement,
%   converts the wavelength support into the Psychtoolbox `S` format, and
%   then queries the human photoreceptor model on the same sampling grid.
%
% Inputs:
%   camera                   - String. Camera sensitivity model to use:
%                              `"standard"` loads the spreadsheet-based
%                              average phone sensitivities and `"imx219"`
%                              loads the measured IMX219 sensitivities.
%
% Outputs:
%   T_receptors              - 3-by-N matrix of human L-, M-, and S-cone
%                              sensitivities sampled on the shared
%                              wavelength grid.
%   T_camera                 - 3-by-N matrix of red, green, and blue
%                              camera sensor sensitivities sampled on the
%                              same grid.
%
% Examples:
%{
    [T_receptors, T_camera] = generate_LMS_transformation_info("imx219");
%}

%{ 

% Table that contains the average spectral sensitivity functions of the
% camera sensors in commercial mobile phones, reported in:
%   Tominaga S, Nishi S, Ohtera R. Measurement and estimation of spectral
%   sensitivity functions for mobile phone cameras. Sensors. 2021 Jul
%   22;21(15):4985.
% 
%}

    arguments
        camera {mustBeMember(camera, ["imx219", "standard"])}
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


    return ; 
end
