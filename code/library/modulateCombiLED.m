function modulateCombiLED(frequency, cal_path, color_profile, contrast)
% Modulate the combiLED indefinitely at a desired frequency (hz)
%
% Syntax:
%    modulateCombiLED(frequency, cal_path)
%
% Description:
%   Set the CombiLED to modulate indefinitely at a given frequency.
%   The CombiLED is initialized with the cal_path calibration file. 
%
% Inputs:
%   frequency             - Float. Represents the frequency at which to modulate
%
%   cal_path              - String. Represents the path to the light source
%                           calibration file.  
%   color_profile         - Vector. Represents the color of the light used in
%                           the modulation. Default is pure white. 
%   contrast              - Float. Represents the contrast of the modulation.     
%
% Outputs:
%   None
%
% Examples:
%{
    [~, calFileName, calDir] = selectCal();
    cal_path = fullfile(calDir,calFileName);
    frequency = 5; 
    color_profile = [1,1,1,1,1,1,1,1];
    modulateCombiLED(frequency, cal_path, color_profile);
%}

    % Validate arguments
    arguments
        frequency (1,1) {mustBeNumeric} % Frequency of the modulation 
        cal_path (1,:) {mustBeText} % Path to calibration file
        color_profile (8,1) {mustBeVector} = [1,1,1,1,1,1,1,1]; % Default color profile is pure white
        contrast (1,1) {mustBeNumeric} = 0.5; % The contrast of the modulation
    end
    
    % First, get the path to where the calibration files are
    calibration_filepath = fullfile(getpref("lightLogger", "light_logger_analysis_path"), "cal", "Sphere", "CombiLED-A_cassette-ND0_sphere.mat");

    % Load the required stuff for the cal file
    tbUseProject('combiExperiments');
    
    % Load in the calibration file for the CombiLED
    load(calibration_filepath,'cals');
    cal = cals{end};

    % Initialize the combiLED
    disp('Opening connection to CombiLED...')
    CL = CombiLEDcontrol(); % Initialize CombiLED Object
    CL.setGamma(cal.processedData.gammaTable);  % Update the combiLED's gamma table

    % Collect information to compose flicker profile
    observerAgeInYears = 30;
    pupilDiameterMm = 3;
    photoreceptors = photoreceptorDictionaryHuman('observerAgeInYears',observerAgeInYears,'pupilDiameterMm',pupilDiameterMm);

    % Compose flicker profile
    modResult = designModulation('LightFlux',photoreceptors,cal);

    % Set the color profile of the modResult
    modResult.settingsHigh = color_profile; 
    modResult.settingsBackground = color_profile./2; 
    modResult.settingsLow = color_profile.*0; 

    % Initialize other settings of the wave
    CL.setSettings(modResult);
    CL.setWaveformIndex(1);
    CL.setContrast(contrast);

    % Set the CL flicker to desired frequency
    CL.setFrequency(frequency);

    % Start flickering 
    disp('Beginning modulation...')
    CL.startModulation();

    % Pause while desired flicker
    disp('Press any key to stop modulation.')
    pause()
    
    % Close connection to the combiLED
    CL.serialClose();


end