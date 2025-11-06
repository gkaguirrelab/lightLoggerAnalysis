
% Local function that will handle initalizing the CombiLED 
function [CombiLED, cal] = initialize_combiLED(combiExperiments_path, NDF) 
    arguments 
        combiExperiments_path {mustBeText}; 
        NDF {mustBeNumeric}; 
    end 

    % First, get the preference for where the light logger analysis dir is. 
    % We will use this to find the cal files 
    light_logger_analysis_dir = getpref("lightLogger", "light_logger_analysis_path"); 

    % Ensure we are using the correct library 
    tbUse("combiLEDToolbox");

    % Connect to the CombiLED and initialize
    % the initial background state
    CombiLED = CombiLEDcontrol();
    CombiLED.verbose = 0; 

    % Adjust the calibration file for this NDF level
    baseCal_path = fullfile(light_logger_analysis_dir, "cal", "Sphere", "CombiLED-A_cassette-ND0_sphere.mat");
    baseCals = load(baseCal_path, "cals").cals;  % First, load in the base cal
    baseCal = baseCals{end}; 

    baseCalsMaxSpectrum_path = fullfile(light_logger_analysis_dir, "cal", "Sphere", "CombiLED-A_cassette-ND0_sphere" + "_maxSpectrum" + ".mat");
    baseCalsMaxSpectrum = load(baseCalsMaxSpectrum_path, "cals").cals; % Load in the corresponding max spectrum 
    baseCalMax = baseCalsMaxSpectrum{end}; 

    targetCalsMaxSpectrum_path = fullfile(light_logger_analysis_dir, "cal", "Sphere", "CombiLED-A_cassette-ND" + string(NDF) + "_sphere" + "_maxSpectrum" + ".mat");
    targetCalsMaxSpectrum = load(targetCalsMaxSpectrum_path, "cals").cals; % Load in the maxs spectrum corresponding to the target ND level
    targetCalMax = targetCalsMaxSpectrum{end}; 

    % Obtain the transmittance for this ND filter setting
    transmittance = targetCalMax.rawData.gammaCurveMeanMeasurements ./ baseCalMax.rawData.gammaCurveMeanMeasurements;
    
    % basecal is my original one in cal_path
    % Modify this cal file to account for ND transmittance
    cal = baseCal;
    for ii = 1:size(cal.processedData.P_device, 2)
        cal.processedData.P_device(:,ii) = ...
            cal.processedData.P_device(:,ii) .* transmittance;
    end
    cal.processedData.P_ambient = cal.processedData.P_ambient .* transmittance;

    CombiLED.setGamma(cal.processedData.gammaTable);

    % By default, gamma correction is not performed in directMode. Here we turn
    % gamma correction on.
    CombiLED.setDirectModeGamma(true);

    return ; 
end