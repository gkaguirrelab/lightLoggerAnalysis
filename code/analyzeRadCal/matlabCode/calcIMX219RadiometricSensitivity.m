% Estimate IMX219 radiometric sensitivity by channel

% We need to obtain the spectrum of a source light and the IMX219 chip
% spectral sensitivity functions.

% Here we load a file that Zach made that contains a calibration file and
% derive the SPD of the light source. This is all a bit idiosyncratic, but
% that's OK for this one-time measurement
load("/Users/aguirre/Desktop/world_ttf_data.mat");

% We will use a particular measurement (NDF1)
measIdx = 2; % This is the 2nd measurement which was NDF=1.
cal = world_ttf_data.metadata.temporal_sensitivity.cal_files{measIdx};

% This cal file had already been adjusted for the presence of the NDF. We
% need to account for the contrast that was used to present the stimulus.
% We assume that the gamma correction was in place, so that the contrast
% may be applied in a linear manner to the spd
contrast = world_ttf_data.metadata.temporal_sensitivity.contrast_levels;
sourceSPD = sum(cal.processedData.P_device,2) * contrast;

% Get the wavelength support for the source
sourceS = cal.rawData.S;

% Load the IMX219 spectral sensitivity functions
filePath = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'),'data','IMX219_spectralSensitivity.mat');
load(filePath,'T');
sensorS = WlsToS(T.wls);

% Interpolate the source SPD to be in the same domain as the sensor
% sensitivity functions. Need to correct for the change from 2 to 1 nm
% sampling
sourceSPDInterp = interp1(SToWls(sourceS),sourceSPD,SToWls(sensorS),"linear","extrap");
sourceSPDInterp = sourceSPDInterp * (sensorS(2) / sourceS(2));

% Obtain the weights for the three channels
channelLabels = {'blue','green','red'};
for cc = 1:length(channelLabels)
    sensorSPD = T.(channelLabels{cc});
    weightBGR(cc) = sourceSPDInterp' * sensorSPD;
end

% This is an example set of observed weights from a video recording of this
% source spectrum
observedBGR = [ 1.14006268, 0.7269785 , 1]; 

% We wish to calculate the factors by which we need to multiple the
% observed BGR weights to result in values that match the predicted BGR
% weights:
correctionBGR = weightBGR ./ observedBGR;

% The absolute value of the correctionBGR vector is not important; instead
% we care about the relative correction weights. Set the mean of the
% correction weights to unity.
correctionBGR = correctionBGR / mean(correctionBGR);

