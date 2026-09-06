% The purpose of this script is to define the relationship between the
% custom AGC settings we use to control the sensitivity of the IMX219
% camera, the irradiance to which the camera is exposed, and the true 
% broad-band radiance of the environment.
%
% We observe a slight difference between the color channels in their
% effective radiometric sensitivity. These differences can arise from 
% several sources (interaction of the source spectrum with the particular 
% sensitivity functions of the channels; differences in radiometric 
% sensitivity between the channels). However, because we must map the AGC
% settings to the true energy-conserved radiance of the scene, we calculate
% the true integrated broad-band radiance directly from the source spectrum.

% Housekeeping
clear

% Use the full-precision fixed AGC settings for the 0.25 contrast
% calibration. These values should remain:
%   NDF       = [0,          1,          2,          3,          4,          5]
%   AGain     = [1,          1.80281687, 10.666,     10.666,     10.666,     10.666]
%   DGain     = [1,          1,          1.58156674, 5.96331605, 7.97561560, 10]
%   Exposure  = [1466,       8333,       8333,       8333,       8333,       8333]
%   AGC score = [1466,       15022.87298, 140569.30074, 530018.20667, 708870.94394, 888797.78]
agcData.ndf = 0:5;
agcData.AGain = [1, 1.80281687, 10.666, 10.666, 10.666, 10.666];
agcData.DGain = [1, 1, 1.58156674, 5.96331605, 7.97561560, 10];
agcData.Exposure = [1466, 8333, 8333, 8333, 8333, 8333];

% Derive a "camera score" by obtaining the product of the AGC settings
cameraScore = agcData.DGain .* agcData.AGain .* agcData.Exposure;

% Load the IMX219 sensitivity functions.
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');
wlsSensor = T.wls;
channelNames = {'red','green','blue'};
channelCodes = {'r','g','b'};

% Preallocate array for the true broad-band radiance
avgSceneRadiance = zeros(length(agcData.ndf), 1);

% Next, load the "maxSpectrum" calibration files and extract for each the
% radiance of the sphere interior. Account for the settings level of the
% combiLED during the measurement.
settings = 0.25;
for ii = 1:length(agcData.ndf)
    thisCalFile = sprintf('CombiLED-A_cassette-ND%d_sphere_maxSpectrum.mat',agcData.ndf(ii));
    load(thisCalFile,'cals');
    cal = cals{end};
    S = cal.rawData.S;
    wlsSource = SToWls(S);
    % This is the average radiance in units of Watts/m2/sr/[S(2)*nm],
    % adjusted for the settings level of the CombiLED.
    spdSource = cal.rawData.gammaCurveMeanMeasurements * settings;
    
    % Calculate the true integrated broad-band radiance by 
    % summing the spectral power distribution across all wavelength bands.
    avgSceneRadiance(ii) = sum(spdSource);
    
    % Loop over the channels to calculate the sensor-weighted effective 
    % radiance (retained exclusively for plotting comparisons).
    for cc = 1:length(channelNames)
        % Spline the sensor sensitivity to match the source SPD
        sensitivitySensor = SplineRaw(wlsSensor,T.(channelNames{cc}),wlsSource);
        sensitivitySensor = sensitivitySensor./max(sensitivitySensor);
        effectiveRadiance(ii,cc) = spdSource' * sensitivitySensor;
    end
end

% Plot the measurements
figure;
for cc = 1:3
    loglog(cameraScore, effectiveRadiance(:,cc),['-.' channelCodes{cc}],'LineWidth',1,'MarkerSize',10); 
    hold on
end
loglog(cameraScore, avgSceneRadiance,'-*k','LineWidth',2,'MarkerSize',10); 
a = gca();
a.XScale = 'log';
a.YScale = 'log';
a.TickDir = 'out';
hold on; grid off; box off;

% Clean up, label, legend
xlabel('Log camera sensitivity score');
ylabel('Log mean radiance (W/m2/sr)');
title('Average radiance vs. camera AGC sensitivity');

% Save the values that relate camera score to true broad-band radiance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'cameraScoreToAverageRadiance.mat');
readme = ['Created by defineAGCToMeanRadiance.\n'...
    'A linear interpolation between these values (in log10 space) maps AGC values to radiance.\n',...
    'cameraScore -- the product of the AGC settings (analog gain, digital gain, exposure).\n',...
    'avgSceneRadiance -- the true broad-band radiance (W/m2/sr) of the scene viewed by the camera.\n'];
save(saveFileName,'readme','cameraScore','avgSceneRadiance');