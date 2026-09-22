% The purpose of this script is to define the relationship between the
% custom AGC settings we use to control the sensitivity of the IMX219
% camera and the true effective integrated radiance of the environment
% as seen by the R, G, and B channels independently.

% Housekeeping
clear

% Load the AGC settings for each ND level
agcData.ndf = 0:4;
for ii = 1:length(agcData.ndf)
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'worldCamera',...
        sprintf('NDF%d',agcData.ndf(ii)),...
        sprintf('NDF%d_AGCandMS_01.mat',agcData.ndf(ii)));
    load(dataFileName,'AGCSettings');
    agcData.AGain(ii) = AGCSettings.Again;
    agcData.DGain(ii) = AGCSettings.Dgain;
    agcData.Exposure(ii) = AGCSettings.exposure;
end

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

% Preallocate an Nx3 array for the sensor-weighted integrated radiance
integratedRadiance = zeros(length(agcData.ndf), 3);

% Next, load radiance spectrum associated with each ND level
for ii = 1:length(agcData.ndf)

    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'PR670',...
        sprintf('AGCSettingsMeasure%dNDF.mat',agcData.ndf(ii)));
    load(dataFileName,'measurement','S');
    wlsSource = SToWls(S);

    % This is the average radiance in units of Watts/m2/sr/[S(2)*nm],
    spdSource = mean(measurement,1);

    % Loop over the channels to calculate the sensor-weighted effective radiance.
    for cc = 1:length(channelNames)
        % Spline the sensor sensitivity to match the source SPD
        sensitivitySensor = SplineRaw(wlsSensor,T.(channelNames{cc}),wlsSource);

        % Scale sensor sensitivity so maximum value is unity
        sensitivitySensor = sensitivitySensor ./ max(sensitivitySensor);

        % Calculate absolute integrated radiance for this channel
        integratedRadiance(ii,cc) = spdSource * sensitivitySensor;
    end
end

% Plot the measurements
figure;
for cc = 1:3
    loglog(cameraScore, integratedRadiance(:,cc),['-' channelCodes{cc}],'LineWidth',2,'MarkerSize',10); 
    hold on
end
a = gca();
a.XScale = 'log';
a.YScale = 'log';
a.TickDir = 'out';
hold on; grid off; box off;

% Clean up, label, legend
xlabel('Log camera sensitivity score');
ylabel('Log integrated radiance (W/m^2/sr)');
title('Integrated Radiance vs. Camera AGC Sensitivity');
legend('Red Channel', 'Green Channel', 'Blue Channel', 'Location', 'northwest');

% Save the values that relate camera score to channel-specific effective radiance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'cameraScoreToIntegratedRadiance.mat');
readme = ['Created by defineAGCToMeanRadiance.\n'...
    'A linear interpolation between these values (in log10 space) maps AGC values to integrated radiance.\n',...
    'cameraScore -- the product of the AGC settings (analog gain, digital gain, exposure).\n',...
    'effectiveRadiance -- the integrated radiance (W/m2/sr) seen by the R, G, and B channels (Nx3 matrix).\n'];
save(saveFileName,'readme','cameraScore','integratedRadiance');