% The purpose of this script is to define the relationship between the
% custom AGC settings we use to control the sensitivity of the IMX219
% camera and the true effective integrated radiance of the environment
% as seen by the R, G, and B channels independently.

% Housekeeping
clear

% Define the common wavelength domain (380 to 730 nm with 1 nm spacing)
commonS = [380, 1, 352];
commonWls = SToWls(commonS);

% Load the AGC settings for each ND level
agcData.ndf = 0:3;
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

% Spline the sensor sensitivities to the common 1 nm wavelength domain
% and scale so maximum value is unity
sensorSensitivities = zeros(length(commonWls), 3);
for cc = 1:length(channelNames)
    sens = SplineRaw(wlsSensor, T.(channelNames{cc}), commonWls);
    sensorSensitivities(:, cc) = sens ./ max(sens);
end

% Preallocate an Nx3 array for the sensor-weighted integrated radiance
integratedRadiance = zeros(length(agcData.ndf), 3);

% Prepare to plot the chromaticity diagrams for each ND level spectrum
figure
tiledlayout

% Next, load radiance spectrum associated with each ND level
for ii = 1:length(agcData.ndf)

    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'AGCSettingsByNDF',...
        'PR670',...
        sprintf('AGCSettingsMeasure%dNDF.mat',agcData.ndf(ii)));
    load(dataFileName,'measurement','S');
    
    % This is the average radiance in units of Watts/m2/sr/[S(2)*nm]
    spdSource_raw(ii,:) = mean(measurement,1);
    
    % Resample the source SPD to 1 nm spacing. SplineSpd handles the
    % power adjustment from 2 nm bands to 1 nm bands automatically.
    spdSource = SplineSpd(SToWls(S), spdSource_raw(ii,:)', commonWls)';

    % Report the spd chromaticity and luminance in a figure
    nexttile
    luminance(ii) = plotChromLum(spdSource',WlsToS(commonWls));

    % Loop over the channels to calculate the sensor-weighted effective radiance.
    for cc = 1:length(channelNames)
        % Calculate absolute integrated radiance for this channel using the 1 nm dot product
        integratedRadiance(ii,cc) = spdSource * sensorSensitivities(:, cc);
    end
end

% Plot the measurements
figure;
yyaxis left
for cc = 1:3
    loglog(cameraScore, integratedRadiance(:,cc),['o-' channelCodes{cc}],'LineWidth',2,'MarkerSize',10); 
    hold on
end
ylabel('Log integrated radiance (W/m^2/sr)');

% Show a couple validation measurements
loglog([1.6808e+05,3.0757e+04],[0.0576,0.4555],'*k');

% Add the luminance values to the right y-axis
yyaxis right
loglog(cameraScore, luminance,'.','MarkerSize',1); 
ylabel('Log luminance (cd/m^2)');

% General plot properties
a = gca();
a.XScale = 'log';
a.YScale = 'log';
a.TickDir = 'out';
hold on; grid off; box off;

% Clean up, label, legend
xlabel('Log camera sensitivity score');
title('Integrated Radiance vs. Camera AGC Sensitivity');
legend('Red Channel', 'Green Channel', 'Blue Channel', 'Location', 'northwest');


% Save the values that relate camera score to channel-specific effective radiance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'cameraScoreToIntegratedRadiance.mat');
readme = ['Created by defineAGCToIntegratedRadiance.\n'...
    'A linear interpolation between these values (in log10 space) maps AGC values to integrated radiance.\n',...
    'cameraScore -- the product of the AGC settings (analog gain, digital gain, exposure).\n',...
    'effectiveRadiance -- the integrated radiance (W/m2/sr) seen by the R, G, and B channels (Nx3 matrix).\n'];
save(saveFileName,'readme','cameraScore','integratedRadiance');