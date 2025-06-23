%% measureCombiLEDEffectiveRefreshRate
%FIx this name
% Add some comments here



fullPath = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_admin/combiLED_flicker-BT_measures/CombiA_FlickerCalibrate_20250623.csv';

%[fileName, fileDir] = uigetfile('.csv');
%fullPath = fullfile(fileDir,fileName);

% Load the contrast response data set
for ii = 14:23
    [meta, signal(:,ii)] = extractLabSphereData(fullPath,ii);
    fHz(ii) = round(meta.FundamentalFrequencyHz);
    dt = 1/meta.SampleRateHz;
    tSecs = 0:dt:(length(signal(:,ii))-1)*dt;
    [~,amplitudes(ii)] = fourierRegression( signal(:,ii), tSecs, fHz(ii) );
end    

%% Make a plot of contrast vs. amplitude,and fit with linear model

% Report the params of the linear fit
fprintf('The effective refresh rate of the source is %2.2f Hz\n',sourceRefreshRate);

