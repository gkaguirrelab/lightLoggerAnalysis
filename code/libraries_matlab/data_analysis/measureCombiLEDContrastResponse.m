%% measureCombiLEDContrastResponse
% 
% Analyzes a Labsphere contrast response dataset from a combiLED calibration session. 
% Extracts the fundamental response amplitude for modulations at different temporal frequencies
% using Fourier regression, and then fits a linear model to the amplitude-contrast relationship.


% fullPath = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_admin/combiLED_flicker-BT_measures/CombiA_FlickerCalibrate_20250623.csv';

[fileName, fileDir] = uigetfile('.csv');
fullPath = fullfile(fileDir,fileName);

% Make sure to clear and preallocate arrays
contrasts = zeros(1,23);
amplitudes = zeros(14,23);
fHz = zeros(14,23);

% Load the contrast response data set (scans 14-23)
for ii = 14:23
    [meta, signal(:,ii)] = extractLabSphereData(fullPath,ii);
    fHz(ii) = round(meta.FundamentalFrequencyHz);

    dt = 1/meta.SampleRateHz;
    tSecs = 0:dt:(length(signal(:,ii))-1)*dt;
    
    [~,amplitudes(ii)] = fourierRegression( signal(:,ii), tSecs, fHz(ii) );

    % Obtain contrasts from comment
    comment = meta.Comment;
    contrast = regexpi(comment, '(0\.\d+) contrast', 'tokens');
    if ~isempty(contrast)
        contrasts(ii) = str2double(contrast{1}{1});
    else
        warning('No contrast value found in comment for scan %d.', ii);
    end

end    

% Find linear model
p = polyfit(contrasts(14:23), amplitudes(14:23), 1);
slope = p(1);
intercept = p(2);

%% Make a plot of contrast vs. amplitude, and fit with linear model
plot(contrasts(14:23), amplitudes(14:23), '*b')
hold on
plot(contrasts(14:23), polyval(p, contrasts(14:23)))

xlabel('Contrast')
ylabel('Amplitude')
title('CombiLED Contrast Response Function')
legend('Data', 'Linear Fit')

% Report the params of the linear fit
fprintf('Slope: %.4f\n', slope);
fprintf('Intercept: %.4f\n', intercept);
