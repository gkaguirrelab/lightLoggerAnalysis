%% measureCombiLEDEffectiveRefreshRate
%
% Add some comments here



fullPath = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_admin/combiLED_flicker-BT_measures/CombiA_FlickerCalibrate_20250623.csv';

%[fileName, fileDir] = uigetfile('.csv');
%fullPath = fullfile(fileDir,fileName);

% Load the frequency measurement data set
for ii = 1:13
    [meta, signal(:,ii)] = extractLabSphereData(fullPath,ii);
    fHz(ii) = round(meta.FundamentalFrequencyHz);
    dt = 1/meta.SampleRateHz;
    tSecs = 0:dt:(length(signal(:,ii))-1)*dt;
    [~,amplitudes(ii)] = fourierRegression( signal(:,ii), tSecs, fHz(ii) );
end    

% Convert to relative amplitude
amplitudes = amplitudes / amplitudes(1);

% Fit the amplitude data with a model of an ideal discrete device producing
% sinusoidal modulations. The filterAmp value should be close to unity.
[sourceRefreshRate,filterAmp] = approxFreqFilter(fHz,amplitudes);

% Plot the data and the ideal device fit
figure
semilogx(fHz,amplitude,'*k');
hold on
fitFreqs = logspace(log10(1),log10(150));
filterProfile = idealDiscreteSampleFilter(fitFreqs,1/sourceRefreshRate);
semilogx(fitFreqs,filterProfile,'-r');

% Report the refresh rate of the source
fprintf('The effective refresh rate of the source is %2.2f Hz\n',sourceRefreshRate);



%% LOCAL FUNCTIONS

function [filterFreqHz,filterAmp] = approxFreqFilter(frequencies,amplitudes)

myObj = @(p) norm(amplitudes - p(1)*idealDiscreteSampleFilter(frequencies,1/p(2)));
x0 = [1,205];
lb = [0, 75];
ub = [1, 300];
f = fmincon(myObj,x0,[],[],[],[],lb,ub);
filterFreqHz = f(2);
filterAmp = f(1);


end