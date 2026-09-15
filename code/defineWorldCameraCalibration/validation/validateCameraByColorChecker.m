% Validation

% Housekeeping
clear

% Get the list of spectral radiometric measurements of checks
dirName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'macbethColorCheck',...
    'PR670',...
    '*.mat');
fileList = dir(dirName);

% Load the spectral measurements
for ii = 1:length(fileList)
    fileName = fullfile(fileList(ii).folder,fileList(ii).name);
    load(fileName,'measurement','S')
    myIndex = int32(sscanf(fileList(ii).name, 'Index-%d'));
    [col, row] = ind2sub([6 4], myIndex);
    spectralRadiance{col, row} = mean(measurement,1);
end
spectralRadianceS = S;

% Get the table of reflectance spectra of the macbeth color checker
[spectralReflectance,spectralReflectanceS] = loadMacbethReflectance();

% Define the common wavelength domain (380 to 730 nm with 2 nm spacing).
% Number of samples: (730 - 380)/2 + 1 = 176 samples.
commonS = [380, 2, 176]; 

% Initialize an array to hold the estimated illuminant spectrum from each measured patch
illuminantEstimates = [];

% 1. Estimate the shared illuminant
for c = 1:6
    for r = 1:4
        if ~isempty(spectralRadiance{c, r})
            % Convert to W/m2/sr/nm and ensure column vector format
            rad_nm = spectralRadiance{c, r}(:) / 2;
            ref_patch = spectralReflectance{c, r}(:);

            % Resample both radiance and reflectance to the common wavelength domain
            rad_common = SplineSpd(SToWls(spectralRadianceS), rad_nm, SToWls(commonS));
            ref_common = SplineSpd(SToWls(spectralReflectanceS), ref_patch, SToWls(commonS));

            % Estimate effective illuminant
            illuminantEstimates(:, end+1) = rad_common ./ ref_common;
        end
    end
end

% Average across the 5 measured patches
sharedIlluminant = mean(illuminantEstimates, 2);

% 2. Predict spectral radiance for unmeasured patches
for c = 1:6
    for r = 1:4
        if isempty(spectralRadiance{c, r})
            ref_patch = spectralReflectance{c, r}(:);
            ref_common = SplineSpd(SToWls(spectralReflectanceS), ref_patch, SToWls(commonS));

            % Predict radiance in W/m2/sr/nm, convert back to measurement units (* 2)
            spectralRadiance{c, r} = (sharedIlluminant .* ref_common) * 2;
        else
            rad_nm = spectralRadiance{c, r}(:) / 2;
            rad_common = SplineSpd(SToWls(spectralRadianceS), rad_nm, SToWls(commonS));
            spectralRadiance{c, r} = rad_common * 2;
        end
    end
end

% Update the radiance S vector to reflect the new common wavelength domain
spectralRadianceS = commonS;



% UNUSED
% 
% % Load the data for the "close" camera acquisition of the color checker 
% dataFileName = fullfile(...
%     tbLocateProjectSilent('lightLoggerAnalysis'),...
%     'data',...
%     'macbethColorCheck',...
%     'lightLogger',...
%     'close_AGCandMS_01.mat');
% load(dataFileName,'worldFrame','AGCSettings');
%
% % Load the world camera channel spectral sensitivity functions. This is a
% % table with the first column providing the wavelength support.
% dataFileName = fullfile(...
%     tbLocateProjectSilent('lightLoggerAnalysis'),...
%     'data',...
%     'IMX219_spectralSensitivity.mat');
% load(dataFileName,'T');
