function MS2illum_lux = msCounts2Illuminance(ms_counts)
% msCountsToIlluminance - Converts MiniSpectrometer readings to illuminance
% (lux) using fits saved from linear fit of illuminance to MS readings
% during calibration. Averages across all channels.
%
% Syntax:
%   MS2illum_lux = msCounts2Illuminance(ms_counts, illum_to_MS)
%
% Inputs:
%   ms_counts    - [nSamples x nChannels] matrix of MS counts (can be linear or log10)
%
% Output:
%   illum_lux    - [nSamples x nChannels] estimated illuminance in lux
%

% load coefficients for converting MS chip values to illuminance
combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');
MS_illumfile = [combiExperiments_path, '/data/PR670_illum_to_MS_fits'];
load(MS_illumfile, 'illum_to_MS');

% Ensure data is in log10 space
if any(ms_counts(:) <= 0)
    error('ms_counts must be strictly positive to apply log10.');
end
ms_log = log10(ms_counts);

% Prepare output
[nSamples, nChannels] = size(ms_log);

if size(illum_to_MS, 1) >= nChannels
    minChannels = min(size(illum_to_MS,1), nChannels);
    warning('Size mismatch: ms_counts has %d channels, but illum_to_MS has %d. Using channels 1- %d.',...
        nChannels, size(illum_to_MS,1), minChannels);
end

MS2illum_lux_by_channel = nan(nSamples, nChannels);

% Invert the regression per channel
for ch = 1:minChannels
    m = illum_to_MS(ch, 1);
    b = illum_to_MS(ch, 2);

    % Invert the linear relationship
    log10_illum = (ms_log(:, ch) - b) ./ m;

    % Convert back to lux
    MS2illum_lux_by_channel(:, ch) = 10.^log10_illum;
end
% average across all the channels
MS2illum_lux = mean(MS2illum_lux_by_channel,2);
end