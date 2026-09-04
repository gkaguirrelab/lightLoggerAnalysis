function [spectralRadiance,S] = estimateRadianceSpectrumFromMiniSpect(minispectValue)
% This function estimates the mean environmental radiance from the
% minispect values

miniSpectWeights = [1 2 1 1 2 2 2 2];

persistent T
if isempty(T)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'ASM7341_spectralSensitivity.mat');
    load(paramFileName,'T');
end
miniSpectWls = T.wl;

% Loop through the channels, find the lambda max and the max sensitivity
% value
for ii = 1:8
    vec(ii,:) = T.(sprintf('F%d',ii));
    [peakSensitivityValue(ii),idxVals(ii)] = max(vec(ii,:));
    lambdaMax(ii) = miniSpectWls(idxVals(ii));
end

% Assemble an estimated radiance spectrum. The estimate starts with the
% weighted sensitivity spectra of the first and last channel
spectralRadiance = miniSpectWeights(1) * vec(1,:) / peakSensitivityValue(1);
spectralRadiance = spectralRadiance + miniSpectWeights(8) * vec(8,:) / peakSensitivityValue(8);

% Create a spline between weights at the lambda max values
x = lambdaMax;
y = miniSpectWeights;
values = csaps(x,y,1,lambdaMax(1):lambdaMax(8));
spectralRadiance(idxVals(1):idxVals(8)) = values;

end