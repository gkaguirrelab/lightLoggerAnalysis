function [minispectValues, meanSpectralRadiance] = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS)
% ESTIMATEMINISPECTVALUESFROMRADIANCEMODEL Computes expected ASM7341 sensor counts
% from a hemispheric radiance model and its wavelength sampling S.
%
%   miniSpectValues = estimateMinispectValuesFromRadianceModel(radianceModel, radianceModelS)

persistent miniSpectT miniSpectWls miniSpectKVals

if isempty(miniSpectT)
    % Load spectral sensitivity
    sensFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'ASM7341_spectralSensitivity.mat');
    load(sensFileName, 'T');
    miniSpectWls = T.wl;
    rawT = table2array(T(:, ["F1" "F2" "F3", "F4", "F5", "F6", "F7", "F8", "Clear"]))';
    % Max-normalize sensitivity matrix matching calibration dot product
    miniSpectT = rawT ./ max(rawT, [], 2);
end

if isempty(miniSpectKVals)
    % Load calibration weights / coefficients
    weightFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'minispectRadianceWeights.mat');
    load(weightFileName, 'fitObj');
    nChannels = size(miniSpectT, 1);
    miniSpectKVals = fitObj.coeff(1:nChannels, 2);
end

% Get wavelengths from radianceModelS
modelWls = SToWls(radianceModelS);

% Interpolate sensitivity matrix A if model wavelengths differ from minispect wavelengths
if length(modelWls) ~= length(miniSpectWls) || any(modelWls ~= miniSpectWls)
    nChannels = size(miniSpectT, 1);
    A = nan(nChannels, length(modelWls));
    for i = 1:nChannels
        A(i, :) = interp1(miniSpectWls, miniSpectT(i, :), modelWls, 'linear', 0);
    end
else
    A = miniSpectT;
end

% 1. Evaluate hemispheric mean spectral radiance from radianceModel
gridRes = 200;
azGrid = linspace(-pi, pi, gridRes * 2);
elGrid = linspace(-pi / 2, pi / 2, gridRes);
[AZ, EL] = meshgrid(azGrid, elGrid);

% Evaluate radianceModel on the spatial grid ([numWls, size(AZ)])
L_grid = radianceModel(AZ, EL);

% 2. Calculate cosine-weighted mean spectral radiance over the hemisphere
weights = cos(EL);
sumWeights = sum(weights(:));

numWls = length(modelWls);
L_grid_reshaped = reshape(L_grid, [numWls, numel(weights)]);
weights_row = reshape(weights, [1, numel(weights)]);

meanSpectralRadiance = (L_grid_reshaped * weights_row') / sumWeights; % [numWls, 1]

% 3. Calculate integrated radiance (dot product with max-normalized sensitivity A)
y = A * meanSpectralRadiance; % [nChannels, 1]

% 4. Convert integrated radiance to expected sensor counts using calibration weights k
k = miniSpectKVals;
minispectValues = y .* (10 .^ k);

% Return as an integer row vector matching standard sensor reading formats
minispectValues = round(minispectValues(:)');
end