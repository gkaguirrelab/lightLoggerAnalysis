function I_radiance = imputeSaturatedPixelValues(I,minispectValue)
% This function imputes absolute values for saturated pixels

persistent deltaSteradians
if isempty(deltaSteradians)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'deltaSteradians.mat');
    load(paramFileName,'deltaSteradians');
end

persistent T
if isempty(T)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(paramFileName,'T');
end
cameraWls = T.wls;

% Derive an estimate of the environmental SPD from the minispect
%% just assuming a uniform distribution for now
msSPD = ones(size(cameraWls));

% Obtain an estimate of the mean scene radiance as observed by the
% minispect
avgSceneRadiance = 4;

% Obtain the predicted camera channel weights given the environmental SPD
channels = {'red','green','blue'};
for cc = 1:3
    sensSPD = T.(channels{cc});
    sensSPD = sensSPD / max(sensSPD);
Wrgb(cc) = sensSPD' * msSPD;
end

% Obtain the visual angle area (in steradians) of the saturated pixels
infIdx = isinf(I);
OmegaSaturated = sum(deltaSteradians(infIdx));

% Obtain the total radiance as seen by the unsaturated pixels, weighted by
% the visual angle area of each of these pixels
unsatPartition = sum(I(~infIdx) .* deltaSteradians(~infIdx));

% Calculate the mean radiance that should be present in the saturated pixels


end
