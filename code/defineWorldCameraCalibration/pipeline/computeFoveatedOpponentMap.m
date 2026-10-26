function [postRecepMap, coneWeights, AzGrid, ElGrid] = computeFoveatedOpponentMap(integratedRadianceMap, options)
% COMPUTEFOVEATEDOPPONENTMAP Converts an integrated RGB radiance map into a 
% foveated, post-receptoral opponent channel map using a uniform RGB-to-LMS 
% conversion based on 2-degree cone fundamentals.
%
% Inputs:
%   integratedRadianceMap - H x W x 3 matrix of RGB radiance
%
% Name-Value Arguments:
%   gazeX, gazeY     - Fixation coordinates (defaults to center of image)
%   azimuthGrid      - Grid spacing for azimuth (default: -90:0.5:90)
%   elevationGrid    - Grid spacing for elevation (default: -90:0.5:90)
%   balanceNeutral   - Dynamically adapt weights to the scene mean (default: false)
%   wL_eff, wM_eff, wS_eff - Fixed cone weights when not balancing neutral

arguments
    integratedRadianceMap (:,:,3) double
    options.gazeX (1,1) double = NaN
    options.gazeY (1,1) double = NaN
    options.azimuthGrid (1,:) double = -90:0.5:90
    options.elevationGrid (1,:) double = -90:0.5:90
    options.balanceNeutral (1,1) logical = false
    options.wL_eff (1,1) double = 1
    options.wM_eff (1,1) double = 1.33
    options.wS_eff (1,1) double = 5
end

persistent nativeAzimuthMap nativeElevationMap F_interp rgb2lmsMatrix

% Load the uniform transformation matrix from the derived directory
if isempty(rgb2lmsMatrix)
    matrixFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'sensorToConeMatrix.mat');
    load(matrixFileName, 'rgb2lmsMatrix');
end

% Load mapping grids if not already in memory
if isempty(nativeAzimuthMap)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraToVisualAngles.mat');
    load(mapFileName, 'nativeAzimuthMap', 'nativeElevationMap');
    
    F_interp = scatteredInterpolant(nativeAzimuthMap(:), nativeElevationMap(:), ...
        zeros(numel(nativeAzimuthMap), 1), 'linear', 'none');
end

% Get the image dimensions and handle default gaze
[H_native, W_native, ~] = size(integratedRadianceMap);
if isnan(options.gazeX)
    options.gazeX = H_native / 2;
end
if isnan(options.gazeY)
    options.gazeY = W_native / 2;
end

gx = max(1, min(W_native, options.gazeX));
gy = max(1, min(H_native, options.gazeY));

[AzGrid, ElGrid] = meshgrid(options.azimuthGrid, options.elevationGrid);
[outH, outW] = size(AzGrid);

gazeAz = interp2(nativeAzimuthMap, gx, gy, 'linear');
gazeEl = interp2(nativeElevationMap, gx, gy, 'linear');

QueryAz = AzGrid + gazeAz;
QueryEl = ElGrid + gazeEl;

% Resample camera radiance to visual angle coordinates
resampledRadiance = zeros(outH, outW, 3);
for ch = 1:3
    F_interp.Values = reshape(integratedRadianceMap(:,:,ch), [], 1);
    resampledRadiance(:,:,ch) = F_interp(QueryAz, QueryEl);
end

validCoverageMask = ~isnan(resampledRadiance(:,:,1));
resampledRadiance(isnan(resampledRadiance)) = 0;

% Apply uniform 2-degree RGB-to-LMS conversion
% Flatten to N x 3, multiply by transpose of the 3x3 matrix, and reshape
radianceFlat = reshape(resampledRadiance, [], 3);
lmsFlat = (rgb2lmsMatrix * radianceFlat')';
isomerizationMap = reshape(lmsFlat, outH, outW, 3);

% Extract the individual LMS cone class maps
L = isomerizationMap(:,:,1);
M = isomerizationMap(:,:,2);
S = isomerizationMap(:,:,3);

% Identify out-of-bounds or invalid pixels
invalidMask = isinf(L) | isinf(M) | isinf(S) | isnan(L) | isnan(M) | isnan(S) | ~validCoverageMask;

% Calculate post-receptoral excitations
if options.balanceNeutral 
    % 1. Use the image mean as a "Gray World" neutral reference
    meanL = mean(L(~invalidMask));
    meanM = mean(M(~invalidMask));
    meanS = mean(S(~invalidMask));
    
    % 2. Anchor the L-weight to 1.0 to preserve absolute physical scale
    wL_eff = 1.0; 
    
    % 3. Calculate wM so that L and M are balanced (RG = 0 at neutral)
    wM_eff = (wL_eff * meanL) / meanM;
    
    % 4. Compute the adapted Achromatic (Luminance) channel
    Ach = wL_eff .* L + wM_eff .* M;
    
    % 5. Calculate wS so that S balances Ach (BY = 0 at neutral)
    meanAch = wL_eff * meanL + wM_eff * meanM;
    wS_eff = meanAch / meanS; 
else
    % Fall back to user-provided or default fixed population weights
    wL_eff = options.wL_eff;
    wM_eff = options.wM_eff;
    wS_eff = options.wS_eff;
    Ach = wL_eff .* L + wM_eff .* M;
end

% Compute the opponent channels using the effective weights
RG = wL_eff .* L - wM_eff .* M;
BY = wS_eff .* S - Ach;

% Explicitly assign NaN to the saturated/invalid and out-of-bounds pixels
Ach(invalidMask) = NaN;
RG(invalidMask)  = NaN;
BY(invalidMask)  = NaN;

% Reconstruct the spatial map
postRecepMap = zeros(outH, outW, 3);
postRecepMap(:,:,1) = Ach;
postRecepMap(:,:,2) = RG;
postRecepMap(:,:,3) = BY;

% Package the weights into an output variable
coneWeights = [wL_eff, wM_eff, wS_eff];

end