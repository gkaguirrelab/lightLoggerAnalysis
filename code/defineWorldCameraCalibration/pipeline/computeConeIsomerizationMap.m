function isomerizationMap = computeConeIsomerizationMap(integratedRadianceMap, coneMapVar, options)
% COMPUTECONEISOMERIZATIONMAP Converts an RGB integrated radiance map into
% an LMS cone isomerization rate map, allowing dynamic fixation and pupil
% adjustments.
%
% Inputs:
%   integratedRadianceMap - H x W x 3 RGB matrix in W/m^2/sr
%   coneMapVar            - Struct containing the 1D LUT transform table
%
% Name-Value Arguments:
%   fixationAzimuth   - Azimuth offset of fixation in degrees (default: 0)
%   fixationElevation - Elevation offset of fixation in degrees (default: 0)
%   pupilDiameterMm   - Pupil diameter for this specific image in mm (default: 3.0)

arguments
    integratedRadianceMap (:,:,3) double
    coneMapVar (1,1) struct
    options.fixationAzimuth (1,1) double = 0
    options.fixationElevation (1,1) double = 0
    options.pupilDiameterMm (1,1) double = 3.0
end

persistent unitDirections

% Load the unit vectors defining the camera geometry
if isempty(unitDirections)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'deltaSteradians.mat');
    
    % Note: ensure defineDeltaSteradians saves 'unitDirections'
    load(mapFileName, 'unitDirections');
end

[H, W, ~] = size(integratedRadianceMap);

% 1. Establish the pupil area scalar relative to the canonical mapping base
pupilScalar = (options.pupilDiameterMm / coneMapVar.basePupilDiameterMm)^2;

% 2. Calculate the dynamic eccentricity map based on fixation offset
az = deg2rad(options.fixationAzimuth);
el = deg2rad(options.fixationElevation);

% Create the fixation unit vector using the established spatial convention
fixationVector = [cos(el)*sin(az), -sin(el), cos(el)*cos(az)];
fixationVector = reshape(fixationVector, 1, 1, 3);

% The dot product of the unit directions and the fixation vector yields the cosine
% of the angular distance (eccentricity) from the fovea
dotProducts = sum(unitDirections .* fixationVector, 3);
dotProducts = min(max(dotProducts, -1), 1); % Clamp to prevent acos precision errors
dynamicEccMap = rad2deg(acos(dotProducts));

% 3. Interpolate the 3x3 matrices from the 1D LUT
% This maps the HxW dynamic eccentricities into an (H*W) x 9 matrix
T_flat = interp1(coneMapVar.eccGrid, coneMapVar.transformTable, dynamicEccMap(:), 'linear', 'extrap');

% Reshape back into the H x W x 3 x 3 spatial matrix format
T_map = reshape(T_flat, H, W, 3, 3);

% 4. Perform spatially varying matrix multiplication
isomerizationMap = zeros(H, W, 3);

for coneClass = 1:3 % 1=L, 2=M, 3=S
    % Extract the spatially varying 1x3 vector for this cone class[cite: 5]
    transformWeights = squeeze(T_map(:, :, coneClass, :));

    % Multiply weights against radiance map, sum across channels, and apply pupil scalar
    isoChannel = sum(integratedRadianceMap .* transformWeights, 3) .* pupilScalar;

    isomerizationMap(:,:,coneClass) = isoChannel;
end

% 5. Threshold non-physiologic negative values[cite: 5]
isomerizationMap(isomerizationMap < 0) = 0;

end