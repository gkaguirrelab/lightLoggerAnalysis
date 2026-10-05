function [isomerizationMap, AzGrid, ElGrid] = computeFoveatedConeIsomerizationMap(integratedRadianceMap, options)
% Converts an RGB integrated radiance map into an LMS cone isomerization
% rate map on an evenly sampled visual angle grid. The output units are
% "pooled" isomerization rates, in units of R*/deg^2/second.
%
% The output is also subject virtual foveation by taking the gaze location 
% as pixel coordinates (gazeX, gazeY) within the IMX camera image. The 
% native radiance is reprojected onto a fovea-centric visual angle grid, 
% assigning the designated gaze [0,0] azimuth and elevation.

arguments
    integratedRadianceMap (:,:,3) double
    options.gazeX (1,1) double = NaN
    options.gazeY (1,1) double = NaN
    options.pupilDiameterMm (1,1) double = 3.0
    options.azimuthGrid (1,:) double = -90:0.5:90
    options.elevationGrid (1,:) double = -90:0.5:90
    options.fovealTritanopiaFlag (1,1) logical = false;
end

persistent nativeAzimuthMap nativeElevationMap coneMapVar F_interp

% Load only the specific visual angle maps needed for geometry
if isempty(nativeAzimuthMap)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraToVisualAngles.mat');
    load(mapFileName, 'nativeAzimuthMap', 'nativeElevationMap');
    
    % Initialize the scatteredInterpolant object once to compute the 
    % Delaunay triangulation. This will persist across image frames.
    F_interp = scatteredInterpolant(nativeAzimuthMap(:), nativeElevationMap(:), ...
        zeros(numel(nativeAzimuthMap), 1), 'linear', 'none');
end

if isempty(coneMapVar)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radianceToConeRateSupport.mat');
    load(mapFileName, 'coneMapVar');
end

% Get the image dimensions
[H_native, W_native, ~] = size(integratedRadianceMap);

% Default to the image center if gaze coordinates are not provided
if isnan(options.gazeX)
    options.gazeX = H_native/2;
end
if isnan(options.gazeY)
    options.gazeY = W_native/2;
end

% Establish the pupil area scalar relative to the canonical mapping base
pupilScalar = (options.pupilDiameterMm / coneMapVar.basePupilDiameterMm)^2;

% Constrain the gaze coordinates within the image bounds
gx = max(1, min(W_native, options.gazeX));
gy = max(1, min(H_native, options.gazeY));

% Create the returned evenly sampled visual angle grid (fovea-centric)
[AzGrid, ElGrid] = meshgrid(options.azimuthGrid, options.elevationGrid);
[outH, outW] = size(AzGrid);

% Find the visual angle of the gaze center in absolute camera coordinates
gazeAz = interp2(nativeAzimuthMap, gx, gy, 'linear');
gazeEl = interp2(nativeElevationMap, gx, gy, 'linear');

% Create internal query grids shifted by the gaze position to sample the
% absolute camera space
QueryAz = AzGrid + gazeAz;
QueryEl = ElGrid + gazeEl;

% Reproject the native integrated radiance map onto the regular grid
resampledRadiance = zeros(outH, outW, 3);
for ch = 1:3
    % Update only the values array rather than recreating the interpolant
    F_interp.Values = reshape(integratedRadianceMap(:,:,ch), [], 1);
    
    % Sample using the shifted query coordinates
    resampledRadiance(:,:,ch) = F_interp(QueryAz, QueryEl);
end

% Mask for valid areas inside the fisheye projection
validCoverageMask = ~isnan(resampledRadiance(:,:,1));

% Temporarily set NaNs to 0 to prevent propagation errors during matrix math
resampledRadiance(isnan(resampledRadiance)) = 0;

% Convert the ABSOLUTE query angles into 3D unit vectors to reflect camera
% geometry
gridUnitDirs_1 = cosd(QueryEl) .* sind(QueryAz);
gridUnitDirs_2 = -sind(QueryEl);
gridUnitDirs_3 = cosd(QueryEl) .* cosd(QueryAz);

% Convert the absolute gaze visual angle into a 3D fixation unit vector
fixVec_1 = cosd(gazeEl) * sind(gazeAz);
fixVec_2 = -sind(gazeEl);
fixVec_3 = cosd(gazeEl) * cosd(gazeAz);

% Calculate the dynamic eccentricity map on the evenly sampled grid
dotProducts = gridUnitDirs_1 .* fixVec_1 + gridUnitDirs_2 .* fixVec_2 + gridUnitDirs_3 .* fixVec_3;
dotProducts = min(max(dotProducts, -1), 1); % Clamp to prevent acos precision errors
dynamicEccMap = rad2deg(acos(dotProducts));

% Interpolate the 3x3 matrices from the 1D LUT. This maps the dynamic
% eccentricities into an (outH*outW) x 9 matrix
T_flat = interp1(coneMapVar.eccGrid, coneMapVar.transformTable, dynamicEccMap(:), 'linear', 'extrap');

% Reshape back into the spatial matrix format matching the new grid
T_map = reshape(T_flat, outH, outW, 3, 3);

% Perform spatially varying matrix multiplication
isomerizationMap = zeros(outH, outW, 3);

for coneClass = 1:3 % 1=L, 2=M, 3=S
    % Extract the spatially varying 1x3 vector for this cone class
    transformWeights = squeeze(T_map(:, :, coneClass, :));
    
    % Multiply weights against resampled radiance map, sum across channels, and apply pupil scalar
    isoChannel = sum(resampledRadiance .* transformWeights, 3) .* pupilScalar;
    isomerizationMap(:,:,coneClass) = isoChannel;
end

% Place NaNs outside the fisheye boundary
isomerizationMap(~repmat(validCoverageMask, 1, 1, 3)) = NaN;

% Apply structural spatial constraint: Foveal Tritanopia (S-cone free zone)
% S-cones are biologically absent in the central foveola (approx. 0.35
% degree diameter) This automatically shifts with the dynamicEccMap to stay
% centered on the gaze.
if options.fovealTritanopiaFlag
sConeMask = dynamicEccMap >= 0.175;
isomerizationMap(:,:,3) = isomerizationMap(:,:,3) .* sConeMask;
end

% Interpolate the spatial cone density map (cones/degree^2)
densityFlat = interp1(coneMapVar.eccGrid, coneMapVar.densityTable, dynamicEccMap(:), 'linear', 'extrap');
densityMap = reshape(densityFlat, outH, outW, 3);

% Multiply: (R*/cone/sec) * (cones/deg^2) = (R*/deg^2/sec)
isomerizationMap = isomerizationMap .* densityMap;

end