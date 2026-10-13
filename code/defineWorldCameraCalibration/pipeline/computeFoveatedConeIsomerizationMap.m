function [isomerizationMap, AzGrid, ElGrid] = computeFoveatedConeIsomerizationMap(integratedRadianceMap, options)
% Converts an RGB integrated radiance map into an LMS cone isomerization
% rate map on an evenly sampled visual angle grid. The output units are
% "pooled" isomerization rates, in units of R*/deg^2/second.
%
% Pre-receptoral optics are dynamically applied by convolving the maps with 
% pre-computed point spread functions loaded from the derived support file.

arguments
    integratedRadianceMap (:,:,3) double
    options.gazeX (1,1) double = NaN
    options.gazeY (1,1) double = NaN
    options.observerAge (1,1) double = 25
    options.pupilDiameterMm (1,1) double = 3.0
    options.azimuthGrid (1,:) double = -90:0.5:90
    options.elevationGrid (1,:) double = -90:0.5:90
    options.fovealTritanopiaFlag (1,1) logical = false
end

persistent nativeAzimuthMap nativeElevationMap coneMapVar opticsSupport F_interp

% Load mapping grids
if isempty(nativeAzimuthMap)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraToVisualAngles.mat');
    load(mapFileName, 'nativeAzimuthMap', 'nativeElevationMap');
    
    F_interp = scatteredInterpolant(nativeAzimuthMap(:), nativeElevationMap(:), ...
        zeros(numel(nativeAzimuthMap), 1), 'linear', 'none');
end

% Load both the cone and optics support variables
if isempty(coneMapVar) || isempty(opticsSupport)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radianceToConeRateSupport.mat');
    load(mapFileName, 'coneMapVar', 'opticsSupport');
end

% Match pupil size to the closest precomputed optics struct
pupilDiffs = abs([opticsSupport.pupilDiameterMm] - options.pupilDiameterMm);
[~, bestPupilIdx] = min(pupilDiffs);
matchedOptics = opticsSupport(bestPupilIdx);

% Get the age-specific coneMapVar
thisConeMapVar = coneMapVar(options.observerAge);
if isempty(thisConeMapVar.observerAge)
    error('computeFoveatedConeIsomerizationMap: Subject age not defined for the coneMapVar');
end

% Get the image dimensions
[H_native, W_native, ~] = size(integratedRadianceMap);

if isnan(options.gazeX)
    options.gazeX = H_native/2;
end
if isnan(options.gazeY)
    options.gazeY = W_native/2;
end

% Pupil scalar relative to mapping base
pupilScalar = (options.pupilDiameterMm / thisConeMapVar.basePupilDiameterMm)^2;

gx = max(1, min(W_native, options.gazeX));
gy = max(1, min(H_native, options.gazeY));

[AzGrid, ElGrid] = meshgrid(options.azimuthGrid, options.elevationGrid);
[outH, outW] = size(AzGrid);

% Calculate target spatial resolution in arcmin for PSF resampling
gridSpacingDeg = abs(options.azimuthGrid(2) - options.azimuthGrid(1));
gridSpacingArcMin = gridSpacingDeg * 60;

gazeAz = interp2(nativeAzimuthMap, gx, gy, 'linear');
gazeEl = interp2(nativeElevationMap, gx, gy, 'linear');

QueryAz = AzGrid + gazeAz;
QueryEl = ElGrid + gazeEl;

resampledRadiance = zeros(outH, outW, 3);
for ch = 1:3
    F_interp.Values = reshape(integratedRadianceMap(:,:,ch), [], 1);
    resampledRadiance(:,:,ch) = F_interp(QueryAz, QueryEl);
end

validCoverageMask = ~isnan(resampledRadiance(:,:,1));
resampledRadiance(isnan(resampledRadiance)) = 0;

gridUnitDirs_1 = cosd(QueryEl) .* sind(QueryAz);
gridUnitDirs_2 = -sind(QueryEl);
gridUnitDirs_3 = cosd(QueryEl) .* cosd(QueryAz);

fixVec_1 = cosd(gazeEl) * sind(gazeAz);
fixVec_2 = -sind(gazeEl);
fixVec_3 = cosd(gazeEl) * cosd(gazeAz);

dotProducts = gridUnitDirs_1 .* fixVec_1 + gridUnitDirs_2 .* fixVec_2 + gridUnitDirs_3 .* fixVec_3;
dotProducts = min(max(dotProducts, -1), 1);
dynamicEccMap = rad2deg(acos(dotProducts));

T_flat = interp1(thisConeMapVar.eccGrid, thisConeMapVar.transformTable, dynamicEccMap(:), 'linear', 'extrap');
T_map = reshape(T_flat, outH, outW, 3, 3);

isomerizationMap = zeros(outH, outW, 3);

for coneClass = 1:3 % 1=L, 2=M, 3=S
    transformWeights = squeeze(T_map(:, :, coneClass, :));
    isoChannel = sum(resampledRadiance .* transformWeights, 3) .* pupilScalar;
    
    % Extract the appropriate precomputed PSF and its native spacing
    psf = matchedOptics.psf{coneClass};
    psfSpacing = matchedOptics.psfSpacingArcMin(coneClass);
    
    % Resize the PSF kernel to match the actual grid spacing
    resizeFactor = psfSpacing / gridSpacingArcMin;
    if resizeFactor ~= 1
        psfResampled = imresize(psf, resizeFactor, 'bilinear');
        psfResampled = max(psfResampled, 0); 
        psfResampled = psfResampled / sum(psfResampled(:)); 
    else
        psfResampled = psf;
    end
    
    % Convolve the isomerization map with the scaled PSF
    isoChannelBlurred = imfilter(isoChannel, psfResampled, 'replicate', 'same');
    isomerizationMap(:,:,coneClass) = isoChannelBlurred;
end

isomerizationMap(~repmat(validCoverageMask, 1, 1, 3)) = NaN;

if options.fovealTritanopiaFlag
    sConeMask = dynamicEccMap >= 0.175;
    isomerizationMap(:,:,3) = isomerizationMap(:,:,3) .* sConeMask;
end

densityFlat = interp1(thisConeMapVar.eccGrid, thisConeMapVar.densityTable, dynamicEccMap(:), 'linear', 'extrap');
densityMap = reshape(densityFlat, outH, outW, 3);

isomerizationMap = isomerizationMap .* densityMap;

end