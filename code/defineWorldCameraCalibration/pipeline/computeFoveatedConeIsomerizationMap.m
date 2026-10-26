function [isomerizationMap, AzGrid, ElGrid] = computeFoveatedConeIsomerizationMap(integratedRadianceMap, options)
% Converts an RGB integrated radiance map into an LMS cone isomerization
% rate map on an evenly sampled visual angle grid. The output units are
% "pooled" isomerization rates, in units of R*/deg^2/second.
%
% Pre-receptoral optics are dynamically applied by convolving the maps with 
% spatially varying, pre-computed point spread functions loaded from the 
% derived support file.

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

persistent nativeAzimuthMap nativeElevationMap coneMapVar opticsSupport F_interp F_transform F_density lastAge

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
        'radianceToConeSupport.mat');
    load(mapFileName, 'coneMapVar', 'opticsSupport');
end

% Match pupil size to the closest precomputed optics struct
pupilDiffs = abs(opticsSupport.pupilDiameterMmRange - options.pupilDiameterMm);
[~, bestPupilIdx] = min(pupilDiffs);

% Get the age-specific coneMapVar
thisConeMapVar = coneMapVar(options.observerAge);
if isempty(thisConeMapVar.observerAge)
    error('computeFoveatedConeIsomerizationMap: Subject age not defined for the coneMapVar');
end

% Setup persistent 2D interpolants if needed (avoids rebuilding on every frame)
if isempty(lastAge) || lastAge ~= options.observerAge
    eccX = thisConeMapVar.eccGrid(1, :)';
    eccY = thisConeMapVar.eccGrid(2, :)';
    
    F_transform = cell(1, 9);
    for k = 1:9
        F_transform{k} = scatteredInterpolant(eccX, eccY, thisConeMapVar.transformTable(:, k), 'linear', 'nearest');
    end
    
    F_density = cell(1, 3);
    for k = 1:3
        F_density{k} = scatteredInterpolant(eccX, eccY, thisConeMapVar.densityTable(:, k), 'linear', 'nearest');
    end
    lastAge = options.observerAge;
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

% Great circle distance for dynamic eccentricity map (used for tritanopia mask)
gridUnitDirs_1 = cosd(QueryEl) .* sind(QueryAz);
gridUnitDirs_2 = -sind(QueryEl);
gridUnitDirs_3 = cosd(QueryEl) .* cosd(QueryAz);

fixVec_1 = cosd(gazeEl) * sind(gazeAz);
fixVec_2 = -sind(gazeEl);
fixVec_3 = cosd(gazeEl) * cosd(gazeAz);

dotProducts = gridUnitDirs_1 .* fixVec_1 + gridUnitDirs_2 .* fixVec_2 + gridUnitDirs_3 .* fixVec_3;
dotProducts = min(max(dotProducts, -1), 1);
dynamicEccMap = rad2deg(acos(dotProducts));

% Evaluate T_map and densityMap based on 2D coordinates relative to gaze.
% We wrap AzGrid in abs() because the derived maps are now defined symmetrically for X >= 0
T_map = zeros(outH, outW, 9);
for k = 1:9
    T_map(:,:,k) = F_transform{k}(abs(AzGrid), ElGrid);
end
T_map = reshape(T_map, outH, outW, 3, 3);

if any(isnan(T_map(:)))
    warning('computeFoveatedConeIsomerizationMap:NanTransform', ...
        'NaN values detected in the interpolated transform map. Scrubbing NaNs to 0.');
    T_map(isnan(T_map)) = 0;
end

densityMap = zeros(outH, outW, 3);
for k = 1:3
    densityMap(:,:,k) = F_density{k}(abs(AzGrid), ElGrid);
end

if any(isnan(densityMap(:)))
    warning('computeFoveatedConeIsomerizationMap:NanDensity', ...
        'NaN values detected in the interpolated density map. Scrubbing NaNs to 0.');
    densityMap(isnan(densityMap)) = 0;
end

% Robustly extract PSF native spacing
psfSpacing = 1.0; 
if isfield(opticsSupport, 'theBasePSF')
    if isfield(opticsSupport.theBasePSF, 'spatialSupportMicrons')
        vals = unique(opticsSupport.theBasePSF.spatialSupportMicrons(:));
        if length(vals) > 1
            psfSpacing = abs(vals(2) - vals(1)) / (300/60);
        end
    elseif isfield(opticsSupport.theBasePSF, 'supportX')
        vals = unique(opticsSupport.theBasePSF.supportX(:));
        if length(vals) > 1
            psfSpacing = abs(vals(2) - vals(1));
        end
    end
end
if isnan(psfSpacing) || psfSpacing <= 0
    psfSpacing = 1.0;
end
resizeFactor = psfSpacing / gridSpacingArcMin;

% Pre-calculate spatial blending weights for the spatially varying PSFs
numEccs = size(thisConeMapVar.eccGrid, 2);
eccWeights = zeros(outH, outW, numEccs);
for iEcc = 1:numEccs
    eccX = thisConeMapVar.eccGrid(1, iEcc);
    eccY = thisConeMapVar.eccGrid(2, iEcc);
    % Use abs(AzGrid) to symmetrically map blending distances on both sides of the vertical meridian
    distSq = (abs(AzGrid) - eccX).^2 + (ElGrid - eccY).^2;
    eccWeights(:,:,iEcc) = 1 ./ (distSq + 1e-4); % Inverse distance squared
end
% Normalize blending weights globally
weightSum = sum(eccWeights, 3);
for iEcc = 1:numEccs
    eccWeights(:,:,iEcc) = eccWeights(:,:,iEcc) ./ weightSum;
end

isomerizationMap = zeros(outH, outW, 3);

for coneClass = 1:3 % 1=L, 2=M, 3=S
    transformWeights = squeeze(T_map(:, :, coneClass, :));
    isoChannel = sum(resampledRadiance .* transformWeights, 3) .* pupilScalar;
    
    isoChannelBlurred = zeros(outH, outW);
    
    for iEcc = 1:numEccs
        eccX = thisConeMapVar.eccGrid(1, iEcc);
        eccY = thisConeMapVar.eccGrid(2, iEcc);
        
        % Extract the cone-weighted PSF localized to this eccentricity
        thePSF = opticsSupport.thePSFs{iEcc, bestPupilIdx};
        if isempty(thePSF) || ~isfield(thePSF, 'coneWeightedData') || isempty(thePSF.coneWeightedData)
            continue;
        end
        
        psf = squeeze(thePSF.coneWeightedData(coneClass, :, :));
        
        if any(isnan(psf(:)))
            warning('computeFoveatedConeIsomerizationMap:NanPSF', ...
                'NaN values detected in PSF for cone class %d at eccentricity [%.1f, %.1f]. Scrubbing NaNs to 0.', ...
                coneClass, eccX, eccY);
            psf(isnan(psf)) = 0;
        end
        
        % Resize the PSF kernel to match the user's requested visual angle grid spacing
        if resizeFactor ~= 1 && resizeFactor > 0
            psfResampled = imresize(psf, resizeFactor, 'bilinear');
            psfResampled = max(psfResampled, 0); 
        else
            psfResampled = max(psf, 0);
        end
        
        % Guard against division by zero if the downsampled PSF zeros out
        psfSum = sum(psfResampled(:), 'omitnan');
        if psfSum > 0
            psfResampled = psfResampled / psfSum;
        else
            warning('computeFoveatedConeIsomerizationMap:ZeroSumPSF', ...
                'Resampled PSF summed to zero for cone class %d at eccentricity [%.1f, %.1f]. Applying identity impulse fallback.', ...
                coneClass, eccX, eccY);
            psfResampled = zeros(size(psfResampled));
            psfResampled(ceil(end/2), ceil(end/2)) = 1; % Fallback identity impulse
        end
        
        % Convolve and blend based on distance to the eccentricity region
        blurred = imfilter(isoChannel, psfResampled, 'replicate', 'same');
        
        if any(isnan(blurred(:)))
            warning('computeFoveatedConeIsomerizationMap:NanConvolution', ...
                'NaN values detected after spatial convolution for cone class %d at eccentricity [%.1f, %.1f]. Scrubbing residual NaNs to 0.', ...
                coneClass, eccX, eccY);
            blurred(isnan(blurred)) = 0; 
        end
        
        isoChannelBlurred = isoChannelBlurred + blurred .* eccWeights(:,:,iEcc);
    end
    
    isomerizationMap(:,:,coneClass) = isoChannelBlurred;
end

isomerizationMap(~repmat(validCoverageMask, 1, 1, 3)) = NaN;

if options.fovealTritanopiaFlag
    sConeMask = dynamicEccMap >= 0.175;
    isomerizationMap(:,:,3) = isomerizationMap(:,:,3) .* sConeMask;
end

isomerizationMap = isomerizationMap .* densityMap;

end