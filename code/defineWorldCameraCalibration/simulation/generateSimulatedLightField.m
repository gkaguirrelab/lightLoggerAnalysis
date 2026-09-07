function [radianceModel, radianceModelS ] = generateSimulatedLightField(wls)
% GENERATESIMULATEDLIGHTFIELD Creates a continuous hemisphere radiance model.
%   [radianceModel, debugVisuals] = generateSimulatedLightField(wls) 
%   returns a function handle 'radianceModel' that computes spectral radiance 
%   (W/m2/sr/nm) for arbitrary azimuth and elevation coordinates. 

    if nargin < 1
        wls = 380:1:780; % Default to visual spectrum with 2nm bins
    end
    wls = wls(:);

    radianceModelS = WlsToS(wls);
    
    %% 1. Define Spectral Bounds (Daylight Proxy)
    c1 = 3.74177e-16; 
    c2 = 1.4388e-2;   
    wls_m = wls * 1e-9;
    
    calcBlackbody = @(T) (c1 ./ (wls_m.^5)) ./ (exp(c2 ./ (wls_m .* T)) - 1);
    
    spdCool = calcBlackbody(6500);
    spdCool = spdCool ./ sum(spdCool); 
    
    spdWarm = calcBlackbody(5000);
    spdWarm = spdWarm ./ sum(spdWarm);
    
    %% 2. Define Spatial Background Variation
    intensityMap = @(az, el) 1 + (1/3) * sin(2*az) .* cos(el); 
    spectralMix  = @(az, el) 0.5 + 0.5 * sin(az - 2*el);
    
    %% 3. Inject 20 Random Hotspots across Full Elevation Range
    nHotspots = 20;
    hAz = (rand(nHotspots, 1) * 2 * pi) - pi;         % -pi to pi
    hEl = (rand(nHotspots, 1) * pi) - (pi/2);         % -pi/2 to pi/2
    hSigma = 0.08;                            % ~4.5 degrees spatial width
    
    hotspotMultiplier = @(az, el) 1 + 4 * computeHotspotField(az, el, hAz, hEl, hSigma);
    
    %% 4. Numerically Balance Mean Radiance to 5 W/m2/sr
    gridRes = 500;
    azGrid = linspace(-pi, pi, gridRes*2);
    elGrid = linspace(-pi/2, pi/2, gridRes);          % Full elevation range
    [AZ, EL] = meshgrid(azGrid, elGrid);
    
    relIntensity = intensityMap(AZ, EL) .* hotspotMultiplier(AZ, EL);
    solidAngleWeights = cos(EL);
    meanRelativeRadiance = sum(relIntensity(:) .* solidAngleWeights(:)) / sum(solidAngleWeights(:));
    
    targetMean = 5.0; 
    nmStep = mean(diff(wls));
    scaleFactor = (targetMean / meanRelativeRadiance) / nmStep;

    %% 5. Construct the Final Radiance Function Handle
    radianceModel = @(az, el) computeSpectralRadiance(...
        az, el, wls, spdCool, spdWarm, intensityMap, spectralMix, hotspotMultiplier, scaleFactor);
    
    %% 6. Export Diagnostics
    debugVisuals.AZ = AZ;
    debugVisuals.EL = EL;
    debugVisuals.BroadbandRadiance = relIntensity * scaleFactor * nmStep; 
end

%% Local Functions

function hField = computeHotspotField(az, el, hAz, hEl, sigma)
    hField = zeros(size(az));
    for i = 1:length(hAz)
        distSq = (az - hAz(i)).^2 + (el - hEl(i)).^2;
        hField = hField + exp(-distSq / (2 * sigma^2));
    end
    hField = min(hField, 1); 
end

function spdOut = computeSpectralRadiance(az, el, wls, spdCool, spdWarm, intMap, mixMap, hsMap, scale)

    numWls = length(wls);
    
    azVec = az(:)';
    elVec = el(:)';
    
    I = intMap(azVec, elVec);
    M = mixMap(azVec, elVec);
    H = hsMap(azVec, elVec);
    
    blendedSpectra = spdCool * M + spdWarm * (1 - M);
    spdOut2D = blendedSpectra .* (I .* H .* scale);
    
    if isscalar(az)
        spdOut = spdOut2D;
    else
        outSize = [numWls, size(az)];
        spdOut = reshape(spdOut2D, outSize);
    end
end