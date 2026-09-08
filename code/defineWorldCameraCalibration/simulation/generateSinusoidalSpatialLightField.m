function [radianceModel, radianceModelS] = generateSinusoidalSpatialLightField(wls)
% GENERATESINUSOIDALSPATIALLIGHTFIELD Creates a continuous hemisphere radiance model
% with a high spatial frequency horizontal and vertical sinusoidal variation, 
% ensuring center/periphery mean balance while introducing fine-scale spatial structure.

    if nargin < 1
        wls = 380:1:980; % Default to visual spectrum with 2nm bins
    end
    wls = wls(:);

    radianceModelS = WlsToS(wls);
    
    %% 1. Define Spectral Bounds (Daylight Proxy - matching generateSimulatedLightField)
    c1 = 3.74177e-16; 
    c2 = 1.4388e-2;   
    wls_m = wls * 1e-9;
    
    calcBlackbody = @(T) (c1 ./ (wls_m.^5)) ./ (exp(c2 ./ (wls_m .* T)) - 1);
    
    spdCool = calcBlackbody(6500);
    spdCool = spdCool ./ sum(spdCool); 
    
    spdWarm = calcBlackbody(5000);
    spdWarm = spdWarm ./ sum(spdWarm);
    
    %% 2. Define High-Frequency Sinusoidal Spatial Variation
    % High frequency factor ensures rapid oscillations across both azimuth and elevation,
    % so local spatial averages in the camera FOV converge to the global hemispheric mean.
    freq = 12; 
    intensityMap = @(az, el) 1 + 0.25 * sin(freq * az) .* cos(freq * el); 
    spectralMix  = @(az, el) 0.5 + 0.2 * sin(2 * az - el);
    
    %% 3. Numerically Balance Mean Radiance to 5 W/m2/sr
    gridRes = 500;
    azGrid = linspace(-pi, pi, gridRes*2);
    elGrid = linspace(-pi/2, pi/2, gridRes);          
    [AZ, EL] = meshgrid(azGrid, elGrid);
    
    relIntensity = intensityMap(AZ, EL);
    solidAngleWeights = cos(EL);
    meanRelativeRadiance = sum(relIntensity(:) .* solidAngleWeights(:)) / sum(solidAngleWeights(:));
    
    targetMean = 5.0; 
    nmStep = mean(diff(wls));
    scaleFactor = (targetMean / meanRelativeRadiance) / nmStep;

    %% 4. Construct the Final Radiance Function Handle
    radianceModel = @(az, el) computeSinusoidalSpectralRadiance(...
        az, el, wls, spdCool, spdWarm, intensityMap, spectralMix, scaleFactor);
end

%% Local Functions

function spdOut = computeSinusoidalSpectralRadiance(az, el, wls, spdCool, spdWarm, intMap, mixMap, scale)

    numWls = length(wls);
    
    azVec = az(:)';
    elVec = el(:)';
    
    I = intMap(azVec, elVec);
    M = mixMap(azVec, elVec);
    
    blendedSpectra = spdCool * M + spdWarm * (1 - M);
    spdOut2D = blendedSpectra .* (I .* scale);
    
    if isscalar(az)
        spdOut = spdOut2D;
    else
        outSize = [numWls, size(az)];
        spdOut = reshape(spdOut2D, outSize);
    end
end