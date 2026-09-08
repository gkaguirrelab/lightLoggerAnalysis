function [radianceModel, radianceModelS ] = generateUniformLightField(wls)
% GENERATESIMULATEDLIGHTFIELD Creates a uniform continuous hemisphere radiance model.
%   [radianceModel, debugVisuals] = generateSimulatedLightField(wls) 
%   returns a function handle 'radianceModel' that computes uniform spectral radiance 
%   (W/m2/sr/nm) for arbitrary azimuth and elevation coordinates. 

    if nargin < 1
        wls = 380:1:980; % Default to visual spectrum with 2nm bins
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
    
    % Define a uniform spectral mix across the hemisphere (e.g., 50/50 blend)
    blendedSpectrum = spdCool * 0.5 + spdWarm * 0.5;
    
    %% 2. Numerically Balance Mean Radiance to 5 W/m2/sr uniformly
    targetMean = 5.0; 
    nmStep = mean(diff(wls));
    
    % For a uniform field, scaling factor directly maps the unit blend to the target mean
    scaleFactor = targetMean / (sum(blendedSpectrum) * nmStep);

    %% 3. Construct the Final Radiance Function Handle
    radianceModel = @(az, el) computeUniformSpectralRadiance(...
        az, el, wls, blendedSpectrum, scaleFactor);
    
    %% 4. Export Diagnostics
    gridRes = 500;
    azGrid = linspace(-pi, pi, gridRes*2);
    elGrid = linspace(-pi/2, pi/2, gridRes);          
    [AZ, EL] = meshgrid(azGrid, elGrid);
    
    debugVisuals.AZ = AZ;
    debugVisuals.EL = EL;
    debugVisuals.BroadbandRadiance = ones(size(AZ)) * targetMean; 
end

%% Local Functions

function spdOut = computeUniformSpectralRadiance(az, el, wls, blendedSpectrum, scale)
    numWls = length(wls);
    
    azVec = az(:)';
    elVec = el(:)';
    
    % Replicate the uniform spectrum across all spatial points
    spdOut2D = (blendedSpectrum * scale) * ones(size(azVec));
    
    if isscalar(az)
        spdOut = spdOut2D(:, 1);
    else
        outSize = [numWls, size(az)];
        spdOut = reshape(spdOut2D, outSize);
    end
end