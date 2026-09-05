function I = imputeSaturatedPixelValues(I,minispectData)
% This function imputes absolute values for saturated pixels
% by balancing energy across the camera FOV independently for the 
% R, G, and B channels, using the reconstructed spectrum from the ASM7341.

persistent deltaSteradians
if isempty(deltaSteradians)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'deltaSteradians.mat');
    load(paramFileName,'deltaSteradians');
end

persistent cameraT cameraWls
if isempty(cameraT)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(paramFileName,'T');
    cameraT = table2array(T(:,["red" "green" "blue"]))';
    cameraWls = T.wls;
end

% Derive an estimate of the environmental SPD from the minispect
[miniSpectSPD,miniSpectS] = estimateRadianceSpectrumFromMinispect(minispectData.AS(1:9));
miniSpectSPD = SplineRaw(SToWls(miniSpectS),miniSpectSPD,cameraWls);
miniSpectS = WlsToS(cameraWls);

% Obtain the expected mean scene radiance for each channel (R, G, B) independently.
% This yields a 1x3 vector of predicted radiances based on the spectral shape.
expectedChannelRadiance = miniSpectSPD' * (cameraT ./ max(cameraT,[],2))';

% Get Bayer indices for the 2D image array to isolate the color channels
bayerPattern = "BGGR";
[rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(I, bayerPattern);

% Loop through each color channel (R=1, G=2, B=3) to perform independent imputation
for cc = 1:3
    thisChannelIdx = rgbIdx{cc};
    
    % Extract values and solid angles for this specific color channel
    channelVals = I(thisChannelIdx);
    channelSteradians = deltaSteradians(thisChannelIdx);
    
    % Identify saturated pixels strictly within this channel
    infMask = isinf(channelVals);
    
    channelSatSteradians = sum(channelSteradians(infMask));
    channelTotalSteradians = sum(channelSteradians);
    
    % Obtain the total radiant intensity of the unsaturated pixels for this channel
    unsatPartition = sum(channelVals(~infMask) .* channelSteradians(~infMask));
    
    % 1. Total energy over the camera's FOV for this specific color channel
    totalCameraEnergy = expectedChannelRadiance(cc) * channelTotalSteradians;
    
    % 2. Isolate the residual energy belonging to the saturated pixels in this channel
    energySaturated = totalCameraEnergy - unsatPartition;
    
    % 3. Calculate the mean radiance and assign it if saturated pixels exist
    if channelSatSteradians > 0
        cameraSaturatedAvgRadiance = energySaturated / channelSatSteradians;
        
        % 4. Bounding: Ensure imputed radiance isn't lower than the brightest unsaturated pixel
        if any(~infMask)
            maxUnsatRadiance = max(channelVals(~infMask));
            if cameraSaturatedAvgRadiance < maxUnsatRadiance
                cameraSaturatedAvgRadiance = maxUnsatRadiance;
            end
        end
        
        % Assign the imputed saturated radiance back to the main image array
        satIndices = thisChannelIdx(infMask);
        I(satIndices) = cameraSaturatedAvgRadiance;
    end
end

end