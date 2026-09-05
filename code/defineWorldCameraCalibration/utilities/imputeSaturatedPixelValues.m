function I = imputeSaturatedPixelValues(I,minispectData)
% This function imputes absolute values for saturated (Inf) or 
% floor (0) pixels by balancing energy across the camera FOV independently 
% for the R, G, and B channels, using the reconstructed spectrum from the ASM7341.
%
% If the image contains both Inf and 0 pixels, only the Inf pixels are imputed.

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
expectedChannelRadiance = miniSpectSPD' * (cameraT ./ max(cameraT,[],2))';

% Get Bayer indices for the 2D image array to isolate the color channels
bayerPattern = "BGGR";
[rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(I, bayerPattern);

% Determine the global imputation state for the image
globalHasInf = any(isinf(I(:)));

% Loop through each color channel (R=1, G=2, B=3) to perform independent imputation
for cc = 1:3
    thisChannelIdx = rgbIdx{cc};
    
    % Extract values and solid angles for this specific color channel
    channelVals = I(thisChannelIdx);
    channelSteradians = deltaSteradians(thisChannelIdx);
    
    % Define the target mask based on the global state
    if globalHasInf
        % If there is a mixture (or only Inf), we strictly target Inf pixels.
        % Any 0 values will remain 0 and contribute 0 to the valid partition.
        targetMask = isinf(channelVals);
        imputingCeiling = true;
    else
        % If there are no Inf pixels, we target the floor (0) pixels.
        targetMask = (channelVals == 0);
        imputingCeiling = false;
    end
    
    targetSteradians = sum(channelSteradians(targetMask));
    channelTotalSteradians = sum(channelSteradians);
    
    % Proceed only if there are actually pixels to impute in this channel
    if targetSteradians > 0
        % 1. Total energy over the camera's FOV for this specific color channel
        totalCameraEnergy = expectedChannelRadiance(cc) * channelTotalSteradians;
        
        % 2. Obtain the total radiant intensity of the valid (non-target) pixels
        validPartition = sum(channelVals(~targetMask) .* channelSteradians(~targetMask));
        
        % 3. Isolate the residual energy belonging to the target pixels
        energyTarget = totalCameraEnergy - validPartition;
        
        % 4. Calculate the mean radiance for the target pixels
        targetAvgRadiance = energyTarget / targetSteradians;
        
        % 5. Bounding logic
        if any(~targetMask)
            if imputingCeiling
                % Ceiling bounding: Imputed radiance cannot be lower than the brightest valid pixel
                maxValidRadiance = max(channelVals(~targetMask));
                if targetAvgRadiance < maxValidRadiance
                    targetAvgRadiance = maxValidRadiance;
                end
            else
                % Floor bounding: Imputed radiance cannot be negative or higher than the dimmest valid pixel
                minValidRadiance = min(channelVals(~targetMask));
                if targetAvgRadiance < 0
                    targetAvgRadiance = 0;
                elseif targetAvgRadiance > minValidRadiance
                    targetAvgRadiance = minValidRadiance;
                end
            end
        end
        
        % Assign the imputed radiance back to the main image array
        targetIndices = thisChannelIdx(targetMask);
        I(targetIndices) = targetAvgRadiance;
    end
end

end