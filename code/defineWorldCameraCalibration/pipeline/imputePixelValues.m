function I = imputePixelValues(I,minispectData)
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

persistent cameraT cameraWls cameraS
if isempty(cameraT)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(paramFileName,'T');
    cameraT = table2array(T(:,["red" "green" "blue"]))';
    cameraWls = T.wls;
    cameraS = WlsToS(cameraWls);
end

% Derive an estimate of the environmental SPD from the minispect, and
% spline this to match the cameraT
[miniSpectSPD,miniSpectS] = estimateRadianceSpectrumFromMinispect(minispectData.AS);
miniSpectSPD = SplineRaw(SToWls(miniSpectS),miniSpectSPD,cameraWls);

% Get Bayer indices for the 2D image array to isolate the color channels
bayerPattern = "BGGR";
[rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(I, bayerPattern);

% Determine the global imputation state for the image
globalHasInf = any(isinf(I(:)));

% Loop through each color channel (R=1, G=2, B=3) to perform independent imputation
for cc = 1:3

    % Extract values and solid angles for this specific color channel
    thisChannelIdx = rgbIdx{cc};
    channelVals = I(thisChannelIdx);
    channelSteradians = deltaSteradians(thisChannelIdx);

    % Define the target mask based on the global state
    if globalHasInf
        % If there is a mixture (or only Inf), we strictly target Inf pixels.
        % Any 0 values will remain 0 and contribute 0 to the valid partition.
        imputeMask = isinf(channelVals);
    else
        % If there are no Inf pixels, we target the floor (0) pixels.
        imputeMask = (channelVals == 0);
    end

    % What are the steradians of the to-be-imputed pixels?
    imputeSteradians = sum(channelSteradians(imputeMask));

    % Do we have any pixels to impute?
    if imputeSteradians > 0

        % The normed spectral sensitivity for this IMX219 chanel
        thisSensitivity = cameraT(cc,:);
        thisSensitivityNormed = thisSensitivity ./ max(thisSensitivity);

        % The total energy expected for this IMX219 channel based upon the mean
        % radiance spectrum as observed by the minispect
        channelRadianceEst = (thisSensitivityNormed * miniSpectSPD) * cameraS(2);
        channelSolidAngleSum = sum(channelSteradians);
        totalChannelEnergy = channelRadianceEst * channelSolidAngleSum;

        % The energy present in the non-saturated pixels
        unsaturatedChannelEnergy = sum(channelVals(~imputeMask) .* channelSteradians(~imputeMask));

        % Distribute the remaining energy amongst the imputable pixels
        imputableChannelEnergy = totalChannelEnergy - unsaturatedChannelEnergy;
        imputableChannelEnergyPerPixel = imputableChannelEnergy / imputeSteradians;
        I(thisChannelIdx(imputeMask)) = imputableChannelEnergyPerPixel;

    end
end

end