function I = imputePixelValues(I, minispectData, equalizeChannels)
% This function imputes absolute values for saturated (Inf) or
% floor (0) pixels by balancing energy across the camera FOV,
% using the reconstructed spectrum from the ASM7341.
%
% If the image contains both Inf and 0 pixels, only the Inf pixels are imputed.
%
% Inputs:
%   I                 - 2D image array
%   minispectData     - structure containing ASM7341 measurements
%   equalizeChannels  - logical flag (optional, default = true). 
%                       If true, pools the total imputable energy across all 
%                       channels and assigns a uniform radiance value to all 
%                       saturated pixels. If false, imputes independently per channel.

if nargin < 3 || isempty(equalizeChannels)
    equalizeChannels = true;
end

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

% Preallocate storage for channel metrics
totalChannelEnergy = zeros(1, 3);
unsaturatedChannelEnergy = nan(1, 3);
imputableChannelEnergy = zeros(1, 3);
imputeSteradiansArr = zeros(1, 3);
imputeMasks = cell(1, 3);

% Loop through each color channel (R=1, G=2, B=3) to collect metrics
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
    imputeMasks{cc} = imputeMask;

    % What are the steradians of the to-be-imputed pixels?
    imputeSteradians = sum(channelSteradians(imputeMask));
    imputeSteradiansArr(cc) = imputeSteradians;

    % Do we have any pixels to impute in this channel?
    if imputeSteradians > 0

        % The normed spectral sensitivity for this IMX219 channel
        thisSensitivity = cameraT(cc,:);
        thisSensitivityNormed = thisSensitivity ./ max(thisSensitivity);

        % The total energy expected for this IMX219 channel based upon the mean
        % radiance spectrum as observed by the minispect
        channelRadianceEst = (thisSensitivityNormed * miniSpectSPD) * cameraS(2);
        channelSolidAngleSum = sum(channelSteradians);
        totalChannelEnergy(cc) = channelRadianceEst * channelSolidAngleSum;

        % The energy present in the non-saturated pixels
        unsaturatedChannelEnergy(cc) = sum(channelVals(~imputeMask) .* channelSteradians(~imputeMask));

        % The energy available for imputation in this channel
        imputableChannelEnergy(cc) = totalChannelEnergy(cc) - unsaturatedChannelEnergy(cc);

    end
end

if equalizeChannels
    % Estimate overall imputable energy and divide equally across all saturated pixels
    totalExpectedEnergy = sum(totalChannelEnergy);
    totalUnsaturatedEnergy = sum(unsaturatedChannelEnergy, 'omitnan');
    totalImputableEnergy = totalExpectedEnergy - totalUnsaturatedEnergy;
    totalImputableSteradians = sum(imputeSteradiansArr);

    if totalImputableSteradians > 0
        uniformImputableEnergyPerPixel = totalImputableEnergy / totalImputableSteradians;
        for cc = 1:3
            if imputeSteradiansArr(cc) > 0
                imputeMask = imputeMasks{cc};
                I(rgbIdx{cc}(imputeMask)) = uniformImputableEnergyPerPixel;
            end
        end
    end
else
    % Impute each channel independently (original behavior)
    for cc = 1:3
        imputeSteradians = imputeSteradiansArr(cc);
        if imputeSteradians > 0
            imputeMask = imputeMasks{cc};
            imputableChannelEnergyPerPixel = imputableChannelEnergy(cc) / imputeSteradians;
            I(rgbIdx{cc}(imputeMask)) = imputableChannelEnergyPerPixel;
        end
    end
end

fprintf('Total unsaturated camera energy = %2.2f\n', sum(unsaturatedChannelEnergy, 'omitnan'));

end