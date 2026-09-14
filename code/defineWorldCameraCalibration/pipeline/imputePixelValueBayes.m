function imgFixed = imputePixelValueBayes(img)
% IMPUTEPIXELVALUEBAYES Imputes saturated pixels in an RGB image 
% following the Zhang & Brainard (2004) Bayesian algorithm.
%
% Inputs:
%   img - Matrix of dimension [H, W, 3] where saturated values are Inf.
%
% Outputs:
%   imgFixed - Matrix of dimension [H, W, 3] with saturated values imputed.

[H, W, C] = size(img);
if C ~= 3
    error('Input image must have 3 color channels (RGB).');
end

% Reshape image to [H*W, C] for matrix operations
pixels = reshape(img, [], C);

% Identify pixels where any channel is saturated (Inf)
isSaturatedAny = any(isinf(pixels), 2);

% Extract nonsaturated pixels to estimate the prior distribution
nonsatPixels = pixels(~isSaturatedAny, :);

if isempty(nonsatPixels)
    error('No nonsaturated pixels found to estimate the prior distribution.');
end

% Estimate prior mean (mu) and covariance (S) from nonsaturated data
mu = mean(nonsatPixels, 1)'; % [C x 1]
S = cov(nonsatPixels);       % [C x C]

% Determine the saturation level s for each channel
s = zeros(C, 1);
for c = 1:C
    finiteVals = pixels(~isinf(pixels(:, c)), c);
    if isempty(finiteVals)
        s(c) = 1.0; % Fallback if an entire channel is saturated
    else
        s(c) = max(finiteVals);
    end
end

% Rank color channels based on distance between prior mean and saturation level
variances = diag(S);
d = (s - mu) ./ sqrt(variances);
[~, channelOrder] = sort(d);

fixedPixels = pixels;

% --- Sequential Procedure Over Channels ---
for idx = 1:C
    cTarget = channelOrder(idx);
    
    % Find pixels where the target channel is saturated
    satMask = isinf(fixedPixels(:, cTarget));
    if ~any(satMask)
        continue;
    end
    
    otherChannels = setdiff(1:C, cTarget);
    
    muS = mu(cTarget);
    muK = mu(otherChannels);
    
    SS = S(cTarget, cTarget);
    SK = S(otherChannels, otherChannels);
    SSK = S(cTarget, otherChannels); % [1 x (C-1)]
    
    % Regularize SK slightly for numerical stability
    SKInv = inv(SK + 1e-6 * eye(length(otherChannels)));
    
    sVal = s(cTarget);
    satIndices = find(satMask);
    
    % Impute each saturated pixel in the target channel
    for pIdx = satIndices'
        k = fixedPixels(pIdx, otherChannels)'; % Given values of other channels
        
        % Compute conditional mean and variance
        muXs = muS + SSK * SKInv * (k - muK);
        Sxs = SS - SSK * SKInv * SSK';
        
        Sxs = max(Sxs, 1e-8); % Ensure positive variance
        stdXs = sqrt(Sxs);
        
        % Truncated normal expected value formula (Eq. 11 in paper)
        zScore = (sVal - muXs) / stdXs;
        Z = 1 - normcdf(zScore);
        
        if Z < 1e-15
            expectedVal = sVal; % Fallback to threshold if probability is near zero
        else
            numeratorTerm = (stdXs / sqrt(2 * pi)) * exp(-(zScore^2) / 2);
            expectedVal = muXs + numeratorTerm / Z;
        end
        
        fixedPixels(pIdx, cTarget) = expectedVal;
    end
end

% Reshape back to original image dimensions [H, W, 3]
imgFixed = reshape(fixedPixels, [H, W, C]);

end