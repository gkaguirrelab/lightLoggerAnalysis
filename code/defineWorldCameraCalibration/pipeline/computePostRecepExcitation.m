function postRecepMap = computePostRecepExcitation(isomerizationMap, options)
% COMPUTEPOSTRECEXCITATION Transforms an LMS cone isomerization rate map
% into a Linear Opponent Excitation map representing post-receptoral
% pathways.
%
% Inputs:
%   isomerizationMap - H x W x 3 matrix of LMS isomerization rates
%
% Name-Value Arguments:
%   balanceNeutral - Boolean to dynamically adapt weights to the scene mean (default: true)

arguments
    isomerizationMap (:,:,3) double
    options.balanceNeutral (1,1) logical = true
end

% Extract the individual cone class maps
L = isomerizationMap(:,:,1);
M = isomerizationMap(:,:,2);
S = isomerizationMap(:,:,3);

% Identify pixels that contain Inf (unimputed saturation) or NaN
invalidMask = isinf(L) | isinf(M) | isinf(S) | isnan(L) | isnan(M) | isnan(S);

if options.balanceNeutral
    % 1. Use the image mean as a "Gray World" neutral reference
    meanL = mean(L(~invalidMask));
    meanM = mean(M(~invalidMask));
    meanS = mean(S(~invalidMask));
    
    % 2. Anchor the L-weight to 1.0 to preserve absolute physical scale
    wL_eff = 1.0;
    
    % 3. Calculate wM so that L and M are balanced (RG = 0 at neutral)
    wM_eff = (wL_eff * meanL) / meanM;
    
    % 4. Compute the adapted Achromatic (Luminance) channel
    Ach = wL_eff .* L + wM_eff .* M;
    
    % 5. Calculate wS so that S balances Ach (BY = 0 at neutral)
    meanAch = wL_eff * meanL + wM_eff * meanM;
    wS_eff = meanAch / meanS;
else
    % Fall back to fixed population weights (will likely suffer from luminance leakage)
    wL_eff = 0.633;
    wM_eff = 0.317;
    wS_eff = 0.050;
    Ach = wL_eff .* L + wM_eff .* M;
end

% Compute the opponent channels using the effective weights
RG = wL_eff .* L - wM_eff .* M;
BY = wS_eff .* S - Ach;

% Explicitly assign NaN to the saturated/invalid pixels
Ach(invalidMask) = NaN;
RG(invalidMask)  = NaN;
BY(invalidMask)  = NaN;

% Reconstruct the spatial map
[H, W, ~] = size(isomerizationMap);
postRecepMap = zeros(H, W, 3);

% Channel 1: Luminance (Achromatic)
% Channel 2: Red-Green
% Channel 3: Blue-Yellow
postRecepMap(:,:,1) = Ach;
postRecepMap(:,:,2) = RG;
postRecepMap(:,:,3) = BY;

end