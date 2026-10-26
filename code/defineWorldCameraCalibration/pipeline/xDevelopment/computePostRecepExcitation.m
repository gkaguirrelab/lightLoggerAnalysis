function [postRecepMap, coneWeights] = computePostRecepExcitation(isomerizationMap, options)
% COMPUTEPOSTRECEXCITATION Transforms an LMS cone isomerization rate map
% into a Linear Opponent Excitation map representing post-receptoral
% pathways.
%
% We assume a set of cone weights that are close to those empirically
% observed to result from assuming Von Kries adaptation in measurements
% made under outdoor lighting conditions.
%
% Inputs:
%   isomerizationMap - H x W x 3 matrix of LMS isomerization rates
%
% Name-Value Arguments:
%   balanceNeutral - Boolean to dynamically adapt weights to the scene mean (default: false)
%   wL_eff         - Fixed L-cone weight when not balancing neutral (default: 1)
%   wM_eff         - Fixed M-cone weight when not balancing neutral (default: 1.33)
%   wS_eff         - Fixed S-cone weight when not balancing neutral (default: 15)

arguments
    isomerizationMap (:,:,3) double
    options.balanceNeutral (1,1) logical = false
    options.wL_eff (1,1) double = 1
    options.wM_eff (1,1) double = 3
    options.wS_eff (1,1) double = 90
end

% Extract the individual cone class maps
L = isomerizationMap(:,:,1);
M = isomerizationMap(:,:,2);
S = isomerizationMap(:,:,3);

% Identify pixels that contain Inf (unimputed saturation) or NaN
invalidMask = isinf(L) | isinf(M) | isinf(S) | isnan(L) | isnan(M) | isnan(S);

if options.balanceNeutral %
    % 1. Use the image mean as a "Gray World" neutral reference
    meanL = mean(L(~invalidMask));
    meanM = mean(M(~invalidMask));
    meanS = mean(S(~invalidMask));
    
    % 2. Anchor the L-weight to 1.0 to preserve absolute physical scale
    wL_eff = 1.0; %
    
    % 3. Calculate wM so that L and M are balanced (RG = 0 at neutral)
    wM_eff = (wL_eff * meanL) / meanM;
    
    % 4. Compute the adapted Achromatic (Luminance) channel
    Ach = wL_eff .* L + wM_eff .* M;
    
    % 5. Calculate wS so that S balances Ach (BY = 0 at neutral)
    meanAch = wL_eff * meanL + wM_eff * meanM;
    wS_eff = meanAch / meanS; %
else
    % Fall back to user-provided or default fixed population weights
    wL_eff = options.wL_eff;
    wM_eff = options.wM_eff;
    wS_eff = options.wS_eff;
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

% Package the weights into an output variable
coneWeights = [wL_eff, wM_eff, wS_eff];

end