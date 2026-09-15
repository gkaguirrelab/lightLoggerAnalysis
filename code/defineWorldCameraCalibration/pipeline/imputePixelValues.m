function rawFixed = imputePixelValues(radianceMap, bayerPattern)
% Implements a Bayesian estimation of linearized sensor values for pixels
% at ceiling (Inf) or floor (zero). Based upon Zhang & Brainard approach:
%
% Zhang X, Brainard DH. Estimation of saturated pixel values in digital
% color imaging. Journal of the Optical Society of America A. 2004 Dec
% 1;21(12):2301-10.
%
% Modified to add imputation of floor values, and to consider the
% distribution of pixel values in the log transformed space.
%

if nargin < 2 || isempty(bayerPattern)
    bayerPattern = "BGGR";
end

% Interpolate to get cross-channel conditioning data
[H, W] = size(radianceMap);
[X, Y] = meshgrid(1:W, 1:H);
[bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(radianceMap, bayerPattern);

rgbMap = zeros(H, W, 3);
for c = 1:3
    subVal = radianceMap(bayerIdx{c});
    infMask = isinf(subVal);
    floorMask = (subVal == 0);

    workSub = subVal;
    workSub(infMask | floorMask) = NaN;

    validIdx = ~isnan(workSub);
    subX = X(bayerIdx{c}); subY = Y(bayerIdx{c});

    if any(validIdx(:))
        F = scatteredInterpolant(subX(validIdx), subY(validIdx), workSub(validIdx), 'linear', 'nearest');
        rgbMap(:,:,c) = F(X, Y);
    end

    % Restore Inf and 0 at sub-grid locations so the imputation step can find them
    channelGrid = rgbMap(:,:,c);
    channelGrid(bayerIdx{c}(infMask)) = Inf;
    channelGrid(bayerIdx{c}(floorMask)) = 0;
    rgbMap(:,:,c) = channelGrid;
end

% Extract Prior Statistics
pixels = reshape(rgbMap, [], 3);
validPixels = pixels(~any(isinf(pixels) | pixels == 0, 2), :);
logValid = log(validPixels);
mu = mean(logValid, 1)';
S = cov(logValid);

s_log = zeros(3, 1); f_log = zeros(3, 1);
for c = 1:3
    validC = pixels(pixels(:, c) > 0 & ~isinf(pixels(:, c)), c);
    if isempty(validC)
        s_log(c) = log(1.0);
        f_log(c) = log(1e-4);
    else
        s_log(c) = log(max(validC));
        f_log(c) = log(min(validC));
    end
end

% Targeted Imputation and Direct Remosaicing
rawFixed = radianceMap;
logFixedPixels = log(pixels);

for cTarget = 1:3
    targetMask = (isinf(radianceMap) | radianceMap == 0);
    bayerTargetMask = false(H, W);
    bayerTargetMask(bayerIdx{cTarget}) = true;

    activeImputeIndices = find(targetMask & bayerTargetMask);

    for pIdx = activeImputeIndices'
        valid_k_mask = isfinite(logFixedPixels(pIdx, :));
        valid_k_mask(cTarget) = false;
        valid_k_cols = find(valid_k_mask);

        if ~isempty(valid_k_cols)
            muK = mu(valid_k_cols);
            SK = S(valid_k_cols, valid_k_cols);
            SSK = S(cTarget, valid_k_cols);
            SKInv = inv(SK + 1e-6 * eye(length(valid_k_cols)));
            k = logFixedPixels(pIdx, valid_k_cols)';

            muXs = mu(cTarget) + SSK * SKInv * (k - muK);
            Sxs = S(cTarget, cTarget) - SSK * SKInv * SSK';
        else
            muXs = mu(cTarget);
            Sxs = S(cTarget, cTarget);
        end

        Sxs = max(Sxs, 1e-8);
        stdXs = sqrt(Sxs);

        if isinf(radianceMap(pIdx))
            zScore = (s_log(cTarget) - muXs) / stdXs;
            Z = 1 - normcdf(zScore);

            if Z < 1e-15
                expectedVal_log = s_log(cTarget);
            else
                numeratorTerm = (stdXs / sqrt(2 * pi)) * exp(-(zScore^2) / 2);
                expectedVal_log = muXs + numeratorTerm / Z;
            end
        else
            zScore = (f_log(cTarget) - muXs) / stdXs;
            Z = normcdf(zScore);

            if Z < 1e-15
                expectedVal_log = f_log(cTarget);
            else
                numeratorTerm = (stdXs / sqrt(2 * pi)) * exp(-(zScore^2) / 2);
                expectedVal_log = muXs - numeratorTerm / Z;
            end
        end

        rawFixed(pIdx) = exp(expectedVal_log);
    end
end
end