function isomerizationMap = computeConeIsomerizationMap(integratedRadianceMap, coneMapVar, options)
% COMPUTECONEISOMERIZATIONMAP Converts an RGB integrated radiance map into
% an LMS cone isomerization rate map, allowing dynamic fixation and pupil
% adjustments.
%
% Inputs:
%   integratedRadianceMap - H x W x 3 RGB matrix in W/m^2/sr
%   coneMapVar            - Struct containing the 1D LUT transform table
%
% Name-Value Arguments:
%   fixationAzimuth   - Azimuth offset of fixation in degrees (default: 0)
%   fixationElevation - Elevation offset of fixation in degrees (default: 0)
%   pupilDiameterMm   - Pupil diameter for this specific image in mm (default: 3.0)

arguments
    integratedRadianceMap (:,:,3) double
    coneMapVar (1,1) struct
    options.fixationAzimuth (1,1) double = 0
    options.fixationElevation (1,1) double = 0
    options.pupilDiameterMm (1,1) double = 3.0
end

persistent unitDirections

% Load the unit vectors defining the camera geometry
if isempty(unitDirections)
    mapFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'cameraToVisualAngles.mat');
    load(mapFileName, 'unitDirections');
end

[H, W, ~] = size(integratedRadianceMap);

% Establish the pupil area scalar relative to the canonical mapping base
pupilScalar = (options.pupilDiameterMm / coneMapVar.basePupilDiameterMm)^2;

% Calculate the dynamic eccentricity map based on fixation offset
az = deg2rad(options.fixationAzimuth);
el = deg2rad(options.fixationElevation);

% Create the fixation unit vector using the established spatial convention
fixationVector = [cos(el)*sin(az), -sin(el), cos(el)*cos(az)];
fixationVector = reshape(fixationVector, 1, 1, 3);

% The dot product of the unit directions and the fixation vector yields the
% cosine of the angular distance (eccentricity) from the fovea
dotProducts = sum(unitDirections .* fixationVector, 3);
dotProducts = min(max(dotProducts, -1), 1); % Clamp to prevent acos precision errors
dynamicEccMap = rad2deg(acos(dotProducts));

% Interpolate the 3x3 matrices from the 1D LUT. This maps the HxW dynamic
% eccentricities into an (H*W) x 9 matrix
T_flat = interp1(coneMapVar.eccGrid, coneMapVar.transformTable, dynamicEccMap(:), 'linear', 'extrap');

% Reshape back into the H x W x 3 x 3 spatial matrix format
T_map = reshape(T_flat, H, W, 3, 3);

% Perform spatially varying matrix multiplication
isomerizationMap = zeros(H, W, 3);

for coneClass = 1:3 % 1=L, 2=M, 3=S
    % Extract the spatially varying 1x3 vector for this cone class[cite: 5]
    transformWeights = squeeze(T_map(:, :, coneClass, :));

    % Multiply weights against radiance map, sum across channels, and apply pupil scalar
    isoChannel = sum(integratedRadianceMap .* transformWeights, 3) .* pupilScalar;

    isomerizationMap(:,:,coneClass) = isoChannel;
end

% Bayesian Imputation of Negative Isomerization Rates. First flatten the
% map to an N x 3 array for cross-channel statistical modeling
pixels = reshape(isomerizationMap, [], 3);

% Identify strictly valid pixels (all three cone classes > 0)
validMask = all(pixels > 0, 2);
validPixels = pixels(validMask, :);

if ~isempty(validPixels)
    % Extract global prior statistics in log space
    logValid = log(validPixels);
    mu = mean(logValid, 1)';
    S = cov(logValid);

    % Establish the log-floor limit for the truncated normal expectation
    f_log = min(logValid, [], 1)';

    logFixedPixels = log(pixels);

    % Iterate through L, M, and S channels to impute negative values
    for cTarget = 1:3
        % Find all pixels where the target cone class requires imputation
        targetMask = pixels(:, cTarget) <= 0;
        activeImputeIndices = find(targetMask);

        for pIdx = activeImputeIndices'
            % Identify which other cone classes at this specific pixel are valid
            k_cols = find(pixels(pIdx, :) > 0);

            if ~isempty(k_cols)
                % Extract sub-matrices for the conditional expectation
                muK = mu(k_cols);
                SK = S(k_cols, k_cols);
                SSK = S(cTarget, k_cols);

                % Add a tiny regularization term to ensure matrix invertibility
                SKInv = inv(SK + 1e-6 * eye(length(k_cols)));
                k_vals = logFixedPixels(pIdx, k_cols)';

                % Calculate conditional mean and variance based on the valid channels
                muXs = mu(cTarget) + SSK * SKInv * (k_vals - muK);
                Sxs = S(cTarget, cTarget) - SSK * SKInv * SSK';
            else
                % If no channels are valid, fall back to global mean and variance
                muXs = mu(cTarget);
                Sxs = S(cTarget, cTarget);
            end

            Sxs = max(Sxs, 1e-8);
            stdXs = sqrt(Sxs);

            % Calculate the expected value of a lower-truncated normal distribution
            % bounded above by the minimum observable valid log rate
            zScore = (f_log(cTarget) - muXs) / stdXs;
            Z = normcdf(zScore);

            if Z < 1e-15
                expectedVal_log = f_log(cTarget);
            else
                numeratorTerm = (stdXs / sqrt(2 * pi)) * exp(-(zScore^2) / 2);
                expectedVal_log = muXs - numeratorTerm / Z;
            end

            % Exponentiate and replace the negative artifact
            pixels(pIdx, cTarget) = exp(expectedVal_log);
        end
    end
    % Rebuild the spatial map
    isomerizationMap = reshape(pixels, H, W, 3);
else
    % Fallback thresholding if the image contains absolutely no valid pixels
    isomerizationMap(isomerizationMap < 0) = 0;
end

end
