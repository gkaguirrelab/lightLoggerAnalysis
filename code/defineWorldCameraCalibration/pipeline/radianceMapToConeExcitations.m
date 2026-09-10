function lmsMap = radianceMapToConeExcitations(radianceMap)
% RADIANCEMAPTOCONEEXCITATIONS Converts an IMX219 absolute radiance map 
% (containing infs and NaNs) into a [3, height, width] matrix of L, M, and S 
% cone excitations, preserving input Inf/NaN locations in the output.
%
% Inputs:
%   radianceMap - 2D Bayer array containing absolute radiance values 
%                 (may include Inf or NaN values for saturated/invalid pixels)
%
% Outputs:
%   lmsMap      - 3D array of dimensions [3, height, width] where the 
%                 first dimension corresponds to [L; M; S] cone excitations,
%                 with Inf/NaN values preserved from the input.

    % 1. Identify and store the invalid pixel mask from the input
    invalidMask = isinf(radianceMap) | isnan(radianceMap);

    % 2. Load IMX219 spectral sensitivities[cite: 1]
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(paramFileName,'T');
    cameraT = table2array(T(:,["red" "green" "blue"]))';
    cameraWls = T.wls;
    cameraS = WlsToS(cameraWls);

    % 3. Generate the 2-degree cone fundamentals matching camera wavelengths
    fieldSizeDegrees = 30;
    observerAgeInYears = 30;
    pupilDiameterMm = 2;
    T_receptors = GetHumanPhotoreceptorSS(cameraS, ...
        {'LConeTabulatedAbsorbance2Deg', ...
        'MConeTabulatedAbsorbance2Deg', ...
        'SConeTabulatedAbsorbance2Deg'}, ...
        fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], [], [], []);

    % 4. Derive the 3x3 camera-to-cone linear transformation matrix
    M = T_receptors / cameraT;

    % 5. Sanitize a copy of the input map for safe interpolation
    cleanMap = radianceMap;
    cleanMap(invalidMask) = NaN;

    [h, w] = size(cleanMap);
    rgbImages = zeros(h, w, 3);

    % 6. Isolate Bayer indices and interpolate each channel independently[cite: 1]
    bayerPattern = "BGGR";
    [rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(cleanMap, bayerPattern);

    [X, Y] = meshgrid(1:w, 1:h);

    for cc = 1:3
        channelGrid = nan(h, w);
        channelGrid(rgbIdx{cc}) = cleanMap(rgbIdx{cc});
        
        validMask = ~isnan(channelGrid);
        
        if any(validMask(:))
            F = scatteredInterpolant(X(validMask), Y(validMask), channelGrid(validMask), 'linear', 'nearest');
            rgbImages(:,:,cc) = F(X, Y);
        else
            rgbImages(:,:,cc) = zeros(h, w);
        end
    end

    % 7. Reshape interpolated RGB maps to [3 x (pixels)] for matrix multiplication
    rgbFlat = reshape(rgbImages, [h * w, 3])'; 

    % 8. Map the continuous R, G, B radiance maps to L, M, and S cone excitations
    lmsFlat = M * rgbFlat; 
    lmsMap = reshape(lmsFlat, [3, h, w]);

    % 9. Restore the original Inf and NaN values across all channels of the output map
    for cc = 1:3
        channelSlice = squeeze(lmsMap(cc, :, :));
        channelSlice(invalidMask) = radianceMap(invalidMask);
        lmsMap(cc, :, :) = channelSlice;
    end

end