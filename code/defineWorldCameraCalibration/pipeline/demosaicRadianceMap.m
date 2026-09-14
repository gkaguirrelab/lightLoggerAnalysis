function rgbMap = demosaicRadianceMap(radianceMap, bayerPattern)
% DEMOSAICRADIANCEMAP Demosaics a Bayer-pattern radiance map into a 3D RGB image.
%   rgbMap = demosaicRadianceMap(radianceMap, bayerPattern) takes a 2D radiance map 
%   with a specified Bayer pattern and returns an H x W x 3 double matrix containing 
%   the interpolated Red, Green, and Blue channels using returnBayerIndices.
%
%   Pixels assigned an Inf value in the input retain their Inf value in their 
%   respective output channel after demosaicing.

    if nargin < 2 || isempty(bayerPattern)
        bayerPattern = "BGGR"; % Default Bayer pattern matching pipeline usage
    end

    [H, W] = size(radianceMap);
    
    % Create full coordinate grid for the output image
    [X, Y] = meshgrid(1:W, 1:H);
    
    % Obtain channel indices using returnBayerIndices (R=1, G=2, B=3)
    [rgbIdx{1}, rgbIdx{2}, rgbIdx{3}] = returnBayerIndices(radianceMap, bayerPattern);
    
    % Interpolate individual channels while preserving Inf values
    fullGridR = interpolateChannel(radianceMap, rgbIdx{1}, X, Y, H, W);
    fullGridG = interpolateChannel(radianceMap, rgbIdx{2}, X, Y, H, W);
    fullGridB = interpolateChannel(radianceMap, rgbIdx{3}, X, Y, H, W);
    
    % Assemble into a 3D H x W x 3 matrix (R, G, B)
    rgbMap = cat(3, fullGridR, fullGridG, fullGridB);

end

function fullGrid = interpolateChannel(radianceMap, channelIdx, X, Y, H, W)
    % Extract values for this channel from the raw radiance map
    subVal = radianceMap(channelIdx);
    
    % Identify which specific sub-pixel locations are Inf in the raw map
    infMaskSub = isinf(subVal);
    
    % Prepare work values by turning Inf and NaN into NaN for the interpolant
    workSub = subVal;
    workSub(isinf(workSub)) = NaN;
    
    % Filter out NaN values for fitting
    subX = X(channelIdx);
    subY = Y(channelIdx);
    validIdx = ~isnan(workSub);
    
    if ~any(validIdx(:))
        fullGrid = nan(H, W);
        return;
    end
    
    xData = subX(validIdx);
    yData = subY(validIdx);
    vData = workSub(validIdx);
    
    % Perform bilinear interpolation with nearest extrapolation for edge robustness
    F = scatteredInterpolant(xData, yData, vData, 'linear', 'nearest');
    fullGrid = F(X, Y);
    
    % Restore Inf values at their original sub-grid locations
    if any(infMaskSub(:))
        if islogical(channelIdx)
            fullInfMask = false(H, W);
            fullInfMask(channelIdx) = infMaskSub;
            fullGrid(fullInfMask) = Inf;
        else
            infLinearIndices = channelIdx(infMaskSub);
            fullGrid(infLinearIndices) = Inf;
        end
    end
end