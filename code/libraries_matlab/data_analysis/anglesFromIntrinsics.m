function [theta, phi, r] = anglesFromIntrinsics(nRows, nCols, fisheyeIntrinsics)
    [xg, yg] = meshgrid(1:nCols, 1:nRows);
    cx = fisheyeIntrinsics.DistortionCenter(1);
    cy = fisheyeIntrinsics.DistortionCenter(2);
    xC = xg - cx;  yC = yg - cy;

    r = sqrt(xC.^2 + yC.^2);
    c = fisheyeIntrinsics.MappingCoefficients;

    theta = nan(size(r));
    for k = 1:numel(r)
        theta(k) = fzero(@(t) c(1)*t + c(2)*t.^3 + c(3)*t.^4 + c(4)*t.^5 - r(k), [0, pi]);
    end
    % radians
    phi = atan2(yC, xC);
end