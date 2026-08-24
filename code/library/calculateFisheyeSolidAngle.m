function solidAngleSteradians = calculateFisheyeSolidAngle(fisheyeIntrinsics, options)
% Calculate the solid angle covered by a calibrated fisheye image.
%
% Syntax:
%   solidAngleSteradians = calculateFisheyeSolidAngle(fisheyeIntrinsics)
%   solidAngleSteradians = calculateFisheyeSolidAngle(..., NumAzimuthSamples=N)
%
% Description:
%   The rectangular sensor boundary is expressed as a radial limit for each
%   azimuth around the distortion center. The fisheye mapping polynomial is
%   inverted at that boundary to obtain theta_max(phi), and the enclosed
%   spherical area is integrated as
%
%       integral_0^(2*pi) (1 - cos(theta_max(phi))) dphi.
%
%   This calculation uses only the camera intrinsics and image boundary. It
%   does not construct or differentiate a per-pixel visual-angle map.
%
% Inputs:
%   fisheyeIntrinsics   - camera.FisheyeIntrinsics object containing
%                         MappingCoefficients, DistortionCenter, and
%                         ImageSize.
%
% Name-Value Inputs:
%   NumAzimuthSamples   - Odd integer number of uniformly spaced boundary
%                         azimuths used for numerical integration. Default
%                         is 100001.
%
% Outputs:
%   solidAngleSteradians - Scalar field of view in steradians.

    arguments
        fisheyeIntrinsics
        options.NumAzimuthSamples (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.NumAzimuthSamples, 1001)} = 100001
    end

    mappingCoefficients = double(fisheyeIntrinsics.MappingCoefficients(:).');
    distortionCenter = double(fisheyeIntrinsics.DistortionCenter(:).');
    imageSize = double(fisheyeIntrinsics.ImageSize(:).');

    if numel(mappingCoefficients) ~= 4 || any(~isfinite(mappingCoefficients))
        error('calculateFisheyeSolidAngle:InvalidMappingCoefficients', ...
            'MappingCoefficients must contain four finite values.');
    end
    if numel(distortionCenter) ~= 2 || any(~isfinite(distortionCenter))
        error('calculateFisheyeSolidAngle:InvalidDistortionCenter', ...
            'DistortionCenter must contain two finite values.');
    end
    if numel(imageSize) ~= 2 || any(~isfinite(imageSize)) || any(imageSize <= 0)
        error('calculateFisheyeSolidAngle:InvalidImageSize', ...
            'ImageSize must contain two positive finite values.');
    end

    % Pixel centers occupy x=1:columns and y=1:rows, so their physical
    % footprint extends half a pixel beyond the outermost centers.
    xLimits = [0.5, imageSize(2) + 0.5];
    yLimits = [0.5, imageSize(1) + 0.5];
    cx = distortionCenter(1);
    cy = distortionCenter(2);
    if cx <= xLimits(1) || cx >= xLimits(2) || cy <= yLimits(1) || cy >= yLimits(2)
        error('calculateFisheyeSolidAngle:DistortionCenterOutsideImage', ...
            'DistortionCenter must lie inside the image footprint.');
    end

    phi = linspace(0, 2*pi, options.NumAzimuthSamples);
    cosPhi = cos(phi);
    sinPhi = sin(phi);

    % Find where each ray from the distortion center first intersects the
    % rectangular image boundary.
    radialX = inf(size(phi));
    positiveX = cosPhi > 0;
    negativeX = cosPhi < 0;
    radialX(positiveX) = (xLimits(2) - cx) ./ cosPhi(positiveX);
    radialX(negativeX) = (cx - xLimits(1)) ./ -cosPhi(negativeX);

    radialY = inf(size(phi));
    positiveY = sinPhi > 0;
    negativeY = sinPhi < 0;
    radialY(positiveY) = (yLimits(2) - cy) ./ sinPhi(positiveY);
    radialY(negativeY) = (cy - yLimits(1)) ./ -sinPhi(negativeY);
    boundaryRadius = min(radialX, radialY);

    % Invert r(theta) = a0*theta + a2*theta^3 + a3*theta^4 + a4*theta^5
    % with vectorized Newton iterations. The linear term is the natural
    % initial approximation for calibrated fisheye models.
    a0 = mappingCoefficients(1);
    a2 = mappingCoefficients(2);
    a3 = mappingCoefficients(3);
    a4 = mappingCoefficients(4);
    if a0 <= 0
        error('calculateFisheyeSolidAngle:InvalidLinearCoefficient', ...
            'The leading mapping coefficient must be positive.');
    end

    theta = boundaryRadius ./ a0;
    for iteration = 1:12
        mappedRadius = a0.*theta + a2.*theta.^3 + a3.*theta.^4 + a4.*theta.^5;
        mappingDerivative = a0 + 3*a2.*theta.^2 + 4*a3.*theta.^3 + 5*a4.*theta.^4;
        if any(mappingDerivative <= 0 | ~isfinite(mappingDerivative))
            error('calculateFisheyeSolidAngle:NonMonotonicMapping', ...
                'The fisheye mapping is not monotonic over the image footprint.');
        end
        theta = theta - (mappedRadius - boundaryRadius) ./ mappingDerivative;
    end

    residual = a0.*theta + a2.*theta.^3 + a3.*theta.^4 + a4.*theta.^5 - boundaryRadius;
    if any(~isfinite(theta)) || any(theta < 0 | theta > pi) || max(abs(residual)) > 1e-8
        error('calculateFisheyeSolidAngle:MappingInversionFailed', ...
            'Could not reliably invert the fisheye mapping at the image boundary.');
    end

    solidAngleSteradians = trapz(phi, 1 - cos(theta));
end
