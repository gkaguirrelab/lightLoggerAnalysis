function visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics, offset)
    arguments
        sensorPoints
        fisheyeIntrinsics
        offset (1,2) double = [0 0]   % [dx, dy] shift
    end

    % Apply offset to distortion center
    cx = fisheyeIntrinsics.DistortionCenter(1);
    cy = fisheyeIntrinsics.DistortionCenter(2);

    % Shift coords
    xC = sensorPoints(:,1) - cx;
    yC = sensorPoints(:,2) - cy;

    r = sqrt(xC.^2 + yC.^2);
    coeffs = fisheyeIntrinsics.MappingCoefficients;
    a0=coeffs(1); a2=coeffs(2); a3=coeffs(3); a4=coeffs(4);

    theta = nan(size(r));
    for k=1:numel(r)
        func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - r(k);
        theta(k) = fzero(func,[0,pi]);
    end
    
    phi = atan2(yC,xC); 
    R = 1;
    X = R*sin(theta).*cos(phi);
    Y = R*sin(theta).*sin(phi);
    Z = R*cos(theta);

    azi = rad2deg(atan(Y./Z)) + offset(1);
    ele = rad2deg(atan(X./Z)) + offset(2);
    
    visualFieldPoints = [azi, ele];
end
