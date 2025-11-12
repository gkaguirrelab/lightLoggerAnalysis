function visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics)
    arguments
        sensorPoints
        fisheyeIntrinsics
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

    tt = linspace(0, pi, 2000);
    rr = a0*tt + a2*tt.^3 + a3*tt.^4 + a4*tt.^5;
    assert(all(diff(rr)>0), 'Non-monotonic r(θ) on [0,π].');

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

    azi = rad2deg(atan2(X, Z));
    ele = rad2deg(asin(-Y));    % since X^2 + Y^2 + Z^2 = 1
    
    visualFieldPoints = [azi, ele];
end
