function visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics)
% Convert image sensor points to visual-field angles using fisheye intrinsics
%
% Syntax:
%   visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics)
%
% Description:
%   This function converts image coordinates from a fisheye camera sensor
%   into angular visual-field coordinates. Given a set of sensor points
%   and a fisheye intrinsics calibration object, the function first shifts
%   the points relative to the distortion center, then computes each
%   point's radial distance on the sensor. Using the camera mapping
%   coefficients, it numerically inverts the fisheye radial mapping to
%   recover the incident angle theta for each point.
%
%   The recovered spherical coordinates are then converted into 3D unit
%   direction vectors and expressed as azimuth and elevation in degrees.
%   The returned coordinates are suitable for downstream use in retinal
%   remapping and virtual foveation routines.
%
% Inputs:
%   sensorPoints         - Numeric matrix. N-by-2 array of sensor/image
%                          coordinates, where each row is a point in pixel
%                          units [x, y].
%   fisheyeIntrinsics    - Camera intrinsics object / struct. Fisheye
%                          camera calibration object containing at least:
%                              .DistortionCenter
%                              .MappingCoefficients
%
% Outputs:
%   visualFieldPoints    - Numeric matrix. N-by-2 array of visual-field
%                          angles in degrees, where each row is:
%                              [azimuth, elevation]
%
% Examples:
%{
    intr = load('/path/to/intrinsics_calibration.mat');
    intr = intr.camera_intrinsics_calibration.results.Intrinsics;

    sensorPoints = [240 240; 300 220; 180 260];
    visualFieldPoints = anglesFromIntrinsics(sensorPoints, intr);
%}    
    
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

    rmax = a0*pi + a2*pi^3 + a3*pi^4 + a4*pi^5;
    rClamped = min(r, rmax - eps(rmax));  % keep strictly inside

    theta = nan(size(r));
    for k=1:numel(rClamped)
        func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - rClamped(k);
        theta(k) = fzero(func,[0,pi]);
    end

    %theta = interp1(rr, tt, r, 'pchip', NaN);
    
    phi = atan2(yC,xC); 
    R = 1;
    X = R*sin(theta).*cos(phi);
    Y = R*sin(theta).*sin(phi);
    Z = R*cos(theta);

    azi = rad2deg(atan2(X, Z));
    ele = rad2deg(asin(-Y));    % since X^2 + Y^2 + Z^2 = 1
    
    visualFieldPoints = [azi, ele];
end
