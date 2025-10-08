function visualFieldPoints = anglesFromIntrinsics(sensorPoints,fisheyeIntrinsics)
%
%
% Examples:
%{
    % Load the fisheyeIntrinsics first
    % Then we calculate the visual field angle position of every camera
    % sensor position
    nCols = 640; nRows = 480;
    [xg, yg] = meshgrid(1:nCols, 1:nRows);
    sensorPoints = [xg(:),yg(:)];
    visualFieldPoints = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);    
%}

% nRows = fisheyeIntrinsics.ImageSize(1);
% nCols = fisheyeIntrinsics.ImageSize(2);
% 
% [xGrid, yGrid] = meshgrid(1:nCols, 1:nRows);
    cx = fisheyeIntrinsics.DistortionCenter(1);
    cy = fisheyeIntrinsics.DistortionCenter(2);
    xC = sensorPoints(:,1) - cx; yC = sensorPoints(:,2) - cy;

    r = sqrt(xC.^2 + yC.^2);
    coeffs = fisheyeIntrinsics.MappingCoefficients;
    a0=coeffs(1); 
    a2=coeffs(2); 
    a3=coeffs(3); 
    a4=coeffs(4);

    theta = nan(size(r));
    for k=1:numel(r)
        func = @(t) a0*t + a2*t^3 + a3*t^4 + a4*t^5 - r(k);
        theta(k) = fzero(func,[0,pi]);
    end

    phi = atan2(yC,xC); R = 1;
    X = R*sin(theta).*cos(phi);
    Y = R*sin(theta).*sin(phi);
    Z = R*cos(theta);

    azi = rad2deg(atan(Y./Z));
    ele = rad2deg(atan(X./Z));

    visualFieldPoints = [azi,ele];

end