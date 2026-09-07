function channelEnergy = calculateChannelRadiantEnergy(radianceMap)
% CALCULATECHANNELRADIANTENERGY Computes the total radiant energy in R, G, and B 
% channels from a camera radiance map using calibrated pixel solid angles.
%
% Returned values are in Watts / m2.
%   channelEnergy = calculateChannelRadiantEnergy(radianceMap)

    % 1. Load camera intrinsics and pixel solid angle mapping (Stage 6)
    projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');
    intrinsicsPath = fullfile(projectRoot, 'derived', 'arducamB0392cameraIntrinsics.mat');
    
    load(intrinsicsPath, 'arducamB0392cameraIntrinsics');
    fisheyeIntrinsics = arducamB0392cameraIntrinsics.results.Intrinsics;
    
    % 2. Determine grid dimensions and compute/retrieve pixel solid angles (steradians)
    [rows, columns] = size(radianceMap);
    [xCoordinates, yCoordinates] = meshgrid(1:columns, 1:rows);
    sensorPoints = [xCoordinates(:), yCoordinates(:)];
    
    % Map pixel coordinates to visual angles
    visualAngles = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    elevationDeg = reshape(visualAngles(:, 2), rows, columns);
    
    % Compute differential solid angle (steradians) per pixel based on angular spacing 
    % and the cosine of elevation for spherical projection geometry.
    % (Alternatively, if stored directly in the intrinsics structure, load it here).
    dTheta = deg2rad(mean(diff(unique(elevationDeg(:, 1))))); 
    dPhi = deg2rad(mean(diff(unique(elevationDeg(1, :)))));
    solidAngles = abs(dTheta * dPhi .* cos(deg2rad(elevationDeg)));

    % 3. Define Bayer filter masks for IMX219 (RGGB pattern)
    redMask   = false(rows, columns); 
    greenMask = false(rows, columns); 
    blueMask  = false(rows, columns); 
    
    redMask(1:2:end, 1:2:end)       = true; % R
    greenMask(1:2:end, 2:2:end)     = true; % G (row 1)
    greenMask(2:2:end, 1:2:end)     = true; % G (row 2)
    blueMask(2:2:end, 2:2:end)      = true; % B

    % 4. Integrate radiance over pixel solid angles for each channel
    % Radiance units: W/m^2/sr. Multiplying by solid angle (sr) yields W/m^2 (flux density).
    redEnergy   = sum(radianceMap(redMask)   .* solidAngles(redMask));
    greenEnergy = sum(radianceMap(greenMask) .* solidAngles(greenMask));
    blueEnergy  = sum(radianceMap(blueMask)  .* solidAngles(blueMask));

    % Return results in a structured format
    channelEnergy = struct(...
        'Red', redEnergy, ...
        'Green', greenEnergy, ...
        'Blue', blueEnergy ...
    );
end