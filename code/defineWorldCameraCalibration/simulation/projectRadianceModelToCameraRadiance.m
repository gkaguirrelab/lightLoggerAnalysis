function radianceMap = projectRadianceModelToCameraRadiance(radianceModel,radianceModelS)
        
    % Load camera intrinsics to map the sensor field of view
    projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');
    intrinsicsPath = fullfile(projectRoot, 'derived', 'arducamB0392cameraIntrinsics.mat');
    
    intrinsicsData = load(intrinsicsPath, 'arducamB0392cameraIntrinsics');
    fisheyeIntrinsics = intrinsicsData.arducamB0392cameraIntrinsics.results.Intrinsics;
    
    % Define the 640x480 sensor grid
    rows = 480;
    columns = 640;
    [xCoordinates, yCoordinates] = meshgrid(1:columns, 1:rows);
    sensorPoints = [xCoordinates(:), yCoordinates(:)];
    
    % Convert pixel centers to visual angles using the calibrated intrinsics
    visualAngles = anglesFromIntrinsics(sensorPoints, fisheyeIntrinsics);
    visualAngles = reshape(visualAngles, rows, columns, 2);
    
    azimuthMap = deg2rad(visualAngles(:, :, 1));
    elevationMap = deg2rad(visualAngles(:, :, 2));
        
    spectralRadianceMap = radianceModel(azimuthMap, elevationMap);
    radianceMap = squeeze(sum(spectralRadianceMap, 1)) * radianceModelS(2);

end