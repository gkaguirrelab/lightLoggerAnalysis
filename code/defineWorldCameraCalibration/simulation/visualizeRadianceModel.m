function visualizeRadianceModel(radianceModel,radianceModelS)
    
    % 2. Generate a custom full-elevation grid to replace debugVisuals
    azGrid = linspace(-pi, pi, 1000);
    elGrid = linspace(-pi/2, pi/2, 500);
    [plotAZ, plotEL] = meshgrid(azGrid, elGrid);
    
    % Evaluate the radiance model across the full grid
    spectralSurf = radianceModel(plotAZ, plotEL);
    plotBroadband = squeeze(sum(spectralSurf, 1)) * radianceModelS(2);
    
    % 3. Load camera intrinsics to map the sensor field of view
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
    
    % Extract the coordinates for the sensor's outer perimeter
    azPerimeter = [azimuthMap(1, :), azimuthMap(2:end, end)', ...
                   azimuthMap(end, end-1:-1:1), azimuthMap(end-1:-1:1, 1)'];
    elPerimeter = [elevationMap(1, :), elevationMap(2:end, end)', ...
                   elevationMap(end, end-1:-1:1), elevationMap(end-1:-1:1, 1)'];
    
    % 4. Plot 1: Radiance across the full field with the camera FOV outline
    figure('Name', 'Full-Field Radiance Map', 'Color', 'w');
    
    % Plot the full-elevation surface
    surf(plotAZ, plotEL, plotBroadband, 'EdgeColor', 'none');
    view(2);
    axis tight;
    colorbar;
    hold on;
    
    % Overlay the dotted outline slightly above the maximum radiance surface
    zMax = max(plotBroadband(:));
    plot3(azPerimeter, elPerimeter, repmat(zMax * 1.1, size(azPerimeter)), ...
          'w--', 'LineWidth', 2, 'DisplayName', 'Camera FOV Bound');
      
    title('Broadband Radiance (W/m^2/sr)');
    xlabel('Azimuth (radians)');
    ylabel('Elevation (radians)');
    legend('Location', 'best');
    hold off;
    
end