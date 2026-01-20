function plotGazeAccuracyNeon(startFixationId, endFixationId)
    % 1. Setup and Data Loading
    filePath = '/Users/samanthamontoya/Downloads/Timeseries Data + Scene Video-7/2026-01-16_20-48-06-8dd5fc8b/gaze.csv';
    data = readtable(filePath);
    
    ts_ns = data.timestamp_ns_; 
    fId   = data.fixationId;
    az    = data.azimuth_deg_;
    el    = data.elevation_deg_;
    
    % 2. Define Intended Gaze Positions
    baseTargets = [0,0; -15,15; -15,-15; 15,15; 15,-15; 0,15; 0,-15; -15,0; 15,0; ...
                   -7.5,7.5; -7.5,-7.5; 7.5,7.5; 7.5,-7.5; 0,7.5; 0,-7.5; -7.5,0; 7.5,0];
    intendedGaze = [baseTargets; baseTargets];
    numTargetsExpected = 34;
    
    % 3. Filter to the specific Fixation ID range
    startIndex = find(fId == startFixationId, 1, 'first');
    endIndex   = find(fId == endFixationId, 1, 'last');
    
    if isempty(startIndex) || isempty(endIndex)
        error('Specified Fixation IDs not found.');
    end
    
    ts_ns = ts_ns(startIndex:endIndex);
    fId   = fId(startIndex:endIndex);
    az    = az(startIndex:endIndex);
    el    = el(startIndex:endIndex);
    
    % 4. Identify discrete fixation blocks
    isValid = ~isnan(fId);
    diffValid = diff([0; isValid; 0]);
    blockStarts = find(diffValid == 1);
    blockEnds = find(diffValid == -1) - 1;
    
    actualGaze = [];
    if isempty(blockStarts), error('No fixations found in range.'); end
    
    % Initialize first target
    currTargetAz = az(blockStarts(1):blockEnds(1));
    currTargetEl = el(blockStarts(1):blockEnds(1));
    targetStartTime = ts_ns(blockStarts(1));
    
    % 5. Logic: Spatial Jump + Temporal Minimum
    for i = 2:length(blockStarts)
        thisBlockAz = median(az(blockStarts(i):blockEnds(i)), 'omitnan');
        thisBlockEl = median(el(blockStarts(i):blockEnds(i)), 'omitnan');
        
        prevMedAz = median(currTargetAz, 'omitnan');
        prevMedEl = median(currTargetEl, 'omitnan');
        
        dist = sqrt((thisBlockAz - prevMedAz)^2 + (thisBlockEl - prevMedEl)^2);
        timeElapsed = (ts_ns(blockStarts(i)) - targetStartTime) / 1e9; 
        
        % Split only if jump > 6 deg AND target has lasted at least 2.5s
        if dist > 6.0 && timeElapsed > 2.5
            actualGaze = [actualGaze; median(currTargetAz, 'omitnan'), median(currTargetEl, 'omitnan')];
            currTargetAz = az(blockStarts(i):blockEnds(i));
            currTargetEl = el(blockStarts(i):blockEnds(i));
            targetStartTime = ts_ns(blockStarts(i));
        else
            currTargetAz = [currTargetAz; az(blockStarts(i):blockEnds(i))];
            currTargetEl = [currTargetEl; el(blockStarts(i):blockEnds(i))];
        end
    end
    actualGaze = [actualGaze; median(currTargetAz, 'omitnan'), median(currTargetEl, 'omitnan')];
    
    % Ensure we match the intended targets (limit to 34)
    numFound = size(actualGaze, 1);
    if numFound > numTargetsExpected
        actualGaze = actualGaze(1:numTargetsExpected, :);
        numFound = numTargetsExpected;
    end
    intendedGaze = intendedGaze(1:numFound, :);

    % 6. Calculate Offsets for Title
    xOffsets = actualGaze(:,1) - intendedGaze(:,1);
    yOffsets = actualGaze(:,2) - intendedGaze(:,2);
    avgXOffset = mean(xOffsets);
    avgYOffset = mean(yOffsets);
    totalErr = mean(sqrt(xOffsets.^2 + yOffsets.^2));

    % 7. Plotting
    figure('Position', [100 100 800 750]);
    hold on;
    
    % Plot Intended (Red X)
    hInt = scatter(intendedGaze(:,1), intendedGaze(:,2), 150, 'r', 'x', 'LineWidth', 2);
    
    % Plot Actual (Cyan Circle)
    hAct = scatter(actualGaze(:,1), actualGaze(:,2), 100, 'cyan', 'filled', 'MarkerEdgeColor', 'k');
    
    % Error lines and numbering
    for i = 1:numFound
        line([intendedGaze(i,1), actualGaze(i,1)], [intendedGaze(i,2), actualGaze(i,2)], ...
             'Color', [0.5 0.5 0.5], 'LineStyle', '-');
         
        text(actualGaze(i,1) + 0.6, actualGaze(i,2) + 0.6, num2str(i), ...
            'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Construct Title with Offsets
    titleStr = {sprintf('Gaze Accuracy (N=%d)', numFound), ...
                sprintf('Avg Offset: X=%.2f°, Y=%.2f° | Mean Error: %.2f°', ...
                avgXOffset, avgYOffset, totalErr)};
    title(titleStr);
    
    xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)');
    legend([hInt, hAct], {'Intended Target', 'Actual Gaze Median'}, 'Location', 'northeastoutside');
    
    grid on; axis equal;
    xlim([-22 22]); ylim([-22 22]);
    hold off;
    
    % Final Stats to Command Window
    fprintf('Found %d targets.\n', numFound);
    fprintf('Mean X Offset: %.2f°, Mean Y Offset: %.2f°\n', avgXOffset, avgYOffset);
    fprintf('Total Mean Error: %.2f°\n', totalErr);

    % 7. Plot again with offset subtracted
    figure('Position', [100 100 800 750]);
    hold on;
    
    % Plot Intended (Red X)
    hInt = scatter(intendedGaze(:,1), intendedGaze(:,2), 150, 'r', 'x', 'LineWidth', 2);
    
    % Plot Actual - offset (Cyan Circle)
    hAct = scatter(actualGaze(:,1) - avgXOffset, actualGaze(:,2) - avgYOffset, 100, 'cyan', 'filled', 'MarkerEdgeColor', 'k');
    
    % Error lines and numbering
    for i = 1:numFound
        line([intendedGaze(i,1), actualGaze(i,1)- avgXOffset], [intendedGaze(i,2), actualGaze(i,2)- avgYOffset], ...
             'Color', [0.5 0.5 0.5], 'LineStyle', '-');
         
        text(actualGaze(i,1) + 0.6, actualGaze(i,2) + 0.6, num2str(i), ...
            'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Construct Title with Offsets
    titleStr = {sprintf('Gaze Accuracy (N=%d)', numFound), ...
                sprintf('Avg Offset: X=%.2f°, Y=%.2f° | Mean Error: %.2f°', ...
                avgXOffset, avgYOffset, totalErr)};
    title(titleStr);
    
    xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)');
    legend([hInt, hAct], {'Intended Target', 'Actual Gaze Median minus avg offset'}, 'Location', 'northeastoutside');
    
    grid on; axis equal;
    xlim([-22 22]); ylim([-22 22]);
    hold off;
end