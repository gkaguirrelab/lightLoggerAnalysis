function plotSPDs(activityData, options)
    arguments 
        activityData; 
        options.fovDegrees = 120; 
    end 

    % Extract the only 
    field_names = fieldnames(activityData);
    activityName = field_names{1};

    exponentMap = activityData.(activityName).exponentMap;
    varianceMap = activityData.(activityName).varianceMap;
    spdByRegion = activityData.(activityName).spdByRegion;
    frq = activityData.(activityName).frq;
    medianImage = activityData.(activityName).medianImage;
    frameDropVector = activityData.(activityName).frameDropVector;

    % Nan out areas beyond an elliptical field of view
    ellipseTransparentParams = [240, 240, 120000, .75, 0];
    p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
    myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
    [X, Y] = meshgrid(1:480, 1:480);
    mask = double(myEllipse(X,Y) < 1e-9);
    exponentMap(mask == 0) = nan;
    exponentMap(mask == 0) = nan;

    degPerPix = fovDegrees / size(exponentMap,1);

    % Display the maps
    figure
    imagesc(-exponentMap);
    hold on
    plot(240,240,'+k')
    for ii = [5, 10, 20, 40]
        viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
    end
    idxStarts = 1:12:(480 - 24 + 1);
    for ii = 1:39
        plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
        plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
    end
    title(sprintf('Exponent - %s', activityName))
    axis square
    colorbar

    figure
    imagesc(varianceMap, [0.015 0.035]);
    hold on
    plot(240,240,'+k')
    for ii = [5, 10, 20, 40]
        viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
    end
    for ii = 1:39
        plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
        plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
    end
    title(sprintf('Contrast Variance - %s', activityName))
    axis square
    colorbar

    figure
    loglog(frq, squeeze(spdByRegion(20,20,:)), '-k');
    hold on
    loglog(frq, squeeze(spdByRegion(31,20,:)), '-r');
    plot([10^0 10^1.5], [10^-2 10^-5], ':k')
    legend({'center','periphery'});
    ylabel('Power [contrast^2/Hz]');
    xlabel('Frequency [log Hz]');
    title(sprintf('SPDs from the center and periphery - %s', activityName));

end