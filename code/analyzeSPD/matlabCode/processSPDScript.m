subjectIDs = {'2001','2003','2004'};%,'2005'};
dropBoxDir = getpref("lightLoggerAnalysis", "dropboxBaseDir");

for ss = 1:length(subjectIDs)
    thisFileName = ['FLIC_' subjectIDs{ss} '_walkIndoor_virtuallyFoveatedVideo.hdf5'];
    thisVideo = fullfile(dropBoxDir,'FLIC_analysis','lightLogger',thisFileName);

    [slopeMaps(:,:,ss), aucMaps(:,:,ss), spdsByRegion(:,:,:,ss), frq, medianImage(:,:,ss)] = mapSPDs(thisVideo);
end

avgSlopeMap = squeeze(mean(slopeMaps,3,'omitmissing'));
avgAUCMap = squeeze(mean(aucMaps,3,'omitmissing'));

% Correct the image tilt
avgSlopeMap = imrotate(avgSlopeMap,-15,'bicubic','crop');
avgAUCMap = imrotate(avgAUCMap,-15,'bicubic','crop');

% Average across the horizontal
%avgSlopeMap = (avgSlopeMap + fliplr(avgSlopeMap))/2;

% Nan out areas beyond an eliptical field of view
ellipseTransparentParams = [240, 250, 140000, .65, 0];
p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
[X, Y] = meshgrid(1:480, 1:480);
mask = double(myEllipse(X,Y)<1e-9);
avgSlopeMap(mask==0)=nan;
avgAUCMap(mask==0)=nan;

% Display the maps
figure
imagesc(avgSlopeMap,[-2.75 -2.25]);
hold on
plot(240,240,'+k')
viscircles([240 240], 50, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
viscircles([240 240], 100, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
viscircles([240 240], 150, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
title('Average Slope')
axis square
colorbar
figure
imagesc(avgAUCMap,[-4 -3.5]);
hold on
plot(240,240,'+k')
viscircles([240 240], 50, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
viscircles([240 240], 100, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
viscircles([240 240], 150, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
title('Average AUC')
axis square
colorbar
