subjectIDs = {'2001','2003','2004','2006'};
dropBoxDir = getpref("lightLoggerAnalysis", "dropboxBaseDir");

for ss = 1:length(subjectIDs)
    thisFileName = ['FLIC_' subjectIDs{ss} '_walkIndoor_virtuallyFoveatedVideo.hdf5'];
    thisVideo = fullfile(dropBoxDir,'FLIC_analysis','lightLogger',thisFileName);

    [exponentMaps(:,:,ss), interceptMaps(:,:,ss), spdsByRegion(:,:,:,ss), frq, medianImage(:,:,ss)] = mapSPDs(thisVideo);
end

avgExponentMap = squeeze(mean(exponentMaps(:,:,1:4),3,'omitmissing'));
avgVarianceMap = squeeze(mean(interceptMaps(:,:,1:4),3,'omitmissing'));
avgSpdByRegion = squeeze(mean(spdsByRegion(:,:,:,1:4),4,'omitmissing'));

% Correct the image tilt
avgExponentMap = imrotate(avgExponentMap,-12.5,'bicubic','crop');
avgVarianceMap = imrotate(avgVarianceMap,-12.5,'bicubic','crop');

% Define the image resolution in degrees
fovDegrees = 120;
degPerPix = 120 / size(avgExponentMap,1);

% Nan out areas beyond an eliptical field of view
ellipseTransparentParams = [240, 240, 120000, .75, 0];
p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
[X, Y] = meshgrid(1:480, 1:480);
mask = double(myEllipse(X,Y)<1e-9);
avgExponentMap(mask==0)=nan;
avgVarianceMap(mask==0)=nan;

% Display the maps
figure
imagesc(avgExponentMap,[1.15 1.2]);
hold on
plot(240,240,'+k')
% Add circles at 5,10,20, and 40 degrees
for ii = [5, 10, 20, 40]
viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
end
idxStarts = 1:12:(480 - 24 + 1);
for ii = 1:39
    plot([idxStarts(ii),idxStarts(ii)],[1 480],'-','Color',[0.5 0.5 0.5]);
    plot([1 480],[idxStarts(ii),idxStarts(ii)],'-','Color',[0.5 0.5 0.5]);
end
title('Average Exponent')
axis square
colorbar

figure
imagesc(avgVarianceMap,[0.025 0.05]);
hold on
plot(240,240,'+k')
for ii = [5, 10, 20, 40]
viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1,'EnhanceVisibility',false);
end
for ii = 1:39
    plot([idxStarts(ii),idxStarts(ii)],[1 480],'-','Color',[0.5 0.5 0.5]);
    plot([1 480],[idxStarts(ii),idxStarts(ii)],'-','Color',[0.5 0.5 0.5]);
end
title('Average Contrast Variance')
axis square
colorbar

figure
loglog(frq,squeeze(avgSpdByRegion(20,20,:)),'-k');
hold on
loglog(frq,squeeze(avgSpdByRegion(31,20,:)),'-r');
legend({'center','periphery'});
ylabel('Power [contrast^2/Hz]');
xlabel('Frequency [log Hz]');
title('SPDs from the center and periphery');

