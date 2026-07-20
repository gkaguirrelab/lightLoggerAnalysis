%% CompareFieldingMaps.m

clear; close all; clc

oldData = load("fielding_correction_map_RGB_OLD.mat");
newData = load("fielding_correction_map_RGB.mat");

oldMap = double(oldData.correctionMap);
newMap = double(newData.correctionMap);

if ~isequal(size(oldMap), size(newMap))
    error("Correction maps have different dimensions.")
end

differenceMap = newMap - oldMap;
absoluteDifference = abs(differenceMap);

fprintf('\n----- FIELDING MAP COMPARISON -----\n')

fprintf('\nOld map:\n')
fprintf('Size: %d x %d\n', size(oldMap,1), size(oldMap,2))
fprintf('Range: %.6f to %.6f\n', min(oldMap(:)), max(oldMap(:)))
fprintf('Mean: %.6f\n', mean(oldMap(:)))
fprintf('SD: %.6f\n', std(oldMap(:)))

fprintf('\nNew map:\n')
fprintf('Size: %d x %d\n', size(newMap,1), size(newMap,2))
fprintf('Range: %.6f to %.6f\n', min(newMap(:)), max(newMap(:)))
fprintf('Mean: %.6f\n', mean(newMap(:)))
fprintf('SD: %.6f\n', std(newMap(:)))

fprintf('\nPixel-wise differences:\n')
fprintf('Maximum absolute difference: %.6f\n', ...
    max(absoluteDifference(:)))
fprintf('Mean absolute difference: %.6f\n', ...
    mean(absoluteDifference(:)))
fprintf('RMSE: %.6f\n', ...
    sqrt(mean(differenceMap(:).^2)))
fprintf('Mean signed difference: %.6f\n', ...
    mean(differenceMap(:)))

R = corrcoef(oldMap(:), newMap(:));
fprintf('Correlation: %.6f\n', R(1,2))

%% Show old and new maps using the same color scale

sharedCLim = [
    min([oldMap(:); newMap(:)]), ...
    max([oldMap(:); newMap(:)])
];

figure

subplot(1,2,1)
imagesc(oldMap)
axis image
clim(sharedCLim)
colorbar
title('Original fielding map')
xlabel('X Pixel')
ylabel('Y Pixel')

subplot(1,2,2)
imagesc(newMap)
axis image
clim(sharedCLim)
colorbar
title('Gamma-corrected fielding map')
xlabel('X Pixel')
ylabel('Y Pixel')

%% Difference map

figure
imagesc(differenceMap)
axis image
colorbar
xlabel('X Pixel')
ylabel('Y Pixel')
title('New minus old fielding map')

%% Absolute difference map

figure
imagesc(absoluteDifference)
axis image
colorbar
xlabel('X Pixel')
ylabel('Y Pixel')
title('Absolute difference between fielding maps')

%% Scatter comparison

figure
scatter(oldMap(:), newMap(:), 4, '.')
hold on

limits = [
    min([oldMap(:); newMap(:)]), ...
    max([oldMap(:); newMap(:)])
];

plot(limits, limits, 'k--', 'LineWidth', 1.5)

xlim(limits)
ylim(limits)
axis square
grid on

xlabel('Original correction factor')
ylabel('Gamma-corrected correction factor')
title('Pixel-wise comparison of fielding maps')
legend('Pixels', 'Identity line', 'Location', 'best')
