%% BuildCorrectionMap.m

clear; close all; clc

load("hotspot_fit_results_RGB.mat")   % loads results, avg_R, avg_G, avg_B, etc.

% Pull fitted hotspot images for each Bayer channel
fit_R = results(1).hotspot_fit;
fit_G = results(2).hotspot_fit;
fit_B = results(3).hotspot_fit;

% Compute correction maps for each reduced-size Bayer channel
corr_R = max(fit_R(:)) ./ fit_R;
corr_G = max(fit_G(:)) ./ fit_G;
corr_B = max(fit_B(:)) ./ fit_B;

%% Verify correction maps flatten the fitted hotspot

correctedFit_R = fit_R .* corr_R;
correctedFit_G = fit_G .* corr_G;
correctedFit_B = fit_B .* corr_B;

fprintf('\n----- CORRECTION MAP VERIFICATION -----\n')

fprintf('\nR channel:\n')
fprintf('Original fit range: %.3f to %.3f\n', ...
    min(fit_R(:)), max(fit_R(:)))
fprintf('Corrected fit range: %.3f to %.3f\n', ...
    min(correctedFit_R(:)), max(correctedFit_R(:)))
fprintf('Corrected fit std: %.6f\n', ...
    std(correctedFit_R(:)))

fprintf('\nG channel:\n')
fprintf('Original fit range: %.3f to %.3f\n', ...
    min(fit_G(:)), max(fit_G(:)))
fprintf('Corrected fit range: %.3f to %.3f\n', ...
    min(correctedFit_G(:)), max(correctedFit_G(:)))
fprintf('Corrected fit std: %.6f\n', ...
    std(correctedFit_G(:)))

fprintf('\nB channel:\n')
fprintf('Original fit range: %.3f to %.3f\n', ...
    min(fit_B(:)), max(fit_B(:)))
fprintf('Corrected fit range: %.3f to %.3f\n', ...
    min(correctedFit_B(:)), max(correctedFit_B(:)))
fprintf('Corrected fit std: %.6f\n', ...
    std(correctedFit_B(:)))

%% Visual confirmation

figure
surf(correctedFit_R,'EdgeColor','none')
view(3)
colorbar
title('R channel after applying correction map')

figure
surf(correctedFit_G,'EdgeColor','none')
view(3)
colorbar
title('G channel after applying correction map')

figure
surf(correctedFit_B,'EdgeColor','none')
view(3)
colorbar
title('B channel after applying correction map')

%% Sanity checks

channelNames = {'R','G','B'};
fits = {fit_R, fit_G, fit_B};
corrs = {corr_R, corr_G, corr_B};

for ii = 1:3
    fit = fits{ii};
    corr = corrs{ii};

    peakFit = max(fit(:));
    minFit = min(fit(:));

    fprintf('\n%s channel:\n', channelNames{ii})
    fprintf('Peak fitted value: %.4f\n', peakFit)
    fprintf('Min fitted value: %.4f\n', minFit)
    fprintf('Correction at peak should be 1: %.4f\n', min(corr(:)))
    fprintf('Max correction: %.4f\n', max(corr(:)))
    fprintf('Expected max correction = peak/min: %.4f\n', peakFit/minFit)

    % Verify that fit .* correction gives a flat image at peakFit
    correctedFit = fit .* corr;
    fprintf('Corrected fit min: %.4f\n', min(correctedFit(:)))
    fprintf('Corrected fit max: %.4f\n', max(correctedFit(:)))
    fprintf('Corrected fit SD: %.8f\n', std(correctedFit(:)))
end

% Recombine into full-size Bayer correction map
[Hsmall, Wsmall] = size(corr_G);

Hfull = Hsmall * 2;
Wfull = Wsmall * 2;

correctionMap = nan(Hfull, Wfull);

% Must match the Bayer pattern used in fitting: BGGR
correctionMap(1:2:end, 1:2:end) = corr_B;
correctionMap(1:2:end, 2:2:end) = corr_G;
correctionMap(2:2:end, 1:2:end) = corr_G;
correctionMap(2:2:end, 2:2:end) = corr_R;

%% Plot correction map image

figure
imagesc(correctionMap)
axis image
colorbar
title('Full Bayer correction map')

xlabel('X Pixel')
ylabel('Y Pixel')

%% Plot correction map surface

figure
surf(correctionMap, 'EdgeColor', 'none')
view(3)
colorbar

xlabel('X Pixel')
ylabel('Y Pixel')
zlabel('Correction Factor')

title('Full Bayer correction map surface')

%% Plot each channel correction map separately

figure
imagesc(corr_R)
axis image
colorbar
title('R channel correction map')

figure
imagesc(corr_G)
axis image
colorbar
title('G channel correction map')

figure
imagesc(corr_B)
axis image
colorbar
title('B channel correction map')

%% Save final correction map

save("fielding_correction_map_RGB.mat", ...
    "correctionMap", ...
    "corr_R", "corr_G", "corr_B", ...
    "fit_R", "fit_G", "fit_B")