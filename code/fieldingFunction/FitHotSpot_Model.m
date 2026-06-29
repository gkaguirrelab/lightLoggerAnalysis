%% FitHotspot_Model.m

clear; close all; clc

load("planetarium_average_SEC3.mat")   % avg_img, med_img, times_sec, frames_idx, fps

I = double(med_img);   % use median first

[H, W] = size(I);
[X, Y] = meshgrid(1:W, 1:H);

% Initial guesses
baseline0 = min(I(:));
amp0 = max(I(:)) - min(I(:));

[~, maxIdx] = max(I(:));
[y0_guess, x0_guess] = ind2sub(size(I), maxIdx);

sigma0 = min(H,W) / 3;

% 2D Gaussian + baseline
gauss2D = @(p, x, y) ...
    p(1) + p(2) .* exp( ...
    -((x - p(3)).^2 ./ (2*p(5)^2) + ...
    (y - p(4)).^2 ./ (2*p(6)^2)) );

% p = [baseline, amplitude, x0, y0, sigmaX, sigmaY]
modelFun = @(p, xy) gauss2D(p, xy(:,1), xy(:,2));

xy = [X(:), Y(:)];
z = I(:);

p0 = [baseline0, amp0, x0_guess, y0_guess, sigma0, sigma0];

lb = [0, 0, 1, 1, 10, 10];
ub = [Inf, Inf, W, H, W, H];

opts = optimoptions('lsqcurvefit', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5000);

pFit = lsqcurvefit(modelFun, p0, xy, z, lb, ub, opts);

hotspot_fit = reshape(modelFun(pFit, xy), H, W);
residual = I - hotspot_fit;

%% Plot original image

figure
imagesc(I)
axis image
colorbar
title('Original median image')

%% Plot fitted hotspot

figure
imagesc(hotspot_fit)
axis image
colorbar
title('Fitted hotspot model')

figure
surf(hotspot_fit, 'EdgeColor', 'none')
view(3)
colorbar
zlabel('Pixel Intensity')
title('Fitted hotspot surface')

%% Plot residual

figure
imagesc(residual)
axis image
colorbar
title('Residual after subtracting hotspot')

figure
surf(residual, 'EdgeColor', 'none')
view(3)
colorbar
zlabel('Residual Intensity')
title('Residual surface after subtracting hotspot')

%% Print fit parameters

fprintf('\nFitted parameters:\n')
fprintf('Baseline: %.3f\n', pFit(1))
fprintf('Amplitude: %.3f\n', pFit(2))
fprintf('x0: %.3f\n', pFit(3))
fprintf('y0: %.3f\n', pFit(4))
fprintf('sigmaX: %.3f\n', pFit(5))
fprintf('sigmaY: %.3f\n', pFit(6))

save("hotspot_fit_results.mat", ...
    "I", "hotspot_fit", "residual", "pFit")