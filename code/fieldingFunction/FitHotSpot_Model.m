%% FitHotspot_Model.m

clear; close all; clc

load("planetarium_average_SEC3.mat")   % avg_img, med_img, times_sec, frames_idx, fps

% Use median or average image
% I = double(med_img);
I = double(avg_img);

% Get image height and width
[H, W] = size(I);
[X, Y] = meshgrid(1:W, 1:H);

% INITIAL GUESSES FOR THE FIT:

% baseline = darkest approximate intensity in the image
% amplitude = difference between the brightest and darkest parts 
baseline0 = min(I(:));
amp0 = max(I(:)) - min(I(:));

% Brightest pixel in the image (CENTER OF HOTSPOT) 
[~, maxIdx] = max(I(:));
[y0_guess, x0_guess] = ind2sub(size(I), maxIdx);

% Potential width of hotspot
sigma0 = min(H,W) / 3;

% 2D Gaussian + baseline
% ^^ predicted brightness = baseline + hotspot amplitude × smooth hill centered at (x0, y0)
gauss2D = @(p, x, y) ...
    p(1) + p(2) .* exp( ...
    -((x - p(3)).^2 ./ (2*p(5)^2) + ...
    (y - p(4)).^2 ./ (2*p(6)^2)) );

% p = [baseline, amplitude, x0, y0, sigmaX, sigmaY]
modelFun = @(p, xy) gauss2D(p, xy(:,1), xy(:,2));

xy = [X(:), Y(:)];
z = I(:);

p0 = [baseline0, amp0, x0_guess, y0_guess, sigma0, sigma0];

% Set bounds for the fit
lb = [0, 0, 1, 1, 10, 10];
ub = [Inf, Inf, W, H, W, H];

opts = optimoptions('lsqcurvefit', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5000);

% Fitting
pFit = lsqcurvefit(modelFun, p0, xy, z, lb, ub, opts);

% Use fitted parameters to generate full image-sized model of hotspot
hotspot_fit = reshape(modelFun(pFit, xy), H, W);

% Subtract fitted hotspot from real image
residual = I - hotspot_fit;

%% Plot original image

figure
imagesc(I)
axis image
colorbar
title('Original average image')

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

%% Horizontal profile through hotspot center

% Hotspot center from fitted parameters
x0 = pFit(3);
y0 = pFit(4);

% Use nearest image row through hotspot center
rowIdx = round(y0);

xPixels = 1:W;

dataProfile = I(rowIdx, :);
fitProfile = hotspot_fit(rowIdx, :);
residProfile = residual(rowIdx, :);

figure
plot(xPixels, dataProfile, 'w-', 'LineWidth', 2)
hold on
plot(xPixels, fitProfile, 'r--', 'LineWidth', 2)
plot(xPixels, residProfile, 'b-', 'LineWidth', 1)

xlabel('X Pixel')
ylabel('Pixel Intensity')
title(sprintf('Horizontal profile through hotspot center, y = %d', rowIdx))
legend('Observed average image', 'Gaussian hotspot fit', 'Residual', ...
    'Location', 'best')
grid on

%%
figure
plot(xPixels, residProfile, 'b-', 'LineWidth', 2)

yline(0, 'k--')

xlabel('X Pixel')
ylabel('Residual Intensity')
title(sprintf('Residual profile through hotspot center, y = %d', rowIdx))
grid on