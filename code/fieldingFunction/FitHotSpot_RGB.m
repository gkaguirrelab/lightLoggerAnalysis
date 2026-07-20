%% FitHotspot_Model.m

clear; close all; clc

%% SETTINGS

fitRGB = true;   % false = fit avg_img only; true = fit R/G/B Bayer channels

dataFile = "selected_twilight_frames.mat";
avgFile  = "planetarium_average.mat";

bayerPattern = "BGGR";

%% LOAD DATA

if fitRGB

    load(dataFile)   % frames, times_sec, frames_idx, fps

    frames = double(frames);
    nFrames = size(frames, 1);

    [avg_R, avg_G, avg_B] = makeBayerChannelAverages(frames, bayerPattern);

    imageSet = {avg_R, avg_G, avg_B};
    imageNames = {'R channel', 'G channel', 'B channel'};
    saveName = "hotspot_fit_results_RGB.mat";

else

    load(avgFile)

    imageSet = {double(avg_img)};
    imageNames = {'Original average image'};
    saveName = "hotspot_fit_results_average.mat";

end

%% FIT EACH IMAGE

results = struct();

for ii = 1:numel(imageSet)

    I = double(imageSet{ii});
    imageName = imageNames{ii};

    [pFit, hotspot_fit, residual, X, Y] = fitGaussianHotspot(I);

    results(ii).imageName = imageName;
    results(ii).I = I;
    results(ii).pFit = pFit;
    results(ii).hotspot_fit = hotspot_fit;
    results(ii).residual = residual;

    fprintf('\n%s fitted parameters:\n', imageName)
    fprintf('Baseline: %.3f\n', pFit(1))
    fprintf('Amplitude: %.3f\n', pFit(2))
    fprintf('x0: %.3f\n', pFit(3))
    fprintf('y0: %.3f\n', pFit(4))
    fprintf('sigmaX: %.3f\n', pFit(5))
    fprintf('sigmaY: %.3f\n', pFit(6))
    fprintf('k flattening: %.3f\n', pFit(7))

end

%% PLOT ALL RGB PROFILES TOGETHER

if fitRGB
    plotRGBProfilesAllDirections(results)
else
    plotSingleProfilesAllDirections(results(1))
end

%% COMPARE RGB FIT PARAMETERS

if fitRGB

    channelNames = string({results.imageName});
    sigmaX = arrayfun(@(r) r.pFit(5), results);
    sigmaY = arrayfun(@(r) r.pFit(6), results);
    x0 = arrayfun(@(r) r.pFit(3), results);
    y0 = arrayfun(@(r) r.pFit(4), results);
    k = arrayfun(@(r) r.pFit(7), results);

    fitSummary = table(channelNames', x0', y0', sigmaX', sigmaY', k', ...
        'VariableNames', {'Channel','x0','y0','sigmaX','sigmaY','k'});

    disp(fitSummary)

end

%% SAVE

save(saveName, "results", "fitRGB")

if fitRGB
    save(saveName, "results", "fitRGB", "avg_R", "avg_G", "avg_B", ...
        "times_sec", "frames_idx", "fps", "-append")
end

%% LOCAL FUNCTIONS

function [avg_R, avg_G, avg_B] = makeBayerChannelAverages(frames, bayerPattern)

    nFrames = size(frames, 1);

    for k = 1:nFrames

        rawFrame = squeeze(frames(k,:,:));

        [R, G, B] = splitBayerFrame(rawFrame, bayerPattern);

        if k == 1
            sum_R = zeros(size(R));
            sum_G = zeros(size(G));
            sum_B = zeros(size(B));
        end

        sum_R = sum_R + R;
        sum_G = sum_G + G;
        sum_B = sum_B + B;

    end

    avg_R = sum_R ./ nFrames;
    avg_G = sum_G ./ nFrames;
    avg_B = sum_B ./ nFrames;

end

function [R, G, B] = splitBayerFrame(rawFrame, bayerPattern)

    rawFrame = double(rawFrame);

    switch upper(bayerPattern)

        case "BGGR"
            B  = rawFrame(1:2:end, 1:2:end);
            G1 = rawFrame(1:2:end, 2:2:end);
            G2 = rawFrame(2:2:end, 1:2:end);
            R  = rawFrame(2:2:end, 2:2:end);

        case "RGGB"
            R  = rawFrame(1:2:end, 1:2:end);
            G1 = rawFrame(1:2:end, 2:2:end);
            G2 = rawFrame(2:2:end, 1:2:end);
            B  = rawFrame(2:2:end, 2:2:end);

        case "GRBG"
            G1 = rawFrame(1:2:end, 1:2:end);
            R  = rawFrame(1:2:end, 2:2:end);
            B  = rawFrame(2:2:end, 1:2:end);
            G2 = rawFrame(2:2:end, 2:2:end);

        case "GBRG"
            G1 = rawFrame(1:2:end, 1:2:end);
            B  = rawFrame(1:2:end, 2:2:end);
            R  = rawFrame(2:2:end, 1:2:end);
            G2 = rawFrame(2:2:end, 2:2:end);

        otherwise
            error("Unknown Bayer pattern: %s", bayerPattern)

    end

    G = (G1 + G2) ./ 2;

end

function [pFit, hotspot_fit, residual, X, Y] = fitGaussianHotspot(I)

    [H, W] = size(I);
    [X, Y] = meshgrid(1:W, 1:H);

    baseline0 = min(I(:));
    amp0 = max(I(:)) - min(I(:));

    [~, maxIdx] = max(I(:));
    [y0_guess, x0_guess] = ind2sub(size(I), maxIdx);

    sigma0 = min(H,W) / 3;
    k0 = 1;

    flatGauss2D = @(p, x, y) ...
        p(1) + p(2) .* tanh( ...
        p(7) .* exp( ...
        -((x - p(3)).^2 ./ (2*p(5)^2) + ...
          (y - p(4)).^2 ./ (2*p(6)^2)) ) );

    modelFun = @(p, xy) flatGauss2D(p, xy(:,1), xy(:,2));

    xy = [X(:), Y(:)];
    z = I(:);

    p0 = [baseline0, amp0, x0_guess, y0_guess, sigma0, sigma0, k0];

    lb = [0, 0, 1, 1, 10, 10, 0.01];
    ub = [Inf, Inf, W, H, W, H, 100];

    opts = optimoptions('lsqcurvefit', ...
        'Display','iter', ...
        'MaxFunctionEvaluations', 5000);

    pFit = lsqcurvefit(modelFun, p0, xy, z, lb, ub, opts);

    hotspot_fit = reshape(modelFun(pFit, xy), H, W);
    residual = I - hotspot_fit;

end

function plotRGBProfilesAllDirections(results)

    colors = {'r','g','b'};
    labels = {'R','G','B'};
    directions = {'horizontal','vertical','diagDown','diagUp'};
    directionTitles = {'Horizontal','Vertical','Diagonal down','Diagonal up'};

    for dd = 1:numel(directions)

        direction = directions{dd};

        figure
        hold on

        legendEntries = {};

        for ii = 1:3

            I = results(ii).I;
            hotspot_fit = results(ii).hotspot_fit;
            pFit = results(ii).pFit;

            % Normalize using min/max of fitted 2D image
            fitMin = min(hotspot_fit(:));
            fitMax = max(hotspot_fit(:));

            I_norm = (I - fitMin) ./ (fitMax - fitMin);
            fit_norm = (hotspot_fit - fitMin) ./ (fitMax - fitMin);

            [axisVals, dataProfile] = getProfile(I_norm, pFit, direction);
            [~, fitProfile] = getProfile(fit_norm, pFit, direction);

            plot(axisVals, dataProfile, ...
                'Color', colors{ii}, ...
                'LineStyle', '-', ...
                'LineWidth', 1.5)

            plot(axisVals, fitProfile, ...
                'Color', colors{ii}, ...
                'LineStyle', '--', ...
                'LineWidth', 2)

            legendEntries{end+1} = sprintf('%s observed', labels{ii});
            legendEntries{end+1} = sprintf('%s fit', labels{ii});

        end

        xlabel('Position along profile')
        ylabel('Normalized Intensity')
        title(sprintf('%s profile through hotspot center: normalized observed and fit', directionTitles{dd}))
        legend(legendEntries, 'Location', 'best')
        grid on

        figure
        hold on

        for ii = 1:3

            I = results(ii).I;
            hotspot_fit = results(ii).hotspot_fit;
            pFit = results(ii).pFit;

            fitMin = min(hotspot_fit(:));
            fitMax = max(hotspot_fit(:));

            I_norm = (I - fitMin) ./ (fitMax - fitMin);
            fit_norm = (hotspot_fit - fitMin) ./ (fitMax - fitMin);

            residual_norm = I_norm - fit_norm;

            [axisVals, residProfile] = getProfile(residual_norm, pFit, direction);

            plot(axisVals, residProfile, ...
                'Color', colors{ii}, ...
                'LineStyle', '-', ...
                'LineWidth', 2)

        end

        yline(0, 'k--')

        xlabel('Position along profile')
        ylabel('Normalized Residual')
        title(sprintf('%s profile through hotspot center: normalized residuals', directionTitles{dd}))
        legend({'R residual','G residual','B residual'}, 'Location', 'best')
        grid on

    end

end

function plotSingleProfilesAllDirections(result)

    directions = {'horizontal','vertical','diagDown','diagUp'};
    directionTitles = {'Horizontal','Vertical','Diagonal down','Diagonal up'};

    for dd = 1:numel(directions)

        direction = directions{dd};

        I = result.I;
        hotspot_fit = result.hotspot_fit;
        residual = result.residual;
        pFit = result.pFit;

        [axisVals, dataProfile] = getProfile(I, pFit, direction);
        [~, fitProfile] = getProfile(hotspot_fit, pFit, direction);
        [~, residProfile] = getProfile(residual, pFit, direction);

        figure
        plot(axisVals, dataProfile, 'k-', 'LineWidth', 2)
        hold on
        plot(axisVals, fitProfile, 'r--', 'LineWidth', 2)
        plot(axisVals, residProfile, 'b-', 'LineWidth', 1)

        xlabel('Position along profile')
        ylabel('Pixel Intensity')
        title(sprintf('%s profile through hotspot center: observed, fit, residual', directionTitles{dd}))
        legend('Observed image','Flattened Gaussian fit','Residual', 'Location','best')
        grid on

        figure
        plot(axisVals, residProfile, 'b-', 'LineWidth', 2)
        yline(0, 'k--')

        xlabel('Position along profile')
        ylabel('Residual Intensity')
        title(sprintf('%s profile through hotspot center: residual', directionTitles{dd}))
        grid on

    end

end

function [axisVals, profile] = getProfile(I, pFit, direction)

    [H, W] = size(I);

    x0 = pFit(3);
    y0 = pFit(4);

    switch direction

        case 'horizontal'

            rowIdx = round(y0);
            xq = 1:W;
            yq = rowIdx .* ones(size(xq));
            axisVals = xq;

        case 'vertical'

            colIdx = round(x0);
            yq = 1:H;
            xq = colIdx .* ones(size(yq));
            axisVals = yq;

        case 'diagDown'

            tMin = ceil(max(1 - x0, 1 - y0));
            tMax = floor(min(W - x0, H - y0));
            t = tMin:tMax;

            xq = x0 + t;
            yq = y0 + t;
            axisVals = t;

        case 'diagUp'

            tMin = ceil(max(1 - x0, y0 - H));
            tMax = floor(min(W - x0, y0 - 1));
            t = tMin:tMax;

            xq = x0 + t;
            yq = y0 - t;
            axisVals = t;

        otherwise
            error('Unknown direction: %s', direction)

    end

    [X, Y] = meshgrid(1:W, 1:H);

    profile = interp2(X, Y, I, xq, yq, 'linear');

end