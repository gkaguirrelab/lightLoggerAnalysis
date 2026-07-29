%% analyzeAGCEmpiricalKernel.m
%
% Uses the AGC logic from adaptiveControlDemo.m to characterize the
% temporal response of the camera AGC and derive an empirical temporal
% kernel from the simulated step response.
%
% The empirical kernel can later be applied to the minispect signal in
% Python and compared against the current single-exponential filter.
%

clear; close all; clc


%% SETTINGS

speedSetting = 0.95;

fps = 120;

secsPerUpdate = 0.25;
framesPerUpdate = fps * secsPerUpdate;

settleDurSec = 30;
responseDurSec = 40;

signalRange = [0 255];


% Test both directions
conditions = {
    'Dark to bright', 0.5, 1.0;
    'Bright to dark', 1.0, 0.5
    };


results = struct();


%% LOOP THROUGH CONDITIONS

for cc = 1:size(conditions,1)

    conditionName = conditions{cc,1};
    sourceBefore = conditions{cc,2};
    sourceAfter = conditions{cc,3};

    fprintf('\n====================================\n')
    fprintf('%s\n', conditionName)
    fprintf('====================================\n')


    %% CREATE SOURCE LIGHT STEP

    settle = ones(1, round(fps*settleDurSec)) * sourceBefore;

    responsePeriod = ones(1, round(fps*responseDurSec)) * sourceAfter;

    source = [settle responsePeriod];


    %% INITIAL AGC STATE
    %
    % Same initial values as adaptiveControlDemo.m

    gain = 1;
    exposure = 37;

    signal = nan(size(source));
    gainStore = nan(size(source));
    exposureStore = nan(size(source));


    %% PASS SOURCE THROUGH AGC

    for ii = 1:length(source)

        s = round(source(ii) * exposure * gain);

        s = min(s, signalRange(2));
        s = max(s, signalRange(1));

        signal(ii) = s;


        if mod(ii, framesPerUpdate) == 0
            [gain, exposure] = AGC(s, gain, exposure, speedSetting, fps);
        end


        gainStore(ii) = gain;
        exposureStore(ii) = exposure;

    end


    %% CAMERA SCORE

    cameraScore = 1 ./ (gainStore .* exposureStore);

    time = (0:length(source)-1) ./ fps;


    %% EXTRACT RESPONSE AFTER LIGHT STEP

    stepIdx = time >= settleDurSec;

    responseTime = time(stepIdx) - settleDurSec;

    responseScore = cameraScore(stepIdx);


    %% NORMALIZE STEP RESPONSE
    %
    % Normalize so the response begins at 0 and ends at 1.
    % This makes both dark->bright and bright->dark comparable.

    responseNorm = ...
        (responseScore - responseScore(1)) ./ ...
        (responseScore(end) - responseScore(1));


    %% FIT SINGLE EXPONENTIAL STEP RESPONSE
    %
    % This is the same general response shape implied by the
    % first-order low-pass model currently used in Python.

    singleExpStep = @(tau,t) 1 - exp(-t ./ tau);

    tau0 = 5;

    tauFit = lsqcurvefit(singleExpStep, tau0, responseTime, ...
        responseNorm, 0.01, 100);

    expStepFit = singleExpStep(tauFit, responseTime);

    stepResidual = responseNorm - expStepFit;

    stepRMSE = sqrt(mean(stepResidual.^2));


    fprintf('Best exponential tau: %.3f s\n', tauFit)
    fprintf('Step-response RMSE: %.6f\n', stepRMSE)


    %% DERIVE EMPIRICAL AGC KERNEL AT THE ACTUAL AGC UPDATE RATE
    %
    % The camera runs at 120 Hz, but the AGC only updates every 0.25 s.
    % Therefore, characterize the temporal kernel at the AGC update rate
    % rather than differentiating the 120-Hz staircase response.

    updateIdx = 1:framesPerUpdate:length(responseTime);

    kernelTime = responseTime(updateIdx);
    responseAtUpdates = responseNorm(updateIdx);

    % Derivative of the step response gives the empirical impulse response.
    empiricalKernel = gradient(responseAtUpdates, kernelTime);

    % Normalize total kernel area to 1.
    kernelArea = trapz(kernelTime, empiricalKernel);

    if kernelArea ~= 0
        empiricalKernel = empiricalKernel ./ kernelArea;
    end


    %% EXPONENTIAL KERNEL FOR COMPARISON
    %
    % Impulse response associated with the fitted first-order exponential.

    expKernel = (1 ./ tauFit) .* exp(-kernelTime ./ tauFit);
    expKernel = expKernel ./ trapz(kernelTime, expKernel);


    %% SAVE RESULTS

    results(cc).condition = conditionName;

    results(cc).sourceBefore = sourceBefore;
    results(cc).sourceAfter = sourceAfter;

    results(cc).responseTime = responseTime;
    results(cc).responseNorm = responseNorm;

    results(cc).tauFit = tauFit;
    results(cc).expStepFit = expStepFit;
    results(cc).stepResidual = stepResidual;
    results(cc).stepRMSE = stepRMSE;

    results(cc).kernelTime = kernelTime;
    results(cc).empiricalKernel = empiricalKernel;
    results(cc).expKernel = expKernel;

    results(cc).gain = gainStore;
    results(cc).exposure = exposureStore;
    results(cc).cameraScore = cameraScore;
    results(cc).source = source;
    results(cc).time = time;


    %% PLOT 1: SOURCE STEP

    figure

    plot(time, source, 'k-', 'LineWidth', 2)

    xlabel('Time (s)')
    ylabel('Source intensity')

    title(sprintf('%s: input light step', conditionName))

    xline(settleDurSec, 'k--')

    grid on


    %% PLOT 2: CAMERA SCORE STEP RESPONSE

    figure

    plot(responseTime, responseNorm, 'k-', 'LineWidth', 2)
    hold on
    plot(responseTime, expStepFit, 'b--', 'LineWidth', 2)

    xlabel('Time after light step (s)')
    ylabel('Normalized camera-score response')

    title(sprintf('%s: AGC temporal response', conditionName))

    legend( ...
        'Actual AGC simulation', ...
        sprintf('Single exponential, tau = %.2f s', tauFit), ...
        'Location', 'best' ...
        )

    grid on


    %% PLOT 3: STEP-RESPONSE RESIDUAL

    figure

    plot(responseTime, stepResidual, 'r-', 'LineWidth', 1.5)
    hold on

    yline(0, 'k--')

    xlabel('Time after light step (s)')
    ylabel('Residual')

    title(sprintf('%s: exponential-model residual', conditionName))

    grid on


    %% PLOT 4: EMPIRICAL KERNEL

    figure

    plot(kernelTime, empiricalKernel, 'k-', 'LineWidth', 2)

    xlabel('Time after impulse (s)')
    ylabel('Kernel weight')

    title(sprintf('%s: empirical AGC kernel', conditionName))

    grid on


    %% PLOT 5: EMPIRICAL VS EXPONENTIAL KERNEL

    figure

    plot(kernelTime, empiricalKernel, 'k-', 'LineWidth', 2)
    hold on

    plot(kernelTime, expKernel, 'b--', 'LineWidth', 2)

    xlabel('Time after impulse (s)')
    ylabel('Kernel weight')

    title(sprintf('%s: empirical vs exponential kernel', conditionName))

    legend( ...
        'AGC-derived empirical kernel', ...
        sprintf('Exponential kernel, tau = %.2f s', tauFit), ...
        'Location', 'best' ...
        )

    grid on

end


%% SUMMARY

fprintf('\n\n====================================\n')
fprintf('SUMMARY\n')
fprintf('====================================\n')

for cc = 1:length(results)

    fprintf('\n%s\n', results(cc).condition)

    fprintf('Best exponential tau: %.3f s\n', results(cc).tauFit)

    fprintf('Step-response RMSE: %.6f\n', results(cc).stepRMSE)

    fprintf('Empirical kernel area: %.6f\n', ...
        trapz(results(cc).kernelTime, results(cc).empiricalKernel))

end


%% OPTIONAL: COMPARE BOTH EMPIRICAL KERNELS DIRECTLY

figure

plot(results(1).kernelTime, results(1).empiricalKernel, ...
    'LineWidth', 2)
hold on

plot(results(2).kernelTime, results(2).empiricalKernel, ...
    'LineWidth', 2)

xlabel('Time after impulse (s)')
ylabel('Kernel weight')

title('AGC empirical kernels: direction comparison')

legend( ...
    results(1).condition, ...
    results(2).condition, ...
    'Location', 'best' ...
    )

grid on


%% CREATE MEAN EMPIRICAL KERNEL
%
% Dark->bright and bright->dark responses are similar, so create one
% provisional average kernel for the first Python comparison.

commonKernelTime = results(1).kernelTime;

kernelDarkToBright = results(1).empiricalKernel;
kernelBrightToDark = interp1( ...
    results(2).kernelTime, ...
    results(2).empiricalKernel, ...
    commonKernelTime, ...
    'linear', ...
    'extrap');

meanKernel = mean([kernelDarkToBright; kernelBrightToDark], 1);

% Re-normalize area to exactly 1.
meanKernel = meanKernel ./ trapz(commonKernelTime, meanKernel);

figure
plot(commonKernelTime, kernelDarkToBright, 'LineWidth', 1.5)
hold on
plot(commonKernelTime, kernelBrightToDark, 'LineWidth', 1.5)
plot(commonKernelTime, meanKernel, 'k-', 'LineWidth', 2.5)

xlabel('Time after impulse (s)')
ylabel('Kernel weight')
title('AGC empirical kernels at 4 Hz')

legend( ...
    'Dark to bright', ...
    'Bright to dark', ...
    'Mean kernel', ...
    'Location', 'best')

grid on


%% SAVE FOR PYTHON

save('agc_empirical_kernels.mat', ...
    'results', ...
    'commonKernelTime', ...
    'meanKernel')

fprintf('\nSaved agc_empirical_kernels.mat\n')


%% LOCAL FUNCTION
% Exact AGC logic from adaptiveControlDemo.m

function [gain, exposure] = AGC(s, gain, exposure, speedSetting, fps)

signalTarget = 127;

gainRange = [1 10.666];

exposureRange = [37 floor(1e6/fps)];

signalRange = [0 255];


% Calculate the adjustment

correction = 1 + (signalTarget-s)/signalTarget;


% Set the speed

speed = speedSetting;


% Move quickly if pegged at sensor limits

if s == signalRange(1) || s == signalRange(2)
    speed = speedSetting^3;
end


% Move quickly if close to destination

if abs(correction - 1) < 0.25
    speed = speedSetting^2;
end


% Correct the correction

correction = 1 + ((1-speed) * (correction-1));


% Need to increase sensitivity

if correction > 1

    if exposure < exposureRange(2)

        exposure = exposure * correction;

        exposure = min(exposure, exposureRange(2));

        exposure = max(exposure, exposureRange(1));

    else

        gain = gain * correction;

        gain = min(gain, gainRange(2));

        gain = max(gain, gainRange(1));

    end

end


% Need to decrease sensitivity

if correction < 1

    if gain > gainRange(1)

        gain = gain * correction;

        gain = min(gain, gainRange(2));

        gain = max(gain, gainRange(1));

    else

        exposure = exposure * correction;

        exposure = min(exposure, exposureRange(2));

        exposure = max(exposure, exposureRange(1));

    end

end

end