%% Adaptive gain demo
% Simulates the behavior of camera temporal sensitivity across a range of
% light levels, under the control of an adaptive gain routine.
%
% Recommended settings for our purposes:
%   framesPerUpdate = 50 
%   speedSetting = 0.95
%   settleDurSec = 30
%
% We can update the camera exposure / gain every ~250 msecs without causing
% a loss in the camera frame rate. With a speedSetting of 0.95, we obtain
% fairly flat TTFs across a range of background and contrast levels, with a
% roll-off below 0.5 Hz. It takes up to 30 seconds for the camera to adjust
% to a large change in background light level under these circumstances.
%


clear
close all

% Flag to create time-series plots at each light level
plotAll = true;

% "Speed" parameter
speedSetting = 0.95;

% Properties of the "settle"
settleDurSec = 30;

% Camera fps
fps = 120;

% Update frequency
secsPerUpdate = 0.25;
framesPerUpdate = fps*secsPerUpdate;

% Range of signal values
signalRange = [0,255];

% Properties of the sinusoid
modDurSecs = 10;
f0 = [0.25, 0.5, 1, 3, 6, 12, 25, 50, 100];
contrast = 0.5;
backgroundNDF = [0,-1,-2,-3,-4];

% Frame rate for the plotted fit
fpsModel = 10000;

figHandleTTF = figure();

for bb = 1:length(backgroundNDF)

    if plotAll
        figHandleTime = figure('Name',sprintf('Background 10^%d',backgroundNDF(bb)));
        t = tiledlayout(4,length(f0));
    end

    for ff = 1:length(f0)

        % Create the settle
        settle = ones(1,round(fps*settleDurSec))*10^backgroundNDF(bb);

        % Create the sinusoidal modulation
        modulation = sin(linspace(0,2*pi*(modDurSecs/(1/f0(ff))),round(fps*modDurSecs))) * contrast * 10^backgroundNDF(bb) + 10^backgroundNDF(bb);

        % Assemble the source time series
        source = [settle, modulation];

        % Set the properties of the gain and exposure
        gain = 1;
        exposure = 5000;

        % Clear the signal variable
        signal = nan(size(source));

        % Loop through the source
        for ii = 1:length(source)

            % Obtain the signal (needed for this simulation, in practice, this is
            % the mean value of the sensor array counts)
            s = round(source(ii)*exposure*gain);
            s = min(s,signalRange(2));
            s = max(s,signalRange(1));
            signal(ii) = s;

            % Update the gain and exposure at the update frequency
            if mod(ii,framesPerUpdate)==0
                [gain, exposure] = AGC(s, gain, exposure, speedSetting, fps);
            end

            % Store the gain and exposure values for plotting
            gainStore(ii) = gain;
            exposureStore(ii) = exposure;

        end

        % Define the time domain of the signal
        deltaTsignal = 1/fps;
        tsSignal = 0:deltaTsignal:(length(source)-1)*deltaTsignal;

        % Get the amplitude of the signal at the source frequency
        Y = signal(length(settle)+1:end);
        meanY = mean(signal(length(settle)+1:end));
        Y = (Y - meanY)/meanY;
        [r2(ff,bb),amplitude(ff,bb),~,fit] = fourierRegression( Y, tsSignal(length(settle)+1:end), f0(ff));
        fit = (fit * meanY) + meanY;

        % Define the time domain of the fit
        deltaTmodel = 1/fpsModel;
        offset = (length(settle)+1)*deltaTsignal;
        tsModel = offset:deltaTmodel:(length(fit)-1)*deltaTmodel + offset;

        % Plot this case
        if plotAll
            figure(figHandleTime);
            num = tilenum(t,1,ff);
            nexttile(num);
            plot(tsSignal,log10(source),'-r','LineWidth',1.5);
            ylim([-3.5 0.5]);
            title('source light intensity')
            ylabel('log intensity')

            num = tilenum(t,2,ff);
            nexttile(num);
            plot(tsSignal,exposureStore,'-r','LineWidth',1.5);
            ylim([-500 5500]);
            title('exposure')
            ylabel('exposure [μsecs]')

            num = tilenum(t,3,ff);
            nexttile(num);
            plot(tsSignal,gainStore,'-r','LineWidth',1.5);
            ylim([-0.5 12]);
            title('gain')
            ylabel('gain [a.u.]')

            num = tilenum(t,4,ff);
            nexttile(num);
            plot(tsSignal,signal,'-r','LineWidth',1.5);
            hold on
            plot(tsModel,fit,'-k','LineWidth',1.5);
            ylim([-25 300])
            title('signal')
            xlabel('time [seconds]')
            a = gca();
            a.YTick = [0:50:250];
        end

    end

    % Plot the TTF
    figure(figHandleTTF);
    plot(log10(f0),amplitude(:,bb),'.-','LineWidth',1.5);
    hold on

end

% Clean up the TTF figure
figure(figHandleTTF);
refline(0,contrast)
legendTxt = 'background 10^'+string(backgroundNDF);
legendTxt = [legendTxt 'stimulus contrast'];
legend(legendTxt,'Location','southwest')
xlabel('frequency [Hz]');
ylabel('amplitude [proportion change]');
a = gca();
a.XTick = log10(f0);
a.XTickLabel = string(f0);
ylim([0 contrast*1.05]);



%% Local function that is the AGC
function [gain, exposure] = AGC(s, gain, exposure, speedSetting, fps)

signalTarget = 127;
gainRange = [1 10.666];
exposureRange = [37,floor(1e6/fps)];
signalRange = [0,255];

% Calculate the adjustment
correction = 1+(signalTarget-s)/signalTarget;

% Set the speed
speed = speedSetting;

% Move quickly if we are pegged at the signal range
if s == signalRange(1) || s == signalRange(2)
    speed = speedSetting^3;
end

% Move quickly if we are close to the destination
if abs(correction - 1)<0.25
    speed = speedSetting^2;
end

% Correct the correction
correction = 1 + ((1-speed) * (correction-1));

% If correction > 1, it means we need to turn up gain or exposure.
if correction > 1
    % First choice is to turn up exposure
    if exposure < exposureRange(2)
        exposure = exposure * correction;
        exposure = min([exposure,exposureRange(2)]);
        exposure = max([exposure,exposureRange(1)]);
    else
        gain = gain * correction;
        gain = min([gain,gainRange(2)]);
        gain = max([gain,gainRange(1)]);
    end
end

% If correction < 1, it means we need to turn down gain or exposure.
if correction < 1
    % First choice is to turn down gain
    if gain > gainRange(1)
        gain = gain * correction;
        gain = min([gain,gainRange(2)]);
        gain = max([gain,gainRange(1)]);
    else
        exposure = exposure * correction;
        exposure = min([exposure,exposureRange(2)]);
        exposure = max([exposure,exposureRange(1)]);
    end
end

end