% An analysis of the sensor value non-linearities observed in the IMX219
% camera chip. The goal of this analysis is to transform raw, 8 bit camera
% sensor values to a linearized form, with a physical interpretation of the
% relative sensor values.
%
% Using the light sphere, Zach placed a 0.4 ND filter in the optical path
% (along with the IR filter), set the CombiLED to the half-on background,
% and measured the AGC values the provided a mean image intensity (across
% all classes of pixels) of 127. The camera was then locked to these AGC
% settings. Then, he obtained images of the interior of the sphere with the
% CombiLED set to settings of 0.1:0.1:1. Next, the 0.4 ND filter was
% removed, and the measurement was again made for settings between 0.1 and
% 1 (importantly, with the AGC values still locked to those values obtained
% with the 0.4 ND filter). Zach then processed the data to provide the mean
% R, G, and B channel value at each CombiLED setting level for each of the
% two ND filter conditions (0 and 0.4). These data are provided in the file
% "camera_linearity_ND0_ND0p4_rgb_means.mat".
%
% This routine models the data present in the file. Examination of the
% camera sensor values with increasing light input reveals a soft roll-off
% in the linearity of the response at higher levels. This is not due to a
% traditional "gamma function", as we are saving the data in RAW format
% with no processing. Instead, this likely represents a "full well
% capacity" effect that can happen to photodiodes as they reach saturation.
% After some experimentation, we find that this effect is well-modeled with
% an exponential clipping function that transitions from linear at low
% values to a flat, asymptotic upper limit.
%
% We presume that the nature of this non-linearity is the same for the R,
% G, and B sensor pixels. These pixels derive their spectral sensitivity
% from transmittance filters that are placed before each well. This causes
% variation in the overall sensitivity of the different pixels to the light
% source. We model this effect by estimating the effective neutral density
% filter that is present for the R and G pixels relative to the B pixels.
%
% Running this routine fits the entire dataset with a single non-linear
% function. The exponential ("n") parameter of the fit is the only free
% parameter that is of consequence for the non-linearity. We then call
% a linearization function that transforms raw sensor values to a
% linearized form. In practice, this adjusts the original value range of
% 16-254 to the linearized range of 0-254 (removing the dark value). Raw
% sensor values of >254 are set to 255. Linearization is performed by the
% Python function world_util.linearize_camera_responsivity, which is the
% single source of truth for this operation.
%
% We find a clipping exponent of 5.3918. With linearization, we gain an
% increased ability to represent higher levels of illumination, with some
% reduction in precision for lower levels. This result is saved in the
% top-level "derived" folder for subsequent use in pre-processing of camera
% recordings.

% Housekeeping
clear
close all

% Load the Python world-camera utility library
addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab"));
world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));
numpy = py.importlib.import_module('numpy');

% Load in the calibration data 
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'camera_linearity_ND0_ND0p4_rgb_means.mat');
load(dataFileName)
data = nd04_rgb_means;

% Loop through the R, G, and B sensors and assemble the full measurement
% set.
y = [];
for channel = 1:3
    y = [y; data.ND0p4.rgb_mean(:,channel); ...
        data.ND0.rgb_mean(:,channel)];
end

% Define the objective
myObj = @(p) norm(y-myFitFunc(p));

% Set the initial search parameters
p0 = [4,5.7,0,0,0.4];

% Search
[p,fVal] = fmincon(myObj,p0);

% Create x support that describes multiples of light intensity relative to
% the light intensity that evokes a sensor value of 127.
myObj = @(x) myClippedVal(x,p)-127;
xAbs127 = fzero(myObj,20);
myObj = @(x) myClippedVal(x,p)-254;
xAbs254 = fzero(myObj,100);
[yFit,xAbs] = myFitFunc(p);
xRel = xAbs / xAbs127;

% Create a figure
figure('Position',[1 1 1000 400]);
tiledlayout(1,3)

% Plot the result
nexttile
[~,idxSort] = sort(xRel);
plot(xRel(idxSort),yFit(idxSort),'-','Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold
chanColors = {'r','g','b'};
for channel = 1:3
    chanIdx = (channel-1)*20+1:channel*20;
    plot(xRel(chanIdx),y(chanIdx),['.',chanColors{channel}],'MarkerSize',20);
end
plot([0 1],[127 127],':k');
plot([1 1],[0 127],':k');
ylim([0 255]);
ylabel('8-bit Sensor value');
xlabel('Relative mean luminance [AU]');
legend({'fit',...
    sprintf('R, relative nd = %2.2f',p(3)),...
    sprintf('G, relative nd = %2.2f',p(4)),...
    sprintf('B, relative nd = %2.2f',0)},...
    'Location','southeast');
title(sprintf('Clipping exponent n = %2.4f\n',p(1)));
box off

% Plot the inverse gamma function
nexttile
yLinear = double(world_util.linearize_camera_responsivity(...
    numpy.array(y(idxSort)),...
    pyargs('clipping_exponent',p(1))));
plot(y(idxSort),yLinear,'-k','LineWidth',2);
xlabel('Raw sensor value');
ylabel('Linearized sensor value');
xlim([0 255]);
ylim([0 255]);
a = gca();
a.XTick = 0:50:250;
box off
title('Linearization function');

% Plot the interpretation of luminance vs. linearized sensor value
nexttile
plot(yLinear,xRel(idxSort),'.k','LineWidth',2);
hold on
mdl = fitlm(yLinear,xRel(idxSort), 'y ~ x1 - 1');
b = [mdl.Coefficients.Estimate, 0];
yFit = polyval(b,0:1:255);
plot(0:1:255,yFit,'-r');
[~,idxCen] = min(abs(yFit-1));
plot([idxCen idxCen],[0 1],':k');
plot([0 idxCen],[1 1],':k');
text(idxCen,0.5,sprintf('<-- %d',idxCen));
xlim([0 255]);
xlabel('Linearized sensor value');
ylabel('Relative luminance [AU]');
title('Linear sensor interpretation');
box off
a = gca();
a.XTick = 0:50:250;

% Save the exponent of the fit in the "derived" directory for use in the
% linearization stage, and some other things
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'nonLinearClippingExponent.mat');
clippingExponent = p(1);
linearizedSetPoint = idxCen;
readme = ['Created by defineFullWellCapacityEffect.\n'...
    'clippingExponent -- the exponent of the soft non-linear function that\ndefines the roll-off of sensor values with higher light levels.\n',...
    'linearizedSetPoint -- the sensor value in the linearized values that reflects\nthe set point of the automatic gain control (127 in the original sensor values).\n'];
save(saveFileName,'readme','clippingExponent','linearizedSetPoint');


%% Local function to implement algebraic soft-clipping
function y = myClippedVal(x,p)

% Define the minimum dark value
Smin = 8*2;

% Define the fixed asymptotic value
Smax = 2^8-1 - Smin;

% Unpack the parameters
n = p(1);
a = p(2);

% Calculate
y = Smin + (a .* x) ./ (1 + ((a .* x) ./ Smax).^n).^(1./n);

end


%% Local function to perform the fit
function [yFit,x] = myFitFunc(p)

% Unpack the parameters
ndRed = p(3);
ndGreen = p(4);
ndSet = p(5);

% Assemble the full set of effective ND filters for the chromatic channels,
% expressed relative to Blue.
bayerNDVals = [ndRed, ndGreen, 0];

% Define the range of CombiLED settings that were studied
setVals = 0.1:0.1:1;

% Loop through the channels and define the x domain, accounting for the
% effective ND filter applied to the B and R channels relative to green,
% and the ND filter present for the second measurement set
background = 100;
x = [];
for channel = 1:3
    thisX = [(background * setVals / 10^ndSet),(background * setVals)];
    thisX = thisX / 10^bayerNDVals(channel);
    x = [x, thisX];
end
x = x';

% Calculate the fit. This is an algebraic soft-clipping function
yFit = myClippedVal(x,p);

end
