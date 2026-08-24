% The data we analyze here began as indoor and outdoor recordings from 17
% subjects engaged in 10 activities. The current linear-scale data file was
% generated from all 111 processed subject/activity recordings after removing
% nonfinite and nonpositive AGC product / illuminance samples, removing frames
% with >40% spatial saturation, and excluding the first 100 cached samples
% from every recording. The subjects and activities represented in the
% exported point cloud are:
%
%   Subjects: FLIC_20, FLIC_21, FLIC_28, FLIC_39, FLIC_42, FLIC_51,
%             FLIC_1029, FLIC_1034, FLIC_1038, FLIC_1044, FLIC_1047,
%             FLIC_2001, FLIC_2002, FLIC_2003, FLIC_2004, FLIC_2005,
%             FLIC_2006
%   Activities: chat, gazeCalibration, lunch, phone, read, sitBiopond,
%               walkBiopond, walkIndoor, walkOutdoor, work
%
% The Python diagnostic dashboard also marks this correlation-qualified
% subset of 20 subject/activity recordings (shared empirical-kernel model
% correlation >= 0.9):
%
%   FLIC_20 read, walkBiopond
%   FLIC_21 read
%   FLIC_28 read
%   FLIC_42 read
%   FLIC_51 read
%   FLIC_1029 read
%   FLIC_1034 read, sitBiopond
%   FLIC_1038 read, sitBiopond
%   FLIC_1044 read
%   FLIC_1047 read
%   FLIC_2002 gazeCalibration
%   FLIC_2003 work
%   FLIC_2004 work
%   FLIC_2005 gazeCalibration, walkOutdoor, work
%   FLIC_2006 work
%
% Each recording includes the AGC camera settings at each time point, and
% the illuminance of the environment as measured by the minispect. The
% product of the AGC camera settings provides the "camera score" which
% measures the light sensitivity of the camera at each point in time. In
% separate code, we estimated the properties of a low-pass temporal filter
% that transforms the illuminance values into the camera score. This
% previous stage of analysis is run from
% deriveEmpircalAGCAndIlluminance.py.

% The data that we load here corresponds to the camera score (x) and the
% temporally filtered illuminance values (y, in units of lux).

% We observe that these values have a power law relationship, such that
% there is a locally linear relationship between the values in a log-log
% plot. There are, however, two domains of this linear relationship, with
% different slopes under lower and higher light conditions. We don't have
% a mechanistic explanation for this. The output of this routine is the
% equation and parameters that may be used to convert current camera score
% into implied environmental illuminance.

% We also have available the steady-state camera AGC values obtained when
% the camera was directed into the table-top light integrating sphere and
% exposed the max settings of the combiLED in the presence of different
% neutral density filters. Zach provided these values which are hard-coded
% below and then added to the plot. There is generally good agreement in
% the form of these data.

% Run the shared fit and retrieve the cleaned environmental calibration data.
analysisDir = fileparts(mfilename('fullpath'));
addpath(analysisDir);
[~, best_p, sensitivity, illuminance] = fitCameraAGCToIlluminance([]);

% Transform the data into log10 space
log_x = log10(sensitivity);
log_y = log10(illuminance);

% Define the piecewise continuous model function
% p(1) = m1 (slope 1), p(2) = c1 (intercept 1)
% p(3) = m2 (slope 2), p(4) = xb (breakpoint in log10 space)
piecewise_model = @(p, x) (x < p(4)) .* (p(1) .* x + p(2)) + ...
                          (x >= p(4)) .* (p(1) * p(4) + p(2) + p(3) .* (x - p(4)));

% Extract optimized parameters
m1 = best_p(1); c1 = best_p(2); m2 = best_p(3); xb = best_p(4);

% Report the function that may be used to derive illuminance from the
% camera score
fprintf('log 10 illuminance (lx) = [\n');
fprintf('                   x < %2.2f | %2.2f * x + %2.2f;\n',xb,m1,c1 );
fprintf('                   x >= %2.2f | (%2.2f * %2.2f) + %2.2f + %2.2f * (x - %2.2f);\n',xb,m1,xb,c1,m2,xb );
fprintf('\n');
fprintf('where x is log10 of the product of the camera AGC values.\n');

% Create a function handle for easy read-out in linear space
predict_illuminance = @(sens) 10.^(piecewise_model(best_p, log10(sens)));

% Plotting
figure;
loglog(sensitivity, illuminance, '.','Color',[0.75 0.75 0.75], 'MarkerSize', 5); 
hold on; grid on;

% Generate points for the fit curve across the range
x_fit_curve = logspace(min(log_x), max(log_x), 500);
y_fit_curve = predict_illuminance(x_fit_curve);

% Plot the fitted piecewise line over the original data
loglog(x_fit_curve, y_fit_curve, 'r-', 'LineWidth', 2);

% Mark the breakpoint on the plot
plot(10^xb, predict_illuminance(10^xb), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

% Add the camera scores and illuminance values obtained within the desktop
% light integrating sphere. The camera was directed into the light
% integrating sphere and presented the "max spectrum" of the combiLED after
% passing through different ND filters. The first column is the product of
% the AGC settings at steady state, and the second column is the
% illuminance, obtained by multiplying the measured luminance of the sphere
% interior for that ND filter by pi.
xy = [
1466	54114.26474
15023	6763.340616
140569	666.1113277
530018	73.6074537
708871	7.2570729];
plot(xy(:,1),xy(:,2),'-*k','LineWidth',2,'MarkerSize',10)

% Labels
xlabel('Camera Sensitivity Score (X Data)');
ylabel('Environmental Illuminance (Y Data)');
legend('Environmental Data', 'Robust Piecewise Fit (L1 Norm)', 'Breakpoint','Integrating sphere measures');
title('Robust Log-Log Calibration');
