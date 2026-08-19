% The purpose of this script is to define the relationship between the
% custom AGC settings we use to control the sensitivity of the IMX219
% camera, the illuminance to which the camera is exposed, and the mean
% luminance of each pixel in the camera array.
clear
close all

% Load the full-precision fixed AGC settings used for the 0.25 contrast
% calibration from the shared Python world-camera utility module. Each
% Python tuple is ordered as analog gain, digital gain, and exposure.
% These values should remain:
%   NDF       = [0,          1,          2,          3,          4,          5]
%   AGain     = [1,          1.80281687, 10.666,     10.666,     10.666,     10.666]
%   DGain     = [1,          1,          1.58156674, 5.96331605, 7.97561560, 10]
%   Exposure  = [1466,       8333,       8333,       8333,       8333,       8333]
%   AGC score = [1466,       15022.87298, 140569.30074, 530018.20667, 708870.94394, 888797.78]
addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab"));
world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));
world_agc_target = double(py.getattr(world_util, "WORLD_AGC_DEFAULT_TARGET"));
get_ndf_settings = py.getattr(world_util, "get_world_contrast_level_ndf_settings");
world_settings_by_ndf = get_ndf_settings(0.25, world_agc_target);
get_ndf_setting = py.getattr(world_settings_by_ndf, "__getitem__");

agcData.ndf = 0:5;
agcData.AGain = zeros(size(agcData.ndf));
agcData.DGain = zeros(size(agcData.ndf));
agcData.Exposure = zeros(size(agcData.ndf));
for ii = 1:numel(agcData.ndf)
    fixed_settings = double(get_ndf_setting(py.float(agcData.ndf(ii))));
    agcData.AGain(ii) = fixed_settings(1);
    agcData.DGain(ii) = fixed_settings(2);
    agcData.Exposure(ii) = fixed_settings(3);
end

% Derive a "camera score" by obtaining the product of the AGC settings
cameraScore = agcData.DGain .* agcData.AGain .* agcData.Exposure;

% Load the XYZ fundamentals (these live in the PsychToolbox, which is
% placed on the matlab path when we invoke the CombiLEDToolbox)
load('T_xyz1931.mat','T_xyz1931','S_xyz1931');

% Next, load the "maxSpectrum" calibration files and extract for each the
% luminance of the sphere interior. The combiLED was set to the half-on
% settings, so we account for this as well.
settings = 0.5;
calDataFolder = getpref('lightLoggerAnalysis','CalDataFolder');
for ii = 1:length(agcData.ndf)
    thisCalFile = sprintf('CombiLED-A_cassette-ND%d_sphere_maxSpectrum.mat',agcData.ndf(ii));
    load(thisCalFile,'cals');
    cal = cals{end};
    S = cal.rawData.S;
    spd = cal.rawData.gammaCurveMeanMeasurements;
    T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
    xyYLocus = XYZToxyY(T_xyz);
    sphereAvgLuminance(ii) = settings*T_xyz(2,:)*cal.rawData.gammaCurveMeanMeasurements;
end

% The data we analyze here began as indoor and outdoor recordings from 17
% subjects engaged in 10 activities. The current linear-scale data file was
% generated from all 111 processed subject/activity recordings after
% removing nonfinite and nonpositive AGC product / illuminance samples,
% removing frames with >40% spatial saturation, and excluding the first 100
% cached samples from every recording. The subjects and activities
% represented in the exported point cloud are:
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
% previous stage of analysis is implemented in
% fit_agc_to_illuminance_util.py.

% The minispect illuminance is related to the average scene luminance by a
% factor of pi.

%% I don't understand how this works. I can't find this function
%[~, best_p, sensitivity, illuminance] = fitCameraAGCToIlluminance([]);
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'camera_agc_illuminance_linear_scale.mat');
load(dataFileName);
empiricalCameraScore = linear_scale_data.x_linear_scale;
empiricalMeanLuminance = linear_scale_data.y_linear_scale / pi;

% Plot the empirical measurements in log space
figure;
s = scatter(empiricalCameraScore, empiricalMeanLuminance, 25, 'MarkerEdgeColor','k');
s.Marker = '.';
s.MarkerEdgeAlpha = 0.1;
a = gca();
a.XScale = 'log';
a.YScale = 'log';
a.TickDir = 'out';
hold on; grid off; box off;

% Now superimpose the laboratory measurements of mean sphere luminance
% against the camera score
loglog(cameraScore, sphereAvgLuminance,'-*r','LineWidth',2,'MarkerSize',10); 

% Clean up, label, legend
xlabel('Log camera sensitivity score');
ylabel('Log mean environmental luminance (cd/m2)');
title('Average scene luminance vs. camera AGC sensitivity');
legend('Environmental Data', 'Integrating sphere measures');

% Save the values that relate camera score to average scene luminance
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'cameraScoreToAverageLuminance.mat');
avgSceneLuminance = sphereAvgLuminance;
readme = ['Created by defineAGCToMeanLuminance.\n'...
    'A linear interpolation between these values (in log10 space) maps AGC values to luminance.\n',...
    'cameraScore -- the product of the AGC settings (analog gain, digital gain, exposure).\n',...
    'avgSceneLuminance -- the average luminance (cd/m2) of the scene viewed by the camera.\n'];
save(saveFileName,'readme','cameraScore','avgSceneLuminance');
