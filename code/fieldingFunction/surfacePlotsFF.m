%% SurfacePlotsFF.m
% The code below creates 3D surface plots for each original image taken
% inside the integrating sphere, as well as for the average of these
% images. 

% Load .mat files (FIRST PAIR OF IMAGES)
% A = load('reg_sophia_world_frame.mat');
% B = load('updown_sophia_world_frame.mat');

% Load .mat files (SECOND PAIR OF IMAGES)
A = load('regularRotation.mat');
B = load('invertedRotation.mat');

Ireg = double(A.image);
Iup  = double(B.image);

% 3D surface plot for each original image
figure
surf(Ireg, 'EdgeColor', 'none')
title('Regular image surface')
xlabel('X'); ylabel('Y'); zlabel('Pixel intensity')
view(3); colorbar

figure
surf(Iup, 'EdgeColor', 'none')
title('Updown image surface')
xlabel('X'); ylabel('Y'); zlabel('Pixel intensity')
view(3); colorbar

% Rotate one image 180 degrees
Iup_rot = rot90(Iup, 2);

% Average rotated image with the other image
Iavg = (Ireg + Iup_rot) / 2;

% 3D surface plot of average image
figure
surf(Iavg, 'EdgeColor', 'none')
title('Average of regular image and 180° rotated updown image')
xlabel('X'); ylabel('Y'); zlabel('Pixel intensity')
view(3); colorbar

%% SurfacePlotsFF_NaNAwareAverage.m
% The code prompts the user to draw a mask of the obstructing baffle object 
% in each provided image. A NaN mask corresponding to the identified location
% is applied to both images before creating an average, like done so above.

clear; close all; clc

% A = load('reg_sophia_world_frame.mat');
% B = load('updown_sophia_world_frame.mat');

A = load('regularRotation.mat');
B = load('invertedRotation.mat');

Ireg = double(A.image);
Iup  = double(B.image);

% Plot original images
figure
surf(Ireg, 'EdgeColor', 'none')
title('Regular image surface')
view(3); colorbar

figure
surf(Iup, 'EdgeColor', 'none')
title('Updown image surface')
view(3); colorbar

% Draw baffle mask on regular image
figure
imshow(Ireg, [])
title('Draw baffle in regular image')
roi1 = drawfreehand;
mask_reg = createMask(roi1);

% Draw baffle mask on updown image BEFORE rotation
figure
imshow(Iup, [])
title('Draw baffle in updown image before rotation')
roi2 = drawfreehand;
mask_up = createMask(roi2);

% NaN out baffles
Ireg_nan = Ireg;
Iup_nan = Iup;

Ireg_nan(mask_reg) = NaN;
Iup_nan(mask_up) = NaN;

% Rotate updown image after masking
Iup_nan_rot = rot90(Iup_nan, 2);

% NaN-aware average across images
stack = cat(3, Ireg_nan, Iup_nan_rot);
Iavg_nan = mean(stack, 3, 'omitnan');

% Plot NaN-aware average
figure
imagesc(Iavg_nan)
axis image
colorbar
title('NaN-aware average')

figure
surf(Iavg_nan, 'EdgeColor', 'none')
title('NaN-aware average surface')
view(3); colorbar

% COMPARE IMAGES
figure

subplot(1,2,1)
imagesc(Ireg)
axis image
colorbar
title('Regular')

subplot(1,2,2)
imagesc(Iup)
axis image
colorbar
title('UpDown')



% %% SurfacePlotsBaffleMaskFIRST
% 
% % Load both .mat files
% A = load('reg_sophia_world_frame.mat');
% B = load('updown_sophia_world_frame.mat');
% 
% Ireg = double(A.image);
% Iup  = double(B.image);
% 
% % Draw baffle mask on regular image
% figure; imshow(Ireg, [])
% title('Draw baffle in regular image')
% roi1 = drawfreehand;
% mask_reg = createMask(roi1);
% 
% % Draw baffle mask on updown image BEFORE rotating
% figure; imshow(Iup, [])
% title('Draw baffle in updown image BEFORE rotation')
% roi2 = drawfreehand;
% mask_up = createMask(roi2);
% 
% % NaN out baffles first
% Ireg_nan = Ireg;
% Iup_nan  = Iup;
% 
% Ireg_nan(mask_reg) = NaN;
% Iup_nan(mask_up)   = NaN;
% 
% % Now rotate the updown image and its NaN mask/data
% Iup_nan_rot = rot90(Iup_nan, 2);
% 
% % NaN-aware average
% stack = cat(3, Ireg_nan, Iup_nan_rot);
% Iavg_nan = mean(stack, 3, 'omitnan');
% 
% % Plot
% figure
% imagesc(Iavg_nan)
% axis image
% colorbar
% title('NaN-aware average after masking baffles first')
% 
% figure
% surf(Iavg_nan, 'EdgeColor', 'none')
% view(3)
% colorbar
% title('NaN-aware average surface')
% 
% %% SurfacePlotsBaffleMaskLAST
% 
% % Load both .mat files
% A = load('reg_sophia_world_frame.mat');
% B = load('updown_sophia_world_frame.mat');
% 
% Ireg = double(A.image);
% Iup  = double(B.image);
% 
% % Manually mask baffle/sliver in the averaged image
% figure
% imshow(Iavg, [])
% title('Draw around remaining baffle/sliver in averaged image')
% roi = drawfreehand;
% mask_avg = createMask(roi);
% 
% Iavg_cut = Iavg;
% Iavg_cut(mask_avg) = NaN;
% 
% % Plot top-down
% figure
% imagesc(Iavg_cut)
% axis image
% colorbar
% title('Averaged image with baffle removed last')
% 
% % Plot surface
% figure
% surf(Iavg_cut, 'EdgeColor', 'none')
% view(3)
% colorbar
% title('Surface after removing baffle from averaged image')
% 
