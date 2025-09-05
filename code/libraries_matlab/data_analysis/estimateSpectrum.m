
% Load the spectral sensitivity functions for the AS chip
% Save the path to CombiExperiments. We will use this as a relative
% path to find other files
combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');
load([combiExperiments_path '/data/ASM7341_spectralSensitivity.mat'])

% Define the number of sensor channels
nChan = 8;

% Create a sensor matrix (M) that has just the visible wavelength channels
% and wavelength support between 400 and 750 nm. At some point we will need
% to wrestle with correcting the AS sensor counts for IR intrusion. There
% are some transpose operations so that the source spectrum ends up as a
% column vector.
wlsMS = T{:,1};
idxMS = and(wlsMS>=400,wlsMS<=750);
wlsMS = wlsMS(idxMS);
M = T{idxMS,2:2+nChan-1}';

% Show the sensor channels
figure
plot(wlsMS,M');

% Create a Fourier basis set of the same rank
N = MakeFourierBasis(wlsMS,nChan)';

% Create a synthetic spectral source using the Fourier set; place it in the
% positive domain. Make sure it is a column vector (Nx1)
S = (rand(1,8)*N)';
S = S-min(S)+1;

% Load a chunk of data
chunk = load('/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/sophia_in_wild_converted_for_sam/chunk_0.mat').chunk;
wMS = chunk.M.v.AS(end,1:8)';
tMS = chunk.M.t.AS(end);

% Now obtain the spectral source as estimated by the sensor.
S_hat = pinv(M)*wMS;

% Plot the source S, and the estimated S
figure
plot(wlsMS,S_hat);
xlabel('wavelengths [nm]');
ylabel('power');
legend({'S estimated'});

% Now load the world camera spectral sensitivity functions. The resulting
% matrix J has row vectors for the sensitivity of the Blue, Green, and Red
% channels, in this order.
load([combiExperiments_path '/data/IMX219_spectralSensitivity.mat']);
wlsCam = T{:,1};
idxCam = and(wlsCam>=400,wlsCam<=750);
wlsCam = wlsCam(idxCam);
J = T{idxCam,2:4}';

% Obtain the camera weights for the S_hat spectrum
wCam = J*S_hat;

% Because we only care about the relative weights on the camera sensor
% values, scale by the max. The resulting values describe the relative,
% average sensor count we should get for the B,G,R camera channels.
wCam = wCam / max(wCam);

% now load the world camera data from a frame at the same time as the MS
% readings
[tCamDiff, tCamidx] = min(abs(tMS - chunk.W.t));
tCam = chunk.W.t(tCamidx);
vCam = chunk.W.v(tCamidx,:,:);
vCam = squeeze(vCam);

% make sure this is a decent frame (i.e., not washed out or too dark) 
minCam = min(vCam, [], "all");
maxCam = max(vCam, [], "all");

% what are we looking at?
figure;
imshow(vCam, [0, 255]);
resCam = size(vCam);

% retreive RGB values from the world camera
[rows, cols] = ndgrid(1:resCam(1), 1:resCam(2));

R_mask = mod(rows, 2) == 0 & mod(cols, 2) == 0;
G_mask = (mod(rows, 2) == 0 & mod(cols, 2) == 1) | ...
         (mod(rows, 2) == 1 & mod(cols, 2) == 0);
B_mask = mod(rows, 2) == 1 & mod(cols, 2) == 1;

WORLD_RGB_MASK(R_mask) = 1;  % R=1
WORLD_RGB_MASK(G_mask) = 2;  % G=2
WORLD_RGB_MASK(B_mask) = 3;  % B=3

% calculate mean pixel intensity for RGB and multiply it by the channel
% weights
% RGB scalar multipliers
WORLD_RGB_SCALARS = [1, 0.8223775, 0.95367937]; %from Zach
%WORLD_RGB_SCALARS = [1.0000, 0.6190, 0.8095]; % from Pagnutti et al., 2017
R_mean = mean(vCam(R_mask))*WORLD_RGB_SCALARS(1);
G_mean = mean(vCam(G_mask))*WORLD_RGB_SCALARS(2);
B_mean = mean(vCam(B_mask))*WORLD_RGB_SCALARS(3);

maxRGB = max([R_mean, G_mean, B_mean]);
BGR_norm = [B_mean/maxRGB; G_mean/maxRGB; R_mean/maxRGB];


% Create empty 3D RGB image
rgbImage = zeros(resCam(1), resCam(2), 3);

% Assign pixels to their corresponding color channels
rgbImage(:,:,1) = double(vCam).*R_mask.*WORLD_RGB_SCALARS(1); % Red channel
rgbImage(:,:,2) = double(vCam).*G_mask.*WORLD_RGB_SCALARS(2); % Green channel
rgbImage(:,:,3) = double(vCam).*B_mask.*WORLD_RGB_SCALARS(3); % Blue channel

% Display the image
figure;
imshow(rgbImage / 255); % Normalize to [0,1] for imshow
title('Raw Bayer Pattern Visualization');



