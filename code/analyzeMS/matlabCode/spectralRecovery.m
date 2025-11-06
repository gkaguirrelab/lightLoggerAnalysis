
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

% Obtain the weights on the sensor for this spectral source
wMS = M*S;

% Now obtain the spectral source as estimated by the sensor.
S_hat = pinv(M)*wMS;

% Plot the source S, and the estimated S
figure
plot(wlsMS,S);
hold on
plot(wlsMS,S_hat);
xlabel('wavelengths [nm]');
ylabel('power');
legend({'S source','S estimated'});

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

