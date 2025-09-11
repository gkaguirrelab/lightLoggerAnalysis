% Load the spectral sensitivity functions for the AS chip
% Save the path to CombiExperiments. We will use this as a relative
% path to find other files
lightLoggerAnalysis_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');
load([lightLoggerAnalysis_path '/data/ASM7341_spectralSensitivity.mat'])
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
title('Spectral Sensitivity of MS Sensor Channels');
xlabel('Wavelength (nm)');
ylabel('Relative Sensitivity');
% Create a Fourier basis set of the same rank
N = MakeFourierBasis(wlsMS,nChan)';
% Create a synthetic spectral source using the Fourier set; place it in the
% positive domain. Make sure it is a column vector (Nx1)
S = (rand(1,8)*N)';
S = S-min(S)+1;
% Get all the chunk files from the directory
chunk_dir = '/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/sophia_in_wild_converted_for_sam/';
chunk_files = dir(fullfile(chunk_dir, 'chunk_*.mat'));
% --- Sort chunk files numerically ---
file_names = {chunk_files.name};
% Use regular expressions to extract the number from each filename
file_numbers = cellfun(@(x) sscanf(x, 'chunk_%d.mat'), file_names);
% Sort the numbers and get the new index
[~, sorted_idx] = sort(file_numbers);
% Reorder the chunk_files struct
chunk_files = chunk_files(sorted_idx);
% Now load the world camera spectral sensitivity functions. The resulting
% matrix J has row vectors for the sensitivity of the Blue, Green, and Red
% channels, in this order.
load([lightLoggerAnalysis_path '/data/IMX219_spectralSensitivity.mat']);
wlsCam = T{:,1};
idxCam = and(wlsCam>=400,wlsCam<=750);
wlsCam = wlsCam(idxCam);
J = T{idxCam,2:4}';
% Initialize cumulative arrays for all data
tMS_all_cumulative = [];
wCam_all_cumulative = nan(3,0);
tCam_all_cumulative = [];
BGR_norm_all_cumulative = nan(0,3);

% Option to run all chunks or just a subset
run_all_chunks = false; % Set to true to process all chunks
if run_all_chunks
    num_chunks_to_process = length(chunk_files);
else
    num_chunks_to_process = 3; % Process only the first 3 chunks
end
% Loop through each chunk file
for fileIdx = 1:num_chunks_to_process
    % Load the chunk data
    chunk_path = fullfile(chunk_dir, chunk_files(fileIdx).name);
    chunk = load(chunk_path).chunk;
    
    % Get the timestamps and data for the current chunk
    tMS_chunk = chunk.M.t.AS;
    nMSReadings = size(chunk.M.v.AS, 1);
    
    tCam_chunk = chunk.W.t;
    vCam_chunk = chunk.W.v;
    nFrames_chunk = size(vCam_chunk, 1);
    
    % --------process all MS data for the current chunk-------
    wCam_chunk = nan(3, nMSReadings);
    for MSidx = 1:nMSReadings
        wMS = chunk.M.v.AS(MSidx,1:8)';
        S_hat = pinv(M)*wMS;
        wCam = J*S_hat;
        wCam = wCam / wCam(3);
        wCam_chunk(:,MSidx) = wCam;
    end
    
    % --------process all world camera data for the current chunk-------
    resCam = size(vCam_chunk, 2,3);
    BGR_norm_chunk = nan(nFrames_chunk,3);
    [rows, cols] = ndgrid(1:resCam(1), 1:resCam(2));
    R_mask = mod(rows, 2) == 0 & mod(cols, 2) == 0;
    G_mask = (mod(rows, 2) == 0 & mod(cols, 2) == 1) | ...
             (mod(rows, 2) == 1 & mod(cols, 2) == 0);
    B_mask = mod(rows, 2) == 1 & mod(cols, 2) == 1;
    WORLD_RGB_SCALARS = [1, 0.8223775, 0.95367937]; %from Zach
    for iFrame = 1:nFrames_chunk
        thisFrame = squeeze(vCam_chunk(iFrame,:,:));
        R_mean = mean(thisFrame(R_mask))*WORLD_RGB_SCALARS(1);
        G_mean = mean(thisFrame(G_mask))*WORLD_RGB_SCALARS(2);
        B_mean = mean(thisFrame(B_mask))*WORLD_RGB_SCALARS(3);
        % normalize and save
        BGR_norm_chunk(iFrame,:) = [B_mean/R_mean, G_mean/R_mean, R_mean/R_mean];
    end
    
    % Append data to cumulative arrays with the time offset
    tMS_all_cumulative = [tMS_all_cumulative, tMS_chunk];
    wCam_all_cumulative = [wCam_all_cumulative, wCam_chunk];
    tCam_all_cumulative = [tCam_all_cumulative, tCam_chunk];
    BGR_norm_all_cumulative = [BGR_norm_all_cumulative; BGR_norm_chunk];
end
% Initialize a new figure for plotting all data
figure;
hold on
title('MS-estimated RGB weights vs. avg camera RGB weights');
ylabel('Normalized Weight');
xlabel('Time (s)');
% Plot the connected MS lines using the new symbols
plot(tMS_all_cumulative, wCam_all_cumulative(1,:), 'bo');
plot(tMS_all_cumulative, wCam_all_cumulative(2,:), 'go');
plot(tMS_all_cumulative, wCam_all_cumulative(3,:), 'ro');
% Plot the camera data after the loop using the new symbols
plot(tCam_all_cumulative, BGR_norm_all_cumulative(:,1), 'b.');
plot(tCam_all_cumulative, BGR_norm_all_cumulative(:,2), 'g.');
plot(tCam_all_cumulative, BGR_norm_all_cumulative(:,3), 'r.');
hold off
ylim([0,max(max(wCam_all_cumulative, [], 'all'), max(BGR_norm_all_cumulative, [], 'all'))]);
% Create the legend for all lines
legend('MS estimated B', 'MS estimated G', 'MS estimated R', ...
       'Cam Avg. B', 'Cam Avg. G', 'Cam Avg. R', 'Location', 'best');
