%% SCRIPT TO PLOT MINISPECT DATA ACROSS MULTIPLE CHUNKS
% 
% Adjust settings as necessary.


%% SETTINGS
directory   = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/sophia_in_wild_7-22/sophia_in_wild_7-22_chunks';
files       = dir(fullfile(directory,'chunk_*.mat'));

% Preallocate cell array for storing (number, file) pairs
file_info = [];

for i = 1:numel(files)
    match = regexp(files(i).name, '^chunk_(\d+)\.mat$', 'tokens');
    if ~isempty(match)
        chunk_num = str2double(match{1}{1});
        file_info(end+1).number = chunk_num;
        file_info(end).file = files(i);
    else
        warning('Invalid file: %s', files(i).name);
    end
end

% Convert to struct array and sort by chunk number
file_info = struct2table(file_info);
file_info = sortrows(file_info, 'number');

% Recover sorted files array
files = file_info.file;

N = numel(files);

fsVid     = 120;      % world‐camera frame rate (Hz)
Twin      = 5;        % half‐window (s) → 10 s total
hop       = 2.5;      % shift between windows (s)
fMax      = 90;       % highest frequency to plot (Hz)

% Collect all AS7 values across all chunks
AS7_all = [];

% GET GLOBAL AVG ACROSS CHUNKS
for i = 1:N
    C = load(fullfile(files(i).folder, files(i).name), 'chunk');
    ch = C.chunk;
    AS7 = ch.M.v.AS(:,7);
    AS7_all = [AS7_all; AS7(:)];
end

as_min = min(AS7_all);
as_max = max(AS7_all);
as_norm = (AS7_all - as_min) / (as_max - as_min);

% normalized threshold
level_norm = graythresh(as_norm);  

% convert back
thr = level_norm*(as_max - as_min) + as_min;

%% AS BRIGHTNESS OVER TIME        
figure; hold on;
for i = 1:N
    C   = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch  = C.chunk;
    tAS = ch.M.t.AS(:);        % AS timestamps
    AS7 = ch.M.v.AS(:,7);      % AS brightness
    plot(tAS, AS7);
    set(gca, 'YScale', 'log');
    yline(thr)
end 

hold off;
xlabel('Time (s)');
ylabel('AS channel 7 (brightness)');
title('AS Brightness over Time (log scale)');
legend({files.name},'Interpreter','none','Location','best');

%% BUILD COMMON FREQUENCY GRID
% Use a nominal 10 s block from chunk 1 to get its frequency bins
tmp      = load(fullfile(files(1).folder,files(1).name),'chunk');
Vid0     = tmp.chunk.W.v;
[~, frq0] = calcTemporalSPD( Vid0(1:round(2*Twin*fsVid),:,:), fsVid, 'lineResolution', false );
f_min    = min(frq0(frq0>0));
% bounds
f_start = ceil(f_min); 
f_end   = floor(fsVid / 2) - 1; 
f_int   = (f_start : f_end)';

%% SLIDING‐WINDOW SPD SPLIT & AVERAGE ACROSS CHUNKS

allHi = nan(N, numel(f_int));
allLo = nan(N, numel(f_int));

% chunk loop
for i = 1:N
    % load
    C     = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch    = C.chunk;
    tAS   = ch.M.t.AS(:);
    AS7   = ch.M.v.AS(:,7);
    W_t   = ch.W.t(:);
    Vid   = ch.W.v;
    
    % threshold & AS→frame interpolation
    AS_if   = interp1(tAS, AS7, W_t, 'pchip', NaN);
    t0s     = (W_t(1)+Twin):hop:(W_t(end)-Twin);
    
    % PSD window loop
    spdHi = [];
    spdLo = [];
    for t0 = t0s
        idx = (W_t>=t0-Twin)&(W_t<=t0+Twin);
        if nnz(idx)<2
            continue;
        end

        as_val = mean(AS_if(idx));
        isHi = as_val > thr;

        fprintf('Chunk %d | t0 = %.2f | nnz(idx) = %d\n', i, t0, nnz(idx));

        % compute & trim PSD
        [P,f] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false);
        fprintf(' → PSD size = [%d %d], freq range = [%.2f %.2f]\n', size(P), min(f), max(f));

        keep  = f>0 & f<=fMax;
        fprintf(' → # kept freqs = %d\n', nnz(keep));

        frqLoc = f(keep);
        spdLoc = P(keep);

        % 30HZ CENSOR
        spdLoc(frqLoc==30) = NaN;

        if numel(frqLoc) < 2 || numel(spdLoc) < 2
            fprintf(' → Too few points for interp1, skipping\n');
            continue;
        end

        % interpolate
        Pi = interp1(frqLoc, spdLoc, f_int, 'pchip')';    
        % GUARD AGAINST EMPTY/SHORT VECTORS
        if numel(Pi) ~= numel(f_int)
            fprintf(' → Pi size mismatch (%d vs %d), skipping\n', numel(Pi), numel(f_int));
            continue
        end
        
        % bucket
        if isHi
            spdHi(end+1,:) = Pi;
        else
            spdLo(end+1,:) = Pi;
        end
    end
    
    % store this chunk's mean curves
    if ~isempty(spdHi)
        allHi(i,:) = mean(spdHi, 1);
    else
        allHi(i,:) = nan(1, numel(f_int));  % fill with NaNs for clean averaging
    end

    if ~isempty(spdLo)
        allLo(i,:) = mean(spdLo, 1);
    else
        allLo(i,:) = nan(1, numel(f_int));  % fill with NaNs for clean averaging
    end

    fprintf('[Chunk %d] spdHi windows: %d | spdLo windows: %d\n', ...
    i, size(spdHi, 1), size(spdLo, 1));
end

% global average across chunks
globalHi = nanmean(allHi,1);
globalLo = nanmean(allLo,1);

%% PLOT GLOBAL SPD (high vs low)
figure;
plotSPD(globalHi, f_int);
hold on;
plotSPD(globalLo, f_int, 'r');
hold off;

xlabel('Frequency (Hz)');
ylabel('Spectral power density (contrast^2/Hz)');
title('All Chunks: Avg Windowed SPD split by High vs Low AS (1–90 Hz)');
xlim([f_int(1) f_int(end)]);
nticks = 10;
ticks = round(logspace(log10(0.2),log10(59),nticks), 2);
xticks(ticks)
xticklabels(ticks)
legend('High AS','Low AS','Location','best');

%% SLOPE MAPS

% preallocate
vHiAll = [];
vLoAll = [];

% chunk loop
for i = 1:N
    % load
    C   = load(fullfile(files(i).folder, files(i).name), 'chunk');
    ch  = C.chunk;
    Vid = ch.W.v;
    tAS = ch.M.t.AS(:);
    AS7 = ch.M.v.AS(:,7);
    Wt  = ch.W.t(:);

    % interpolate AS onto each video frame time
    AS_if = interp1(tAS, AS7, Wt, 'pchip', NaN);

    % split into high vs low
    hiIdx = AS_if >  thr;
    loIdx = AS_if <= thr;

    % accumulate frames
    vHiAll = cat(1, vHiAll, ch.W.v(hiIdx,:,:));
    vLoAll = cat(1, vLoAll, ch.W.v(loIdx,:,:));
end

% compute & plot high AS
[Xh,Yh,Zh, slopeHigh, interceptHigh, frq, hFigHighSlope, hFigHighIntercept] = mapSlopeSPD(vHiAll, fsVid, [40,40], 20, fisheyeIntrinsics);

% compute & plot low AS
[Xl,Yl,Zl, slopeLow,  interceptLow, ~, hFigLowSlope, hFigLowIntercept] = mapSlopeSPD(vLoAll, fsVid, [40,40], 20, fisheyeIntrinsics);

% find common limits
vminSlope = min( [slopeHigh(:); slopeLow(:)] );
vmaxSlope = max( [slopeHigh(:); slopeLow(:)] );

vminInt = min( [interceptHigh(:); interceptLow(:)] );
vmaxInt = max( [interceptHigh(:); interceptLow(:)] );

% re-apply caxis and add titles
figure(hFigHighSlope);    caxis([vminSlope vmaxSlope]); title('1/f Slope Map - High AS');
figure(hFigHighIntercept);caxis([vminInt vmaxInt]);     title('1/f Intercept Map - High AS');
figure(hFigLowSlope);     caxis([vminSlope vmaxSlope]); title('1/f Slope Map - Low AS');
figure(hFigLowIntercept); caxis([vminInt vmaxInt]);     title('1/f Intercept Map - Low AS');

%% LOAD & PROJECT CALIBRATION DOTS
% Load video
calibFile = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/SophiaGazeCalib2/W_converted.avi';
vr = VideoReader(calibFile);

% Video parameters
fps = vr.FrameRate;
startSec = 24;    % Start at 2:42 from the end (which is 24s from the start)
endSec   = 177;   % End at 0:09 from the end (which is 177s from start)

% Move to start time
vr.CurrentTime = startSec;

% Prepare max brightness mask
maxBright = [];
while hasFrame(vr) && vr.CurrentTime <= endSec
    frm = readFrame(vr);
    g   = rgb2gray(frm);
    
    % Only keep pixels brighter than a threshold
    % Normalize and threshold
    g = double(g);
    g(g < 220) = 0; % Adjust this threshold as needed (try 220, 240, etc.)

    % Combine max
    if isempty(maxBright)
        maxBright = g;
    else
        maxBright = max(maxBright, g);
    end
end

% Convert to uint8 for display
maxBright = uint8(maxBright);

% Resize to match calibration data
[nR, nC] = size(Xh); % Use size of one of your calibration surfaces
maxBright = imresize(maxBright, [nR, nC]);

% Overlay on calibration figures
figs = [hFigHighSlope, hFigHighIntercept, hFigLowSlope, hFigLowIntercept];
Xs = {Xh, Xh, Xl, Xl};
Ys = {Yh, Yh, Yl, Yl};
Zs = {Zh, Zh, Zl, Zl};

for k = 1:4
    figure(figs(k));
    hold on;
    hTex = surf(Xs{k}, Ys{k}, Zs{k}, ...
        'CData', repmat(maxBright, 1, 1, 3), ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none');
    alpha(hTex, 0.6);
    hold off;
end
