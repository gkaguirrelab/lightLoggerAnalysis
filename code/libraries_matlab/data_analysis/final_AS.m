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

%% AS BRIGHTNESS OVER TIME (red = above thr, blue = below thr)
figure; hold on;

% Legend anchors (so legend doesn’t list every chunk)
hHigh = plot(nan,nan,'r-','LineWidth',1.2,'DisplayName','High (>= thr)');
hLow  = plot(nan,nan,'b-','LineWidth',1.2,'DisplayName','Low (< thr)');

for i = 1:N
    C   = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch  = C.chunk;
    tAS = ch.M.t.AS(:);        % time (s)
    AS7 = ch.M.v.AS(:,7);      % brightness

    % Split into above/below threshold (use NaNs to break lines)
    yHigh = AS7;              yHigh(AS7 <= thr) = NaN;
    yLow  = AS7;              yLow(AS7  > thr) = NaN;

    % Plot; hide handles so legend doesn’t get cluttered
    plot(tAS, yHigh, 'r-', 'LineWidth', 1.0, 'HandleVisibility','off');
    plot(tAS, yLow,  'b-', 'LineWidth', 1.0, 'HandleVisibility','off');
end

% Axes, threshold line, labels, legend
set(gca, 'YScale','log', 'FontSize',14, 'TickLength',[0.02 0.02]);
yline(thr, 'k--', 'HandleVisibility','off');

xlabel('Time (s)', 'FontSize',14);
ylabel('AS channel 7 (brightness)', 'FontSize',14);
title('AS Brightness over Time (log scale)', 'FontSize',16, 'FontWeight','normal');
legend([hHigh hLow], 'Location','best');  % only two entries
set(gcf, 'color','white');
hold off;


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
allCen = nan(N, numel(f_int));
allPer = nan(N, numel(f_int));

% mask specs
centerRows = 121:360;
centerCols = 161:480;

tileSize = 2;
maxTilesPerRegion = 4000;

% chunk loop
for i = 1:N
    % load
    C     = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch    = C.chunk;
    tAS   = ch.M.t.AS(:);
    AS7   = ch.M.v.AS(:,7);
    W_t   = ch.W.t(:);
    Vid   = ch.W.v;

    [nFrames, nRows, nCols] = size(Vid);
    fps = fsVid;

    % create central and peripheral mask
    centerMask = false(nRows, nCols);
    centerMask(centerRows, centerCols) = true;
    peripheryMask = ~centerMask;

    % ---------- BUILD TILE SUBMASKS FULLY INSIDE EACH REGION ----------
    subMasks_center  = {};
    subMasks_periph  = {};
    for r = 1:tileSize:nRows
        rr = r:min(r+tileSize-1, nRows);
        for c = 1:tileSize:nCols
            cc = c:min(c+tileSize-1, nCols);
            if all(centerMask(rr,cc), 'all')
                M = false(nRows,nCols); M(rr,cc) = true;
                subMasks_center{end+1} = M;
            elseif all(peripheryMask(rr,cc), 'all')
                M = false(nRows,nCols); M(rr,cc) = true;
                subMasks_periph{end+1} = M; 
            end
        end
    end
    
    % threshold & AS→frame interpolation
    AS_if   = interp1(tAS, AS7, W_t, 'pchip', NaN);
    t0s     = (W_t(1)+Twin):hop:(W_t(end)-Twin);
    
    % accumulators
    spdHi = [];
    spdLo = [];
    spdCen = [];
    spdPer = [];

    % PSD window loop
    for t0 = t0s
        idx = (W_t>=t0-Twin)&(W_t<=t0+Twin);
        if nnz(idx)<2
            continue;
        end

        as_val = mean(AS_if(idx));
        isHi = as_val > thr;
        fprintf('Chunk %d | t0 = %.2f | nnz(idx) = %d\n', i, t0, nnz(idx));

        % % -------- Whole-frame SPD (existing Hi/Lo path) --------
        [P,f] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false);
        fprintf(' → PSD size = [%d %d], freq range = [%.2f %.2f]\n', size(P), min(f), max(f));

        keep  = f>0 & f<=fMax;
        fprintf(' → # kept freqs = %d\n', nnz(keep));
        frqLoc = f(keep);
        spdLoc = P(keep);

        % 30HZ CENSOR
        spdLoc(frqLoc==30) = NaN;

        if numel(frqLoc) >= 2 && numel(spdLoc) >= 2
            Pi = interp1(frqLoc, spdLoc, f_int, 'pchip')';
            if numel(Pi) == numel(f_int)
                if isHi
                    spdHi(end+1,:) = Pi;
                else
                    spdLo(end+1,:) = Pi;
                end
            end
        end

    % ======== PSD→Avg using regionMatrix with MANY SMALL MASKS ========
        % Center tiles
        cen_tile_curves = [];
        for k = 1:numel(subMasks_center)
            [Pc, fc] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false, ...
                                       'regionMatrix', double(subMasks_center{k}));
            keepC = fc>0 & fc<=fMax;
            if nnz(keepC) < 2, continue; end
            PiC = interp1(fc(keepC), Pc(keepC), f_int, 'pchip')';
            PiC(f_int==30) = NaN;
            cen_tile_curves(end+1,:) = PiC;
        end
        if ~isempty(cen_tile_curves)
            spdCen(end+1,:) = mean(cen_tile_curves, 1, 'omitnan');
        end
        % Periphery tiles
        per_tile_curves = [];
        for k = 1:numel(subMasks_periph)
            [Pp, fp] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false, ...
                                       'regionMatrix', double(subMasks_periph{k}));
            keepP = fp>0 & fp<=fMax;
            if nnz(keepP) < 2, continue; end
            PiP = interp1(fp(keepP), Pp(keepP), f_int, 'pchip')';
            PiP(f_int==30) = NaN;
            per_tile_curves(end+1,:) = PiP;
        end
        if ~isempty(per_tile_curves)
            spdPer(end+1,:) = mean(per_tile_curves, 1, 'omitnan');
        end
    end
    
    % store this chunk's mean curves
    allHi(i,:) = mean(spdHi, 1, 'omitnan');
    allLo(i,:) = mean(spdLo, 1, 'omitnan');
    
    % again for center and periphery
    allCen(i,:) = mean(spdCen, 1, 'omitnan');
    allPer(i,:) = mean(spdPer, 1, 'omitnan');

    fprintf('[Chunk %d] spdHi windows: %d | spdLo windows: %d\n', ...
    i, size(spdHi, 1), size(spdLo, 1), size(spdCen, 1), size(spdPer, 1));
end

% global average across chunks
globalHi = nanmean(allHi,1);
globalLo = nanmean(allLo,1);
globalCen = nanmean(allCen, 1);
globalPer = nanmean(allPer, 1);


%% PLOT GLOBAL SPD (high vs low)
figure;
plotSPD(globalHi, f_int, 'r');
hold on;
plotSPD(globalLo, f_int, 'b');
hold off;

xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Spectral power density (contrast^2/Hz)', 'FontSize',14);
title('All Chunks: Avg Windowed SPD split by High vs Low AS (1–90 Hz)', ...
      'FontSize',16, 'FontWeight','normal');

xlim([f_int(1) f_int(end)]);
nticks = 10;
ticks = round(logspace(log10(0.2),log10(59),nticks), 2);
xticks(ticks)
xticklabels(ticks)
legend('High AS','Low AS','Location','best');

legend({'High AS','Low AS'}, 'Location','best');
set(gca, 'FontSize',14, 'TickLength',[0.02 0.02]);
set(gcf, 'Color','white');


%% PLOT GLOBAL SPD (center vs surround)
figure;

% Center (red), Surround (blue)
plotSPD(globalCen, f_int, 'm');
hold on;
plotSPD(globalPer, f_int, 'g');
hold off;

xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Spectral power density (contrast^2/Hz)', 'FontSize',14);
title('All Chunks: Avg Windowed SPD — Center vs Surround (1–90 Hz)', ...
      'FontSize',16, 'FontWeight','normal');

xlim([f_int(1) f_int(end)]);
nticks = 10;
ticks = round(logspace(log10(0.2), log10(59), nticks), 2);
xticks(ticks);
xticklabels(ticks);

legend({'Center','Surround'}, 'Location','best');
set(gca, 'FontSize',14, 'TickLength',[0.02 0.02]);
set(gcf, 'Color','white');


%% LOAD CALIBRATION DOTS & CALC AFFINE TRANSFORMATION

calibFile = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/SophiaGazeCalib2/W_mjpeg.avi';
vr = VideoReader(calibFile);

% video parameters
fps = vr.FrameRate;
startSec = 82;
endSec   = 276;

% move to start time
vr.CurrentTime = startSec;

% prepare max brightness mask
maxBright = [];
while hasFrame(vr) && vr.CurrentTime <= endSec
    frm = readFrame(vr);
    g   = rgb2gray(frm);
    
    % only keep pixels brighter than a threshold
    % normalize and threshold
    g = double(g);
    g(g < 200) = 0; % Adjust this threshold as needed (try 220, 240, etc.)

    % combine max
    if isempty(maxBright)
        maxBright = g;
    else
        maxBright = max(maxBright, g);
    end
end

% convert to uint8 for display
maxBright = uint8(maxBright);

% extract dot centroids
BW = imbinarize(maxBright, graythresh(maxBright));  
stats = regionprops(BW, 'Centroid');
imgPts = vertcat(stats.Centroid); % k×2 [x y] in pixel coords
x = imgPts(:,1);
y = imgPts(:,2);
r = sqrt(x.^2 + y.^2);

% world point coords
worldPts = [ ...
        -20, -20;   -20, 0;   -20, 20;   -15, -15;  -15, 15; ...
        -10, -10;   -10, 0;   -10, 10;   -5, -5;    -5, 5;   ...
         0, -20;     0, -10;   0, 0;      0, 10;     0, 20;  ...
         5, -5;      5, 5;     10, -10;   10, 0;     10, 10; ...
         15, -15;    15, 15;   20, -20;   20, 0;     20, 20];

% Compute affine transform
tform = fitgeotrans(imgPts, worldPts, 'affine');
affineMat = tform.T;


%% SLOPE AND INTERCEPT MAPS

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

% Use the frame size from the video (or from vHiAll)
[~, nRows, nCols] = size(vHiAll);

% Camera-space angles (radians)
[theta_cam, phi_cam, r] = anglesFromIntrinsics(nRows, nCols, fisheyeIntrinsics);

% Warp to participant space using affine (apply in DEGREES to [az, el])
az_deg = rad2deg(phi_cam);
el_deg = rad2deg(theta_cam);

azelW      = [az_deg(:), el_deg(:), ones(numel(az_deg),1)] * affineMat.';  % N×2 deg
azW_deg    = reshape(azelW(:,1), size(az_deg));
elW_deg    = reshape(azelW(:,2), size(el_deg));
phi_w      = deg2rad(azW_deg);
theta_w    = deg2rad(elW_deg);

% High AS maps on participant-space sphere
[hFigHighSlope, hFigHighIntercept, frq] = ... 
    mapSlopeIntSPD(vHiAll, fsVid, [40,40], 20, theta_w, phi_w, r);

% Low AS maps on participant-space sphere
[hFigLowSlope,  hFigLowIntercept, ~]  = ...
    mapSlopeIntSPD(vLoAll, fsVid, [40,40], 20, theta_w, phi_w, r);

% find common limits
vminSlope = min( [slopeHigh(:); slopeLow(:)] );
vmaxSlope = max( [slopeHigh(:); slopeLow(:)] );
vminInt = min( [interceptHigh(:); interceptLow(:)] );
vmaxInt = max( [interceptHigh(:); interceptLow(:)] );

% re-apply caxis and add titles
figure(hFigHighSlope);    caxis([vminSlope vmaxSlope]); title('1/f Slope Map - High AS', FontSize=16);
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white')
figure(hFigHighIntercept);caxis([vminInt vmaxInt]);     title('1/f Intercept Map - High AS');
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white')
figure(hFigLowSlope);     caxis([vminSlope vmaxSlope]); title('1/f Slope Map - Low AS');
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white')
figure(hFigLowIntercept); caxis([vminInt vmaxInt]);     title('1/f Intercept Map - Low AS');
set(gca, 'FontSize', 14);
set(gcf, 'color', 'white')


%% PROJECT DOTS

img_aug   = [imgPts, ones(size(imgPts,1),1)];
azel_dots = img_aug * affineMat.';          % deg [az el]
azd = deg2rad(azel_dots(:,1));
eld = deg2rad(azel_dots(:,2));
Xd  = sin(eld).*cos(azd);
Yd  = sin(eld).*sin(azd);
Zd  = cos(eld);

for f = [hFigHighSlope, hFigHighIntercept, hFigLowSlope, hFigLowIntercept]
    figure(f); hold on;
    scatter3(Xd, Yd, Zd, 60, 'k', 'filled');
    hold off;
end
