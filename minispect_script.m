%% SETUP

directory = '/Users/sophiamirabal/Desktop/sophia_in_wild_7-22_chunks';
filename  = 'chunk_9.mat';
filepath  = fullfile(directory, filename);

C = load(filepath,'chunk');
chunk = C.chunk;

% AS brightness trace
AS7 = chunk.M.v.AS(:,7);
tAS = chunk.M.t.AS; 

% define AS average               
thr   = mean(AS7);

% world video
Vid = chunk.W.v; 
W_t   = chunk.W.t(:); 

%% AS BRIGHTNESS ON LOG-Y AXIS

figure;
semilogy(t, AS7);
yline(thr)
xlabel('Time (s)')
ylabel('AS channel 7 (brightness)')
title('AS Brightness over Time (log scale)')

%% SLIDING-WINDOW SPLIT-SPD

% define sampling window
win = 5;

% build a grid of window centers every 2.5 sec
t0s = (W_t(1)+Twin) : 2.5 : (W_t(end)-Twin);

% determine midpoint of AS brightness by interpolating
AS_interp = @(t) interp1(tAS, AS7, t, 'pchip', NaN);

% precompute to get full freq‐axis
[~,frq] = calcTemporalSPD( Vid(1:round(2*Twin*fsVid),:,:), fsVid );
f_min   = frq(find(frq>0,1));    % lowest nonzero bin
f_max   = min(90, frq(end));     % cap at 90 Hz

% integer‐Hz grid from f_min to 90 Hz
f_int   = (ceil(f_min):floor(f_max))';

% prepare dynamic storage
spdHiAll = [];  
spdLoAll = [];

for t0 = t0s
  % extract window around t0
  winIdx = (W_t >= t0-Twin) & (W_t <= t0+Twin);
  if nnz(winIdx)<2, continue; end

  [spdWin, frqWin] = calcTemporalSPD(Vid(winIdx,:,:), fsVid);
  keep = frqWin>=f_min & frqWin<=f_max;
  spdLoc = spdWin(keep);
  frqLoc = frqWin(keep);

  % resample onto integer‐Hz grid
  spdInt = interp1(frqLoc, spdLoc, f_int, 'pchip')';

  % decide high vs low by interpolated AS7 at t0
  a = AS_interp(t0);
  if isnan(a), continue; end
  if a > thr
    spdHiAll(end+1,:) = spdInt;
  else
    spdLoAll(end+1,:) = spdInt;
  end
end

% average across all windows
meanHi = nanmean(spdHiAll, 1);
meanLo = nanmean(spdLoAll, 1);


%% PLOT AVERAGED SPDs (1–90 Hz) WITH plotSPD

figure;
plotSPD(meanHi, f_int);
hold on;
plotSPD(meanLo, f_int, 'r');
hold off;

legend('AS > mean','AS ≤ mean','Location','best');
xlabel('Frequency (Hz)');
ylabel('Spectral power density (contrast^2/Hz)');
title('Chunk 9: Windowed SPD split by high/low AS - Shift=1');
xlim([1 90]);