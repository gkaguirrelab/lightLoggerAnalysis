%% SETTINGS
directory = '/Users/sophiamirabal/Desktop/sophia_in_wild_7-22_chunks';
files     = dir(fullfile(directory,'chunk_*.mat'));
N         = numel(files);

fsVid     = 120;      % world‐camera frame rate (Hz)
Twin      = 5;        % half‐window (s) → 10 s total
hop       = 2.5;      % shift between windows (s)
fMax      = 90;       % highest frequency to plot (Hz)

%% 1) AS BRIGHTNESS OVER TIME (all chunks)
figure; hold on;
for i = 1:N
    C   = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch  = C.chunk;
    tAS = ch.M.t(:);           % AS timestamps
    AS7 = ch.M.v.AS(:,7);      % AS brightness
    semilogy(tAS, AS7, '-o','LineWidth',1.2);
end
hold off;
xlabel('Time (s)');
ylabel('AS channel 7 (brightness)');
title('AS Brightness over Time (log scale)');
legend({files.name},'Interpreter','none','Location','best');

%% 2) BUILD COMMON FREQUENCY GRID
% Use a nominal 10 s block from chunk 1 to get its frequency bins
tmp      = load(fullfile(files(1).folder,files(1).name),'chunk');
Vid0     = tmp.chunk.W.v;
[~, frq0] = calcTemporalSPD( Vid0(1:round(2*Twin*fsVid),:,:), fsVid );
f_min    = min(frq0(frq0>0));
f_int    = (ceil(f_min):1:fMax)';   % integer‐Hz grid

%% 3) SLIDING‐WINDOW SPD SPLIT & AVERAGE ACROSS CHUNKS
allHi = nan(N, numel(f_int));
allLo = nan(N, numel(f_int));

for i = 1:N
    % load
    C     = load(fullfile(files(i).folder,files(i).name),'chunk');
    ch    = C.chunk;
    tAS   = ch.M.t(:);
    AS7   = ch.M.v.AS(:,7);
    W_t   = ch.W.t(:);
    Vid   = ch.W.v;
    
    % threshold & AS→frame interpolation
    thr     = mean(AS7);
    AS_if   = interp1(tAS, AS7, W_t, 'pchip', NaN);
    t0s     = (W_t(1)+Twin):hop:(W_t(end)-Twin);
    
    % collect each window's PSD
    spdHi = [];
    spdLo = [];
    for t0 = t0s
        idx = (W_t>=t0-Twin)&(W_t<=t0+Twin);
        if nnz(idx)<2, continue, end
        
        % require window entirely high or low
        winAS = AS_if(idx);
        isHi  = interp1(tAS,AS7,t0,'pchip') > thr;
        if any(isnan(winAS)) ...
           || any(winAS>thr  & ~isHi) ...
           || any(winAS<=thr &  isHi)
            continue
        end
        
        % compute & trim PSD
        [P,f] = calcTemporalSPD(Vid(idx,:,:), fsVid);
        keep  = f>0 & f<=fMax;
        Pi    = interp1(f(keep), P(keep), f_int, 'pchip')';
        
        % bucket
        if isHi
            spdHi(end+1,:) = Pi;
        else
            spdLo(end+1,:) = Pi;
        end
    end
    
    % store this chunk's mean curves
    allHi(i,:) = nanmean(spdHi,1);
    allLo(i,:) = nanmean(spdLo,1);
end

% global average across chunks
globalHi = nanmean(allHi,1);
globalLo = nanmean(allLo,1);

%% 4) PLOT GLOBAL SPD (high vs low)
figure;
plotSPD(globalHi, f_int, 'LineWidth',2);
hold on;
plotSPD(globalLo, f_int, '--','LineWidth',2,'Color',[1 0 0]);
hold off;

xlabel('Frequency (Hz)');
ylabel('Spectral power density (contrast^2/Hz)');
title('All Chunks: Avg Windowed SPD split by High vs Low AS (1–90 Hz)');
xlim([f_int(1) fMax]);
legend('High AS','Low AS','Location','best');
