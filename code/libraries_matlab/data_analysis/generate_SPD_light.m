%% FUNCTION TO OBTAIN SPD AND LIGHT LEVEL INFO FROM RECORDING ACROSS MULTIPLE CHUNKS
% 
% Adjust settings as necessary.
function generate_SPD_light(directory, visualize_results)

    arguments
        directory {mustBeText} = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/sophia_in_wild_7-22/sophia_in_wild_7-22_chunks';
        analyses_to_perform = [false, true, false]; 
        visualize_results = [false, true, false];  
    end 

    % Define some constants about our analysis
    fsVid     = 120;      % world‐camera frame rate (Hz)
    Twin      = 5;        % half‐window (s) → 10 s total
    hop       = 2.5;      % shift between windows (s)
    fMax      = 90;       % highest frequency to plot (Hz)
    centerRows = 121:360; 
    centerCols = 161:480; 

    % Sort the chunks in the directory numerically 
    files = sort_chunk_filenames(directory); 

    % Ensure if we want to run slope intercept maps, 
    % we also have done MS calculation 
    if(analyses_to_perform(end) && ~analyses_to_perform(1))
        error("Must perform MS analysis to do slope map analysis");
    end 

    % Analyze the MS-AS data over the course of the video
    if(analyses_to_perform(1)) 
        [yHigh, yLow, thr] = obtain_ms_high_low(files, N_chunks);
        if(visualize_results(1))
            plot_ms_high_low(yHigh, yLow);
        end     
    end 
    
    % Obtain the SPD data over all the chunks 
    if(analyses_to_perform(2))
        [globalHi, globalLo, globalCen, globalPer, f_int] = obtain_SPD_data(files, N_chunks); 
        if(visualize_results(2))
            plot_SPD_data(globalHi, globalLo, globalCen, globalPer, f_int); 
        end 
    end 
    

    % Iterate over the chunks 
    for ii = 1:N_chunks
        % Load in the given chunk 
        C   = load(fullfile(files(ii).folder, files(ii).name), 'chunk');
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

        % Splice out the desired hi, lo frames 
        vHi = ch.W.v(hiIdx,:,:);
        vLo = ch.W.v(loIdx,:,:);

        % Retrieve the frame size of the video 
        [~, nRows, nCols] = size(Vid);

        % Camera-space angles (radians)
        visualFieldPoints = anglesFromIntrinsics(nRows, nCols, fisheyeIntrinsics);

        % Warp to participant space using affine (apply in DEGREES to [az, el])
        az_deg = rad2deg(visualFieldPoints(:,1));
        el_deg = rad2deg(visualFieldPoints(:,2));

        azelW      = [az_deg(:), el_deg(:), ones(numel(az_deg),1)] * affineMat.';  % N×2 deg
        azW_deg    = reshape(azelW(:,1), size(az_deg));
        elW_deg    = reshape(azelW(:,2), size(el_deg));
        phi_w      = deg2rad(azW_deg);
        theta_w    = deg2rad(elW_deg);

        % High AS maps on participant-space sphere
        [hFigHighSlope, hFigHighIntercept, frq] = ... 
            mapSlopeIntSPD(vHi, fsVid, [40,40], 20, theta_w, phi_w, r);

        % Low AS maps on participant-space sphere
        [hFigLowSlope,  hFigLowIntercept, ~]  = ...
            mapSlopeIntSPD(vLo, fsVid, [40,40], 20, theta_w, phi_w, r);

        % All AS maps on participant-space sphere
        [hFigAllSlope,  hFigAllIntercept, ~]  = ...
            mapSlopeIntSPD(Vid, fsVid, [40,40], 20, theta_w, phi_w, r);



    end

    % find common limits (COMMENTED OUT FOR COORDINATE TRANSFORM USE)
    %{
    vminSlope = min( [slopeHigh(:); slopeLow(:); slopeAll(:)] );
    vmaxSlope = max( [slopeHigh(:); slopeLow(:); slopeAll(:)] );

    vminInt   = min( [interceptHigh(:); interceptLow(:); interceptAll(:)] );
    vmaxInt   = max( [interceptHigh(:); interceptLow(:); interceptAll(:)] );

    % re-apply caxis and add titles
    figure(hFigHighSlope);     caxis([vminSlope vmaxSlope]); title('1/f Slope Map - High AS', 'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white');
    figure(hFigLowSlope);      caxis([vminSlope vmaxSlope]); title('1/f Slope Map - Low AS',  'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white');
    figure(hFigAllSlope);      caxis([vminSlope vmaxSlope]); title('1/f Slope Map - All AS',  'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white');   

    figure(hFigHighIntercept); caxis([vminInt   vmaxInt  ]); title('1/f Intercept Map - High AS', 'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white');
    figure(hFigLowIntercept);  caxis([vminInt   vmaxInt  ]); title('1/f Intercept Map - Low AS',  'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white');
    figure(hFigAllIntercept);  caxis([vminInt   vmaxInt  ]); title('1/f Intercept Map - All AS',  'FontSize',16); set(gca,'FontSize',14); set(gcf,'color','white'); 
    %}
end 

% Local function to sort filenames from the directory 
function files = sort_chunk_filenames(directory)
    % Retrieve the chunk filepaths for the video  (unsorted)
    files = dir(fullfile(directory,'chunk_*.mat'));
    N_chunks = numel(files); 

    % Preallocate cell array for storing (number, file) pairs
    file_info = []; 

    % Iterate over the filepaths and find the number 
    % per chunk 
    for ii = 1:N_chunks
        % Attempt to find the number in the chunk name
        match = regexp(files(ii).name, '^chunk_(\d+)\.mat$', 'tokens');

        % If the match is empty, something's gone really wrong 
        if(isempty(match))
            error('Invalid file: %s', files(ii).name); 
        end 

        % If we did find a numeric match, retrieve the chunk number
        chunk_num = str2double(match{1}{1});

        % Save the filepath and its associated chunk number
        file_info(ii).number = chunk_num+1;
        file_info(ii).file = files(ii);
    end

    % Convert struct array to table, then sort by the chunk number
    % This ensures that the files are sorted numerically. 
    file_info = struct2table(file_info);
    file_info = sortrows(file_info, 'number');

    % Retireve sorted file structs
    files = file_info.file;

    % Clear unneeded information now that we have the files sorted 
    clear file_info; 

end 


% Local function to plot MS high/ligh 
function [yHigh, yLow]  = obtain_ms_high_low(files, N_chunks)
    % Collect all AS7341 values across all chunks
    AS7_all = [];
    AS7_all_t = []; 

    % GET GLOBAL AVG ACROSS CHUNKS
    for ii = 1:N_chunks
        % Load in the chunk 
        C = load(fullfile(files(ii).folder, files(ii).name), 'chunk');
        ch = C.chunk;
        
        % Retrieve all the readings of the 7th channel of the AS chip and concat 
        AS7 = ch.M.v.AS(:,7);
        AS7_t = ch.M.t.AS; 
        AS7_all = [AS7_all; AS7(:)];
        AS7_all_t = [AS7_all_t; AS7_t(:)]; 
    end

    % Find global min, max, and norm 
    as_min = min(AS7_all);
    as_max = max(AS7_all);
    as_norm = (AS7_all - as_min) / (as_max - as_min);

    % Find a global threshold from all of the AS readings 
    level_norm = graythresh(as_norm);  

    % convert back
    thr = level_norm*(as_max - as_min) + as_min;

    %% AS BRIGHTNESS OVER TIME (red = above thr, blue = below thr)
    figure; hold on;

    % Legend anchors (so legend doesn’t list every chunk)
    hHigh = plot(nan,nan,'r-','LineWidth',1.2,'DisplayName','High (>= thr)');
    hLow  = plot(nan,nan,'b-','LineWidth',1.2,'DisplayName','Low (< thr)');

    % Iterate over all the chunks 
    % Plot all the AS7 data across all chunks 
    % distinguish between high and low light level.
    
    % Split into above/below threshold (use NaNs to break lines)
    yHigh = AS7_all;              
    yHigh(AS7_all <= thr) = NaN;
    yLow  = AS7_all;              
    yLow(AS7_all  > thr) = NaN;

    % Now, we do not need the AS7 all information anymore, clear to save memory 
    clear AS7_all; 
    clear AS7_all_t; 
end 

% Local function to plot the MS high/low results 
function plot_ms_high_low(yHigh, yLow)
    figure; 

    % Plot; hide handles so legend doesn’t get cluttered
    plot(AS7_all_t, yHigh, 'r-', 'LineWidth', 1.0, 'HandleVisibility','off');
    hold on; 
    plot(AS7_all_t, yLow,  'b-', 'LineWidth', 1.0, 'HandleVisibility','off');

    % Axes, threshold line, labels, legend
    set(gca, 'YScale','log', 'FontSize',14, 'TickLength',[0.02 0.02]);
    yline(thr, 'k--', 'HandleVisibility','off');

    xlabel('Time (s)', 'FontSize',14);
    ylabel('AS channel 7 (brightness)', 'FontSize',14);
    title('AS Brightness over Time (log scale)', 'FontSize',16, 'FontWeight','normal');
    legend([hHigh hLow], 'Location','best');  % only two entries
    set(gcf, 'color','white');
    hold off;

end 


% Local function to obtain SPD data 
function  [globalHi, globalLo, globalCen, globalPer, f_int] = obtain_SPD_data(files, N_chunks)
    %% BUILD COMMON FREQUENCY GRID
    % Use a nominal 10 s block from chunk 1 to get its frequency bins
    tmp = load(fullfile(files(1).folder,files(1).name),'chunk');
    Vid0 = tmp.chunk.W.v;
    [~, frq0] = calcTemporalSPD( Vid0(1:round(2*Twin*fsVid),:,:), fsVid, 'lineResolution', false );
    f_min = min(frq0(frq0>0)); % bounds
    f_start = ceil(f_min);
    f_end = floor(fsVid / 2) - 1;
    f_int = (f_start : f_end)';

    % We no longer need tmp, so clear it to save memory
    clear tmp; 

    %% SLIDING‐WINDOW SPD SPLIT & AVERAGE ACROSS CHUNKS
    allAll = nan(N_chunks, numel(f_int));
    allHi = nan(N_chunks, numel(f_int));
    allLo = nan(N_chunks, numel(f_int));
    allCen = nan(N_chunks, numel(f_int));
    allPer = nan(N_chunks, numel(f_int));

    % Iterate over the chunks of the video 
    for ii = 1:N_chunks
        % Load in the given chunk 
        C     = load(fullfile(files(ii).folder,files(ii).name),'chunk');
        ch    = C.chunk;
        tAS   = ch.M.t.AS(:);
        AS7   = ch.M.v.AS(:,7);
        W_t   = ch.W.t(:);
        Vid   = ch.W.v;
        [nFrames, nRows, nCols] = size(Vid);
        fps = fsVid;
        
        % Create central and peripheral mask
        centerMask = false(nRows, nCols);
        centerMask(centerRows, centerCols) = true;
        peripheryMask = ~centerMask;
        
        % threshold & AS→frame interpolation
        AS_if   = interp1(tAS, AS7, W_t, 'pchip', NaN);
        t0s     = (W_t(1)+Twin):hop:(W_t(end)-Twin);
        
        % accumulators
        spdAll = [];
        spdHi  = [];
        spdLo  = [];
        spdCen = [];
        spdPer = [];
        
        % PSD window loop. Loop over each t0 in t0s
        for t0 = t0s
            idx = (W_t>=t0-Twin)&(W_t<=t0+Twin);
            if nnz(idx)<2
                continue;
            end
            as_val = mean(AS_if(idx));
            isHi = as_val > thr;
            fprintf('Chunk %d | t0 = %.2f | nnz(idx) = %d\n', i, t0, nnz(idx));
            % -------- Whole-frame SPD (existing Hi/Lo path) --------
            [P,f] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false);
            fprintf(' → PSD size = [%d %d], freq range = [%.2f %.2f]\n', size(P), min(f), max(f));
            keep = f>0 & f<=fMax;
            fprintf(' → # kept freqs = %d\n', nnz(keep));
            frqLoc = f(keep);
            spdLoc = P(keep);
            % 30HZ CENSOR
            spdLoc(frqLoc==30) = NaN;
            if numel(frqLoc) >= 2 && numel(spdLoc) >= 2
                Pi = interp1(frqLoc, spdLoc, f_int, 'pchip')';
                if numel(Pi) == numel(f_int)
                    spdAll(end+1, :) = Pi;
                    if isHi
                        spdHi(end+1,:) = Pi;
                    else
                        spdLo(end+1,:) = Pi;
                    end
                end
            end
            % -------- Center --------
            [Pc, fc] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false, ...
                                    'regionMatrix', double(centerMask));
            keepC = fc>0 & fc<=fMax;
            if nnz(keepC) >= 2
                PiC = interp1(fc(keepC), Pc(keepC), f_int, 'pchip')';
                PiC(f_int==30) = NaN;
                spdCen(end+1,:) = PiC;
            end
            % -------- Periphery --------
            [Pp, fp] = calcTemporalSPD(Vid(idx,:,:), fsVid, 'lineResolution', false, ...
                                    'regionMatrix', double(peripheryMask));
            keepP = fp>0 & fp<=fMax;
            if nnz(keepP) >= 2
                PiP = interp1(fp(keepP), Pp(keepP), f_int, 'pchip')';
                PiP(f_int==30) = NaN;
                spdPer(end+1,:) = PiP;
            end
        end
        % store this chunk's mean curves
        allAll(ii,:)    = mean(spdAll, 1, 'omitnan');
        % again for high light and low light
        allHi(ii,:)     = mean(spdHi, 1, 'omitnan');
        allLo(ii,:)     = mean(spdLo, 1, 'omitnan');
        % again for center and periphery
        allCen(ii,:)    = mean(spdCen, 1, 'omitnan');
        allPer(ii,:)    = mean(spdPer, 1, 'omitnan');
        
        fprintf('[Chunk %d] spdHi windows: %d | spdLo windows: %d\n', ...
        ii, size(spdAll, 1), size(spdHi, 1), size(spdLo, 1), size(spdCen, 1), size(spdPer, 1));
    end
    % global average across chunks
    globalAll = nanmean(allAll,1);
    globalHi = nanmean(allHi,1);
    globalLo = nanmean(allLo,1);
    globalCen = nanmean(allCen, 1);
    globalPer = nanmean(allPer, 1);
end 

% Local function to plot the SPD data 
function plot_SPD_data(globalAll, globalHi, globalLo, globalCen, globalPer, f_int)
    % PLOT GLOBAL SPD
    figure;
    plotSPD(globalAll, f_int);
    xlabel('Frequency (Hz)', 'FontSize',14);
    ylabel('Spectral power density (contrast^2/Hz)', 'FontSize',14);
    title('All Chunks: Avg Windowed SPD (1–90 Hz)', ...
        'FontSize',16, 'FontWeight','normal');
    xlim([f_int(1) f_int(end)]);
    nticks = 10;
    ticks = round(logspace(log10(0.2),log10(59),nticks), 2);
    xticks(ticks);
    xticklabels(ticks);
    set(gca, 'FontSize',14, 'TickLength',[0.02 0.02]);
    set(gcf, 'Color','white');

    % PLOT GLOBAL SPD (high vs low)
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
    xticks(ticks);
    xticklabels(ticks);
    legend('High AS','Low AS','Location','best');
    legend({'High AS','Low AS'}, 'Location','best');
    set(gca, 'FontSize',14, 'TickLength',[0.02 0.02]);
    set(gcf, 'Color','white');

    % PLOT GLOBAL SPD (center vs surround)
    figure;
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

end 


% Local function to obtain slope and intercept maps 
function obtain_slope_and_intercept_maps()

end