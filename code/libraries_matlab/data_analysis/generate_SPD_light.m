%% FUNCTION TO OBTAIN SPD MAP INFO FOR/AND MINISPECT LIGHT LEVELS FROM RECORDING ACROSS MULTIPLE CHUNKS
% 
% Adjust settings as necessary.
function maps = generate_SPD_light(directory, analyses_to_perform, visualize_results)

    arguments
        directory {mustBeText} = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_data/lightLogger/HERO_sm/sophia_in_wild_7-22/sophia_in_wild_7-22_chunks'
        analyses_to_perform (1,3) logical = [false, true, false]
        visualize_results   (1,3) logical = [false, true,  false]
    end

    % Define some constants about our analysis
    fsVid     = 120;      % world‐camera frame rate (Hz)
    Twin      = 5;        % half‐window (s) → 10 s total
    hop       = 2.5;      % shift between windows (s)
    fMax      = 90;       % highest frequency to plot (Hz)
    centerRows = 121:360; 
    centerCols = 161:480; 

    % predictable output shape
    maps = struct('slope', struct('highAS', [], 'lowAS', [], 'allAS', []), ...
              'intercept', struct('highAS', [], 'lowAS', [], 'allAS', []));

    % Sort the chunks in the directory numerically 
    files = sort_chunk_filenames(directory); 
    N_chunks = numel(files);

    % Ensure if we want to run slope intercept maps, 
    % we also have done MS calculation 
    if analyses_to_perform(end) && ~analyses_to_perform(1)
        error("Must perform MS analysis to do slope map analysis.");
    end 

    % Analyze the MS-AS data over the course of the video
    thr = NaN; t_all = [];
    if analyses_to_perform(1) 
        [yHigh, yLow, thr, t_all] = obtain_ms_high_low(files, N_chunks);
        if(visualize_results(1))
            plot_ms_high_low(yHigh, yLow, t_all, thr);
        end
    end
    if isnan(thr)  % Ensure we have a threshold for hi/lo split
        [~,~,thr] = obtain_ms_high_low();
    end 
    
    % Obtain the SPD data over all the chunks 
    if analyses_to_perform(2)
        [globalAll, globalHi, globalLo, globalCen, globalPer, f_int] = ...
            obtain_SPD_data(files, N_chunks, Twin, hop, fsVid, fMax, centerRows, centerCols, thr); 
        if(visualize_results(2))
            plot_SPD_data(globalAll, globalHi, globalLo, globalCen, globalPer, f_int); 
        end 
    end 

    maps = obtain_slope_and_intercept_maps();
 
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
    function [yHigh, yLow, thr, t_all] = obtain_ms_high_low(files, N_chunks, channelIdx)
        if nargin < 3, channelIdx = 7; end
        AS_all  = [];
        t_all   = [];
        
        for ii = 1:N_chunks
            C  = load(fullfile(files(ii).folder, files(ii).name), 'chunk');
            ch = C.chunk;
            AS = ch.M.v.AS(:, channelIdx);
            t  = ch.M.t.AS(:);
            AS_all = [AS_all; AS(:)];
            t_all  = [t_all;  t(:)];
        end
        
        % global Otsu threshold in ORIGINAL SCALE
        as_min   = min(AS_all);
        as_max   = max(AS_all);
        as_range = max(eps, as_max - as_min);
        as_norm  = (AS_all - as_min) / as_range;
        level    = graythresh(as_norm);
        thr      = level * as_range + as_min;
        
        % split into high/low traces
        yHigh = AS_all;  yHigh(AS_all <= thr) = NaN;
        yLow  = AS_all;  yLow(AS_all  > thr)  = NaN;
    end 
    
    % Local function to plot the MS high/low results 
    function plot_ms_high_low(yHigh, yLow)
        figure; hold on;
    
        % Plot; hide handles so legend doesn’t get cluttered
        plot(t_all, yHigh, 'r-', 'LineWidth', 1.0); 
        plot(t_all, yLow,  'b-', 'LineWidth', 1.0);
    
        % Axes, threshold line, labels, legend
        set(gca, 'YScale','log', 'FontSize',14, 'TickLength',[0.02 0.02]);
        yline(thr, 'k--');
    
        xlabel('Time (s)', 'FontSize',14);
        ylabel('AS channel 7 (brightness)', 'FontSize',14);
        title('AS Brightness over Time (log scale)', 'FontSize',16, 'FontWeight','normal');
        legend([hHigh hLow], 'Location','best');
        set(gcf, 'color','white');
        hold off;
    
    end 
    
    % Local function to obtain SPD data 
    function [globalAll, globalHi, globalLo, globalCen, globalPer, f_int] = ...
            obtain_SPD_data(files, N_chunks, Twin, hop, fsVid, fMax, centerRows, centerCols, thr)
        % BUILD COMMON FREQUENCY GRID
        % Use a nominal 10 s block from chunk 1 to get its frequency bins
        tmp = load(fullfile(files(1).folder,files(1).name),'chunk');
        Vid0 = tmp.chunk.W.v;
        [~, frq0] = calcTemporalSPD( Vid0(1:round(2*Twin*fsVid),:,:), fsVid, 'lineResolution', false );
        f_min = min(frq0(frq0>0)); % bounds
        f_start = ceil(f_min);
        f_end = floor(fsVid / 2) - 1;
        f_int = (f_start : f_end)';
        nFreq = numel(f_int);
    
        % We no longer need tmp, so clear it to save memory
        clear tmp; 
    
        % SLIDING‐WINDOW SPD SPLIT & AVERAGE ACROSS CHUNKS
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
            [~, nRows, nCols] = size(Vid);

            % Create central and peripheral mask
            centerMask = false(nRows, nCols);
            centerMask(centerRows, centerCols) = true;
            peripheryMask = ~centerMask;
            
            % threshold & AS→frame interpolation
            AS_if   = interp1(tAS, AS7, W_t, 'pchip', NaN);
            t0s     = (W_t(1)+Twin):hop:(W_t(end)-Twin);
            
            % accumulators
            spdAll = zeros(0, nFreq);
            spdHi  = zeros(0, nFreq);
            spdLo  = zeros(0, nFreq);
            spdCen = zeros(0, nFreq);
            spdPer = zeros(0, nFreq);
            
            % PSD window loop. Loop over each t0 in t0s
            for t0 = t0s
                idx = (W_t>=t0-Twin)&(W_t<=t0+Twin);
                if nnz(idx)<2
                    continue;
                end
                as_val = mean(AS_if(idx));
                isHi = as_val > thr;
                fprintf('Chunk %d | t0 = %.2f | nnz(idx) = %d\n', ii, t0, nnz(idx));
                
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
    function mapsOut = obtain_slope_and_intercept_maps()
        % ===== Running means (lazy init after first chunk) ===== *********
        meanSlopeHigh = [];  cntHigh = [];
        meanIntHigh   = [];
        meanSlopeLow  = [];  cntLow  = [];
        meanIntLow    = [];
        meanSlopeAll  = [];  cntAll  = [];
        meanIntAll    = [];

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

            % --- HIGH ---
            if any(hiIdx)
                vHi = Vid(hiIdx,:,:);
                [sH, iH] = mapSlopeIntSPD(vHi, fsVid, [40 40], 20, 'doPlot', false);
                [meanSlopeHigh, cntHigh] = updateMean(meanSlopeHigh, cntHigh, sH, size(vHi,1));
                [meanIntHigh,   ~      ] = updateMean(meanIntHigh,   cntHigh, iH, size(vHi,1)); % use same count
                clear vHi sH iH
            end
    
            % --- LOW ---
            if any(loIdx)
                vLo = Vid(loIdx,:,:);
                [sL, iL] = mapSlopeIntSPD(vLo, fsVid, [40 40], 20, 'doPlot', false);
                [meanSlopeLow, cntLow] = updateMean(meanSlopeLow, cntLow, sL, size(vLo,1));
                [meanIntLow,   ~     ] = updateMean(meanIntLow,   cntLow, iL, size(vLo,1));
                clear vLo sL iL
            end
    
            % --- ALL ---
            [sA, iA] = mapSlopeIntSPD(Vid, fsVid, [40 40], 20, 'doPlot', false);
            [meanSlopeAll, cntAll] = updateMean(meanSlopeAll, cntAll, sA, size(Vid,1));
            [meanIntAll,   ~     ] = updateMean(meanIntAll,   cntAll, iA, size(Vid,1));
            clear sA iA Vid C ch
        end

        % Final maps (avoid div-by-zero)
        mapsOut.slope.highAS   = finalizeMean(meanSlopeHigh, cntHigh);
        mapsOut.slope.lowAS    = finalizeMean(meanSlopeLow,  cntLow);
        mapsOut.slope.allAS    = finalizeMean(meanSlopeAll,  cntAll);
        mapsOut.intercept.highAS = finalizeMean(meanIntHigh, cntHigh);
        mapsOut.intercept.lowAS  = finalizeMean(meanIntLow,  cntLow);
        mapsOut.intercept.allAS  = finalizeMean(meanIntAll,  cntAll);
    end

    function [meanM, countM] = updateMean(meanM, countM, newM, w)
        if isempty(newM), return; end
        if isempty(meanM), meanM=zeros(size(newM)); countM=zeros(size(newM)); end
        valid = ~isnan(newM);
        countM(valid) = countM(valid) + w;
        delta = zeros(size(newM)); delta(valid) = newM(valid) - meanM(valid);
        meanM(valid) = meanM(valid) + (w .* delta(valid)) ./ countM(valid);
    end
    
    function out = finalizeMean(meanM, countM)
        out = meanM; out(countM==0) = NaN;
    end

end