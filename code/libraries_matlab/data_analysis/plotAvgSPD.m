%% SCRIPT TO EXPRESS AVERAGE SPD ACROSS CHUNKS FOR EACH LMS OPERATION
%
% This script loads video chunks, calculates their individual Temporal SPD
% (for combined RGB, 200 FPS), averages these SPDs, and plots the result.


% Get number of frames
[nFrames, nRows, nCols] = size(v);
fps = 184.74;

% Create central and peripheral matrix
centerMatrix = zeros(nRows, nCols);
centerMatrix(121:360, 161:480) = 1;

peripheryMatrix = ones(nRows, nCols);
peripheryMatrix(121:360, 161:480) = 0;

% Filter empty/invalid chunks
good_chunks = cellfun(@(chunk) ~isempty(chunk) && any(chunk.W.v(:)), chunks);

% Define the post-receptoral channels to analyze
channels = {'LM', 'L-M', 'S'};
% Define colors for plotting (Note: users on light mode adjust 'w' to 'b')
channel_colors = containers.Map({'LM', 'L-M', 'S'}, {'w', 'r', 'b'});

% Prepare storage
avg_spds = cell(1, numel(channels));

    for c = 1:numel(channels)
    
        % Get SPD for center
        [center_spds, frqs] = cellfun(@(chunk) calcTemporalSPD(chunk.W.v(1:nFrames,:,:), fps, 'postreceptoralChannel', channels{c}, 'regionMatrix', centerMatrix), chunks(good_chunks), 'UniformOutput', false);
        % Get SPD for periphery
        [periphery_spds, ~] = cellfun(@(chunk) calcTemporalSPD(chunk.W.v(1:nFrames,:,:), fps, 'postreceptoralChannel', channels{c}, 'regionMatrix', peripheryMatrix), chunks(good_chunks), 'UniformOutput', false);
    
        % Average the SPDs across all chunks for center and periphery
        avg_c_spds{c} = mean(cat(2, center_spds{:}), 2)';
        avg_p_spds{c} = mean(cat(2, periphery_spds{:}), 2)';
    end

% Visualize the result (STRUCK THIS FOR NOW)
% figure;
% hold on;

for c = 1:numel(channels)
    figure;
    hold on;

    nonzero_idx = frqs{c} > 0;
    clean_frq = frqs{c}(nonzero_idx);

    clean_c_spd = avg_c_spds{c}(nonzero_idx);
    clean_p_spd = avg_p_spds{c}(nonzero_idx);

    %ref_c_curve = clean_c_spd(1) * (clean_frq(1) ./ clean_frq).^2;
    %ref_p_curve = clean_p_spd(1) * (clean_frq(1) ./ clean_frq).^2;

    plotSPD(clean_c_spd, clean_frq) %, 'Color', channel_colors(channels{c}));
    plotSPD(clean_p_spd, clean_frq) %, 'Color', channel_colors(channels{c}));
    
    %plotSPD(ref_c_curve, clean_frq, 'Color', channel_colors(channels{c}), 'LineStyle', '--');
    %plotSPD(ref_p_curve, clean_frq, 'Color', channel_colors(channels{c}), 'LineStyle', '--');

    legend({'Center', 'Periphery', 'Center 1/f²', 'Periphery 1/f²'}, 'Location', 'best');
    title(sprintf('Temporal SPD - %s Channel', channels{c}));
    xlabel('Temporal Frequency (Hz)');
    ylabel('Power (contrast² / Hz)');
    set(gca, 'XScale', 'log', 'YScale', 'log');

end

% (STRUCK FOR NOW)
% legend({'LM', '1/f² LM', 'L-M', '1/f² L-M', 'S', '1/f² S'});
% title('SPD and 1/f² Reference for LMS Channel Operations - IMX219 Camera');
% xlabel('Temporal Frequency (Hz)');
% ylabel('Power (contrast^2 / Hz)');
% set(gca, 'XScale', 'log', 'YScale', 'log');
