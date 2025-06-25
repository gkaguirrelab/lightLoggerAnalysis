%% SCRIPT TO EXPRESS AVERAGE SPD ACROSS CHUNKS FOR EACH LMS OPERATION
%
% This script loads video chunks, calculates their individual Temporal SPD
% (for combined RGB, 200 FPS), averages these SPDs, and plots the result.


% Get number of frames
[nFrames, ~] = size(v);
fps = 184.74;

% Filter empty/invalid chunks
good_chunks = cellfun(@(chunk) ~isempty(chunk) && any(chunk.W.v(:)), chunks);

% Define the post-receptoral channels to analyze
channels = {'LM', 'L-M', 'S'};
% Define colors for plotting (Note: users on light mode adjust 'w' to 'b')
channel_colors = containers.Map({'LM', 'L-M', 'S'}, {'w', 'r', 'b'});

% Prepare storage
avg_spds = cell(1, numel(channels));

for c = 1:numel(channels)
    channel_name = channels{c};

    % Apply calcTemporalSPD to each video data element.
    [spds, frqs] = cellfun(@(chunk) calcTemporalSPD(chunk.W.v(1:nFrames,:,:), fps, 'postreceptoralChannel', channel_name), chunks(good_chunks), 'UniformOutput', false);

    % Average the SPDs across chunks
    avg_spds{c} = mean(cat(2, spds{:}), 2)';
end

% Visualize the result 
figure;
hold on;

for c = 1:numel(channels)
    nonzero_idx = frqs{c} > 0;
    clean_frq = frqs{c}(nonzero_idx);
    clean_spd = avg_spds{c}(nonzero_idx);
    plotSPD(clean_spd, clean_frq, 'Color', channel_colors(channels{c}));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
legend(channels);
title('Average Temporal SPD for LMS Channel Operations');
xlabel('Temporal Frequency (Hz)');
ylabel('Power (contrast^2 / Hz)');