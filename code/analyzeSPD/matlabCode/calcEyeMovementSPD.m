function [spd, frq] = calcEyeMovementSPD(eyePos, fs, options)
    arguments
        eyePos (:,2) {mustBeNumeric}
        fs (1,1) {mustBeNumeric} = 120
        options.windowDuration (1,1) {mustBeNumeric} = 5; % 5 seconds
        options.outputType (1,:) char {mustBeMember(options.outputType, {'magnitude', 'separate'})} = 'magnitude'
        options.nan_threshold (1,1) {mustBeNumeric} = 0.1;
    end

    % 1. Handle NaNs/Blinks
    nanFraction = sum(isnan(eyePos(:))) / numel(eyePos);
    if nanFraction >= options.nan_threshold
        frq = nan; spd = nan; return;
    end

    % 2. Process Signal Magnitude
    if strcmp(options.outputType, 'magnitude')
        centeredPos = eyePos - mean(eyePos, 1, 'omitnan');
        signal = sqrt(sum(centeredPos.^2, 2));
    else
        signal = eyePos; 
    end

    % 3. Sliding Window Logic
    windowSize = round(options.windowDuration * fs);
    nSamples = size(signal, 1);
    
    % If the signal is shorter than one window, just process the whole thing
    if nSamples <= windowSize
        [frq, spd] = simplePSD(double(signal - mean(signal, 'omitnan')), fs);
        return;
    end

    % Calculate how many full windows we have (no overlap for simplicity)
    nWindows = floor(nSamples / windowSize);
    spdSum = 0;

    for i = 1:nWindows
        idx = ((i-1)*windowSize + 1) : (i*windowSize);
        chunk = signal(idx, :);
        
        % Skip chunk if it contains too many NaNs
        if sum(isnan(chunk(:))) / numel(chunk) > options.nan_threshold
            continue;
        end
        
        % De-mean the chunk and calculate PSD
        chunk = chunk - mean(chunk, 1, 'omitnan');
        [frq, currentSpd] = simplePSD(double(chunk), fs);
        
        % Accumulate
        spdSum = spdSum + currentSpd;
    end
    
    % Average the SPDs
    spd = spdSum / nWindows;
end