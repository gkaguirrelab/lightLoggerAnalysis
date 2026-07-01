function [spd, frq] = calcEyeMovementSPD(eyePos, fs, options)
% Compute the spectral power density of eye movement position data
%
% Syntax:
%   [spd, frq] = calcEyeMovementSPD(eyePos)
%   [spd, frq] = calcEyeMovementSPD(eyePos, fs)
%   [spd, frq] = calcEyeMovementSPD(eyePos, fs, options)
%
% Description:
%   Computes the power spectral density (PSD) of an eye position time
%   series using non-overlapping sliding windows. The input eye position
%   signal can be returned as either a scalar magnitude (Euclidean
%   distance from the centered mean) or kept as separate x/y channels.
%   If the signal is shorter than one window, the entire signal is
%   processed as a single segment. Segments with NaN fractions exceeding
%   the threshold are skipped, and if the overall NaN fraction is too
%   high the function returns NaN immediately.
%
% Inputs:
%   eyePos                - Numeric matrix (Nx2). Eye position samples
%                           with columns [x, y].
%   fs                    - Scalar. Sampling rate in Hz. Defaults to 120.
%
% Optional key/value pairs:
%   windowDuration        - Scalar. Duration of each analysis window in
%                           seconds. Defaults to 5.
%   outputType            - Char. Either 'magnitude' (default) to compute
%                           the PSD of the Euclidean magnitude, or
%                           'separate' to keep x and y channels as-is.
%   nan_threshold         - Scalar. Fraction of NaN values (0 to 1) above
%                           which the signal or a window is skipped.
%                           Defaults to 0.1.
%
% Outputs:
%   spd                   - Numeric vector. Spectral power density
%                           averaged across windows.
%   frq                   - Numeric vector. Corresponding frequency
%                           values in Hz.
%
% Examples:
%{
    eyePos = randn(6000, 2);
    fs = 120;
    [spd, frq] = calcEyeMovementSPD(eyePos, fs);
    loglog(frq, spd);
%}

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