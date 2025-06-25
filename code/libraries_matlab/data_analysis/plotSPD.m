function plotSPD(spd, frq, varargin)
% Plots the SPD function for a chunk.
%
% Syntax:
%   [spd,frq] = PlotSPD(spd, frq)
%
% Description:
%   xxxx
%
% Inputs:
%   spd                   - 1xf. The spectral power density in units of
%                           contrast^2/Hz.
%   frq                   - 1xf. The frequencies (in Hz) for the spd value.
%   varargin              - Optional plotting style arguments, e.g.:
%                           'Color', [1 0 0], 'LineWidth', 2, etc.
%
% Examples:
%{  
    dataDir = getpref('lightLoggerAnalysis','dataBaseDir');
    test_path = fullfile(dataDir,;
    chunks = parse_chunks(test_path);
    [spd, frq] = calcTemporalSPD(chunks{1}.W.v, 200);
    plotSPD(spd, frq);
%} 

arguments 
        spd (1, :)
        frq (1, :)
end

arguments (Repeating)
        varargin
end


% Plot and label
spd = reshape(spd, 1, []);
frq = reshape(frq, 1, []);

% Fix any zero/negative frequencies
frq(frq <= 0) = eps;

%plot(log10(frq), log10(spsd), varargin{:})
loglog(frq, spd, varargin{:});
