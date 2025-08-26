function plotSPD(spd, frq, varargin)
% Plots the SPD function for a chunk.
%
% Syntax:
%   [spd, frq] = plotSPD(spd, frq)
%
% Description:
%   Plots the SPD function for a chunk.
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

    %--- 1) flatten to column vectors
    spd = spd(:);
    frq = frq(:);

    %--- 2) align lengths
    N = min(numel(spd), numel(frq));
    spd = spd(1:N);
    frq = frq(1:N);

    %--- 3) remove zero or negative freqs
    keep = frq > 0;
    spd  = spd(keep);
    frq  = frq(keep);

    %--- 3) plot the SPD (log-log)
    hMain = loglog(frq, spd, varargin{:});
    holdState = ishold; hold on;

    %--- 4) fit line on log10-log10 domain
    x = log10(frq);
    y = log10(spd);
    p = polyfit(x, y, 1);               % y ≈ p(1)*x + p(2)
    yfit = polyval(p, x);
    spdFit = 10.^yfit;

    %--- 5) plot dotted best-fit, same color as main line
    fitColor = get(hMain, 'Color');
    hFit = loglog(frq, spdFit, '--', 'Color', fitColor, 'LineWidth', 1.2, ...
                  'HandleVisibility','off');

    %--- 6) cosmetics
    xlabel('Frequency (Hz)');
    ylabel('Spectral power density (contrast^2/Hz)');
    title('Temporal SPD');

    if ~holdState, hold off; end
end
