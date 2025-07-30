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

  %--- 4) log–log plot
  loglog(frq, spd, varargin{:});
  xlabel('Frequency (Hz)')
  ylabel('Spectral power density (contrast^2/Hz)')
  title('Temporal SPD')
end
