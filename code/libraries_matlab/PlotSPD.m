function PlotSPD(spd, frq)
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
%
% Examples:
%{  
%   test_path = "./mychunks";
%   chunks = parse_chunks(test_path);
%   [spd, frq] = calcTemporalSPD(chunks{1}.W.v, 200);
%   PlotSPD(spd, frq);
%} 

arguments 
        spd (2999, 1)
        frq (1, 2999)
       
end

% Plot and label 
    figure;
    loglog(frq, spd);
    xlabel('Frequency (Hz)');
    ylabel('Spectral Power Density (contrast^2/Hz)');
    title('Temporal SPD');

end
