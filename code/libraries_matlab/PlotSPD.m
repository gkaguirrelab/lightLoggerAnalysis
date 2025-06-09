chunks = parse_chunks(pathtoexperiment, pathtolightloggerlibraries, pathtopiutil);
[spd, frq] = calcTemporalSPD(chunks{1}.W.v, 200);

%% PlotSPD
figure;
loglog(frq, spd);
xlabel('Frequency (Hz)');
ylabel('Spectral Power Density (contrast^2/Hz)');
title('Temporal SPD for Chunk 1');
