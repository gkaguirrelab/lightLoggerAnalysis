function figHandle = timeFreqSpectrogramFromChunks( chunks, varargin )
% Take a set of light logger chunks and create a summary plot
%
% Syntax:
%   timeFreqSpectrogramFromChunks( chunks )
%
% Description:
%   More to say
%
% Inputs:
%   chunks                - Cell array. Foo foo foo foo foo foo foo foo foo
%                           foo foo foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo
%
% Optional key/value pairs:
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
    load('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_data/LightLogger_Recordings/quickOutdoorIndoor.mat');
    timeFreqSpectrogramFromChunks( chunks);
%}

% Parse the parameters
p = inputParser; p.KeepUnmatched = false;
p.addParameter('fps',200,@isnumeric);
p.addParameter('windowDurSecs',4,@isnumeric);
p.addParameter('windowStepSecs',1,@isnumeric);
p.addParameter('hallEffectThreshold',1350,@isnumeric);
p.addParameter('savePDF',true,@islogical);
p.addParameter('figureSaveDir','~/Desktop',@ischar);
p.parse(varargin{:})

% Extract the parameters
fps = p.Results.fps;
windowDurSecs = p.Results.windowDurSecs;
windowStepSecs = p.Results.windowStepSecs;

% Setup a figure
figHandle = figure('color', 'w');
figuresize(400,600,'pt');
t = tiledlayout(6,1);
t.TileSpacing = 'tight';
t.Padding = 'none';

% Get the first value of the world time; all plots will be with respect to
% this
startTime = chunks{1}.W.t(1);

%% Plot the sensor counts for one MS channel
msCounts = [];
msTime = [];
channel = 7;

% Loop over chunks and concatenate the measures
for cc = 1:length(chunks)
    t = chunks{cc}.M.t;
    v = chunks{cc}.M.v.AS;
    msCounts = [msCounts v(:,channel)'];
    msTime = [msTime t];
end

% Plot the data
nexttile
plot(msTime-startTime,log10(msCounts),'-r','LineWidth',2);
ylabel('log counts');
ylim([0 6]);
xlim([0 600]);
a = gca;
a.TickDir = 'out';
a.XTick = [];
box off

%% Plot the sunglasses detector
shadeCounts = [];
shadeTime = [];

% Loop over chunks and concatenate the measures
for cc = 1:length(chunks)
    t = chunks{cc}.S.t;
    v = chunks{cc}.S.v;
    shadeCounts = [shadeCounts v];
    shadeTime = [shadeTime t];
end

% Apply a threshold to binarize the sensor reading into sunglasses on or
% off. The property of the Hall effect sensor is that the sensor reading
% decreases in value in the presence of the magnets of the clip on shades.
shadesOn = shadeCounts<p.Results.hallEffectThreshold;

% Plot the data
nexttile
plot(shadeTime-startTime,shadesOn,'-k','LineWidth',2)
ylim([-0.5 1.5]);
a=gca();
a.YTick=[0 1];
a.YTickLabel = {'off','on'};
xlim([0 600]);
a = gca;
a.TickDir = 'out';
a.XTick = [];
box off

%% Calculate world camera spectral power density 

% This is an index into the time-frequency spectrogram matrix we will build
lineIdx = 1;

% These are indices for the accumuluation of high and low ambience spds
highCount = 0;
lowCount = 0;

% Loop over chunks
for cc = 1:length(chunks)

    % Get the time base and signal for this chunk
    t =chunks{cc}.W.t;
    v =chunks{cc}.W.v;

    % Check for an empty chunk
    if isempty(t)
        continue
    end

    % Account for time elapsed between chunks in the lineIdx
    if cc>1
        lineIdx = lineIdx + floor( (t(1) - spectotime(end)) / windowStepSecs);
    end

    % Slide a window across the time series and obtain the spectral density
    frameIdx = 1;
    notDone = true;
    while notDone
        endIdx = frameIdx + windowDurSecs*fps;
        if endIdx > length(v)
            notDone = false;
            continue
        end

        % FFT of the signal in contrast units
        signal = v(frameIdx:endIdx);
        signal = (signal - mean(signal))/mean(signal);
        [frq,amp] = simpleFFT(signal,fps);

        % Discard the zero and nyquit frequencies
        frq = frq(2:end-1); amp = amp(2:end-1);

        % Covert from amplitude (contrast) to spectral power density
        % (contrast^2/Hz)
        spd=(amp.^2)./frq;

        % Sum the spds across the windows for later plotting. Do this
        % separately for the high and low ambient illuminance values
        if lineIdx == 1
            avgSPDLow = zeros(size(spd));
            avgSPDHigh = zeros(size(spd));
        end
        [~,msTimeIdx]=min(abs(msTime-t(frameIdx)));
        if log10(msCounts(msTimeIdx))>2.5
            avgSPDHigh = avgSPDHigh + spd;
            highCount = highCount + 1;
        else
            avgSPDLow = avgSPDLow + spd;
            lowCount = lowCount + 1;
        end

        % Store the data and time point
        spectrogram(lineIdx,:) = spd;
        spectotime(lineIdx) = t(frameIdx);

        % Increment the indices
        lineIdx = lineIdx+1;
        frameIdx = frameIdx+windowStepSecs*fps;
    end
end

% Convert the summed SPD into the average
avgSPDHigh = avgSPDHigh/highCount;
avgSPDLow = avgSPDLow/lowCount;

%% Plot the time-frequency spectrogram

% Some hard-coded properties of the plot
axisVals = [0.5,1,3,6,12,25,50,100];
limVals = [-10 0];

% Prepare the spectrogram matrix
Im = spectrogram';
zeroIdx = Im==0;
Im = log10(Im);
Im(zeroIdx)=nan;

% Threshold the matrix for extreme values
Im(Im< limVals(1))=limVals(1); Im(Im>limVals(2))=limVals(2);

% Plot the spectrogram over two rows
nexttile([2,1]);
if p.Results.savePDF
    plot([1,size(Im,2)],[1,size(Im,1)],'-r');
else
    contourf(Im,25,'LineStyle','none');
end

% Clean up the ticks and labels
a = gca;
a.YScale='log';
a.TickDir = 'out';
colormap(turbo)
clim(limVals);
for ff = 1:length(axisVals)
    [~,yTickVals(ff)] = min(abs(axisVals(ff)-frq));
end
a.YTick = yTickVals;
a.YTickLabel = arrayfun(@(x) {num2str(x)},axisVals);
ylim([1 length(frq)+1]);
oneMinute = 60*windowStepSecs;
xlim([1,600]);
a.XTick = 1:oneMinute:oneMinute*ceil((size(Im,2)/oneMinute));
a.XTickLabel = arrayfun(@(x) {num2str(x)},0:1:(1+length(a.XTick)));
ylabel('Freq [Hz]');
xlabel('time [minutes]')
a.YMinorTick = 'off';
a.TickDir = 'out';
box off

%% Show the color bar for the spectrogram

nexttile;
hCB = colorbar('north','AxisLocation','in');
hCB.Label.String = 'Power [log_1_0(contrast^2/Hz)]';
hCB.Ticks = [0 0.25 0.5 0.75 1];
hCB.XTickLabel = arrayfun(@(x) {num2str(x)},linspace(limVals(1),limVals(2),length(hCB.Ticks)));
set(gca,'Visible',false)
a = gca;
a.TickDir = 'out';
a.XTick = [];
colormap(turbo);
axis off


%% Plot the average SPD

nexttile()
plot(log10(frq),log10(avgSPDHigh),'-.r','LineWidth',2);
hold on
plot(log10(frq),log10(avgSPDLow),'-k','LineWidth',1);

% Add a 1/f^2 reference line
vec = log10(1./(frq.^2));
vec = vec-vec(1)+log10(avgSPDHigh(1));
plot(log10(frq),vec,':k')

% Clean up ticks and labels
ylabel({'Power','[log_1_0(contrast^2/Hz)]'});
xlabel('Frequency [Hz]');
xlim(log10([0.25 100]));
a = gca();
a.XTick = log10(axisVals);
a.XTickLabel = arrayfun(@(x) {num2str(x)},axisVals);
ylim(limVals);
box off

% Save the figure


if p.Results.savePDF
    figurePath = fullfile(p.Results.figureSaveDir,'demoSpectrogram.pdf');
    export_fig(figHandle,figurePath,'-Painters','-transparent');
else
    figurePath = fullfile(p.Results.figureSaveDir,'demoSpectrogram.png');
    export_fig(figHandle,figurePath,'-r600','-opengl');
end

end