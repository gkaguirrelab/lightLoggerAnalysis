function [meta, signal] = extractLabSphereData(filename,sel)
% Loads a Labsphere CSV, prompts user to pick a timestamp by number,
% and returns key metadata fields in a struct plus time-series data for that scan.
%
% Syntax:
%   [meta, signal] = extractLabSphereData('Labsphere_SaveAllScans.csv');
%
% Description:
%   This function loads the CSV data from one or more Labsphere scan(s),
%   lists their available timestamps, and prompts the user for a selection.
%   Extracts sample rate, frequency, and flicker metrics into the output
%   struct. Also reads the full time-series data for the
%   chosen scan and converts into a numeric 'signal' vector to be used in
%   determining temporal SPD in units of contrast.
%
% Inputs:
%   filename            - String. CSV file containing Labsphere scan data.
%
% Outputs:
%   meta                - Struct. The desired metadata contained in the CSV
%                         file. This can be manipulated to include/exclude
%                         specific metadata.
%   signal              - Numeric vector.
%
% Examples:
%{
    [meta, signal] = extractLabSphereData('Labsphere_SaveAllScans.csv');

    signal(:,1) = (signal(:,1) - mean(signal(:,1)))/mean(signal(:,1));
    [frq, spd] = simplePSD(signal(:,1), meta.SampleRateHz);
    
    figure;
    loglog(frq(2:end), spd(2:end));
    xlabel('Frequency (Hz)')
    ylabel('Power/Hz')
    title('PSD (loglog)')
%}

% Read CSV into cell array
C = readcell(filename, 'Delimiter', ',');

% Separate row names from raw data
rawData  = C(:,2:end);
nScans   = size(rawData,2);

% Get human‐readable timestamps
stamps = string(rawData(2,:));
comment = string(rawData(8,:));

% Display scans with timestamp and comment
if nargin == 1
    fprintf('Available scans:\n');
    for i = 1:nScans
        if strlength(comment(i)) > 0
            fprintf('  %2d: %s  — %s\n  ', i, char(stamps(i)), char(comment(i)));
        else
            fprintf('  %2d: %s\n  ', i, char(stamps(i)));
        end
    end

    % Prompt user to select scan number
    sel = input(sprintf('Enter scan number (1–%d): ', nScans));
    if ~isscalar(sel) || sel<1 || sel>nScans || sel~=round(sel)
        error('Invalid scan number. Must be an integer between 1 and %d.', nScans);
    end
end

% Initialize output struct and extract desired metadata values
meta = struct();
meta.Comment = rawData{8, sel};
meta.SampleRateHz = rawData{10, sel};
meta.FundamentalFrequencyHz = rawData{29, sel};
meta.PercentFlicker = rawData{30, sel};
% add scale max value (row 12)

% Extract time series data and convert into signal vector
seriesCells = rawData(37:end, sel);
N = numel(seriesCells);
signal = zeros(N,1);
for i = 1:N
    raw = seriesCells{i};
    signal(i) = raw;
end

% Trim trailing NaN values for shorter scans.
validIdx = find(~isnan(signal));
signal = signal(1:validIdx(end));

% Convert to contrast units
signal = (signal - mean(signal))/mean(signal);

end