function [meta, signal] = extractLabSphereData(filename)
% Loads a Labsphere CSV, prompts user to pick a timestamp by number,
% and returns key metadata fields in a struct plus time-series data for that scan.
%
% Syntax:
%   [meta, signal] = extractLabSphereData('Labsphere_SaveAllScans.csv');
%
% Description: 
%   This function loads the CSV data from one or more Labsphere scan(s), 
%   lists their available timestamps, and prompts the user for a selection.
%   Extracts sample rate, min/max values, frequency, and flicker metrics
%   into the output struct. Also reads the full time-series data for the 
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
rowNames = string(C(:,1));
rawData  = C(:,2:end);
nScans   = size(rawData,2);

% Get human‐readable timestamps
tsRow   = find(rowNames == "TimeStamp:", 1);
stamps  = string(rawData(tsRow,:));

% 4) Print available scan options for user
fprintf('Available scans:\n');
for i = 1:nScans
  fprintf('  %2d: %s\n', i, stamps(i));
end

% Prompt user to select scan number
sel = input(sprintf('Enter scan number (1–%d): ', nScans));
if ~isscalar(sel) || sel<1 || sel>nScans || sel~=round(sel)
  error('Invalid scan number. Must be an integer between 1 and %d.', nScans);
end

% Define desired data rows and their field names
  rows = { ...
    'Sample Rate (Hz):'
    'Minimum Value:'
    'Maximum Value:'
    'Fundamental Frequency (Hz):'
    'Percent Flicker:'
  };
  fieldNames = { ...
    'SampleRateHz'
    'MinimumValue'
    'MaximumValue'
    'FundamentalFrequencyHz'
    'PercentFlicker'
  };

% Initialize output struct and extract each metadata value
meta = struct();

for k = 1:numel(rows)
  idx = find(rowNames == rows{k},1);
  raw = rawData{idx, sel};
  meta.(fieldNames{k}) = raw;
end

% Extract time series data and convert into signal vector
tsDataRow = find(rowNames == "Time Series Data:", 1);
  seriesCells = rawData(tsDataRow+1:end, sel);
  N = numel(seriesCells);
  signal = zeros(N,1);
  for i = 1:N
    raw = seriesCells{i};
    signal(i) = raw;
  end

% Trim trailing NaN values for shorter scans. 
validIdx = find(~isnan(signal));
signal = signal(1:validIdx(end));

end