function [spectralReflectance,spectralReflectanceS] = loadMacbethReflectance()
% loadMacbethReflectance Loads spectral reflectance data for the Macbeth ColorChecker.
% Returns a 6 x 4 cell array corresponding to the columns and rows of the chart, respectively.

% Identify the location of the file
filename = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'ColorChecker_RGB_and_spectra.xls');

% 1. Read the raw data using readcell
rawData = readcell(filename, 'Sheet', 'spectral_data');

% 2. Extract patch names and reflectance data
patchNames = string(rawData(3:26, 2));
reflectances = cell2mat(rawData(3:26, 3:38)); % 24 patches x 36 wavelength bins

% 3. Extract the wavelengths
spectralReflectanceS = WlsToS(cell2mat(rawData(2, 3:38))');

% 3. Hardcode the standard Macbeth ColorChecker chart layout (4 rows x 6 columns)
chartLayout = {
    "dark skin",    "light skin",    "blue sky",     "foliage",     "blue flower",   "bluish green";
    "orange",       "purplish blue",  "moderate red", "purple",      "yellow green",  "orange yellow";
    "blue",         "green",         "red",          "yellow",      "magenta",       "cyan";
    "white 9.5 (.05 D)",    "neutral 8 (.23 D)"     "neutral 6.5 (.44 D)",    "neutral 5 (.70 D)",   "neutral 3.5 (1.05 D)",   "black 2 (1.5 D)"
    };

% 4. Initialize the 6 x 4 cell array (6 columns, 4 rows)
spectralReflectance = cell(6, 4);

% 5. Map chart positions to the reflectance data
for r = 1:4
    for c = 1:6
        patchName = chartLayout{r, c};

        % Find the matching row index in the raw data
        idx = find(strcmpi(patchNames, patchName));

        if ~isempty(idx)
            % Store the reflectance spectrum in the [col, row] position
            spectralReflectance{c, r} = reflectances(idx, :)';
        else
            error('foo');
        end
    end
end
end