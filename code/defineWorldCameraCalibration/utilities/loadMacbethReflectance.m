function macbethPatchReflectance = loadMacbethReflectance()
    % loadMacbethReflectance Loads spectral reflectance data for the Macbeth ColorChecker.
    %
    %   macbethPatchReflectance = loadMacbethReflectance() reads the
    %   'spectral_data' sheet from the provided BabelColor Excel file and
    %   returns a formatted MATLAB table.
    
    % Identify the location of the file
    filename = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'ColorChecker_RGB_and_spectra.xls');

    % 1. Read the raw data using readcell
    % This prevents issues with the metadata in the very first row of the sheet.
    rawData = readcell(filename, 'Sheet', 'spectral_data');

    % 2. Extract patch names
    % The names are located in rows 3 through 26, in the 2nd column.
    patchNames = string(rawData(3:26, 2));

    % Convert patch names to valid MATLAB variable names 
    % (e.g., "dark skin" becomes "dark_skin")
    validNames = matlab.lang.makeValidName(patchNames);

    % 3. Extract the wavelengths
    % The wavelengths (380 to 730) are located in row 2, starting at column 3.
    % There are 36 wavelength bins (380:10:730).
    wavelengths = cell2mat(rawData(2, 3:38))';

    % 4. Extract the reflectance data and transpose
    % The numeric data is in rows 3 through 26, spanning columns 3 through 38.
    % We transpose it (') so wavelengths are rows and patches are columns.
    reflectances = cell2mat(rawData(3:26, 3:38))';

    % 5. Build the final formatted table
    macbethPatchReflectance = array2table(reflectances, 'VariableNames', validNames);

    % Append the wavelength vector and move it to be the first column
    macbethPatchReflectance.Wavelength = wavelengths;
    macbethPatchReflectance = movevars(macbethPatchReflectance, 'Wavelength', 'Before', 1);
end