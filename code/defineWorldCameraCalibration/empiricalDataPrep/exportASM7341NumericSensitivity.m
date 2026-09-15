function outputFileName = exportASM7341NumericSensitivity()
% Export the AS7341 spectral sensitivity as a plain numeric MAT file.
%
% data/ASM7341_spectralSensitivity.mat stores the sensitivity as a MATLAB
% table object, which scipy.io.loadmat cannot read. The Python minispect
% code (code/library/sensor_utility/ms_util.py) therefore reads a numeric
% companion file instead. That companion is a generated artifact and is
% listed in .gitignore, so run this function once after cloning the project
% to recreate it.
%
% The exported matrix is max-normalized per channel and transposed to
% nChannels x nWavelengths so that it is identical to the miniSpectT
% variable built inside estimateRadianceSpectrumFromMinispect.m.

projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');

sourceFileName = fullfile(projectRoot, 'data', 'ASM7341_spectralSensitivity.mat');
load(sourceFileName, 'T');

miniSpectWls = T.wl;
miniSpectT = table2array(T(:, ["F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "Clear" "NIR"]))';

% Sensitivity is max-normalized, matching the calibration dot product.
miniSpectT = miniSpectT ./ max(miniSpectT')';

outputFileName = fullfile(projectRoot, 'data', 'ASM7341_spectralSensitivity_numeric.mat');

% Save in a scipy-readable format rather than the HDF5-based v7.3 format.
save(outputFileName, 'miniSpectWls', 'miniSpectT', '-v7');

fprintf('Wrote %s (miniSpectT %dx%d, miniSpectWls %dx1)\n', ...
    outputFileName, size(miniSpectT, 1), size(miniSpectT, 2), numel(miniSpectWls));

end
