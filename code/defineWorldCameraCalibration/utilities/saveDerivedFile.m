function saveDerivedFile(outputPath, README, derivedVariables)
% Save a derived MAT artifact with required human-readable metadata.
%
% Every programmatically generated file in derived/ must be written through
% this function. README must be nonempty, and derivedVariables must be a
% scalar struct whose fields become top-level MAT variables.

arguments
    outputPath {mustBeTextScalar}
    README {mustBeTextScalar}
    derivedVariables (1, 1) struct
end

outputPath = char(outputPath);
README = char(README);
assert(strlength(strtrim(string(README))) >= 20, ...
    'saveDerivedFile:InvalidREADME', ...
    'README must provide a substantive description of the derived data.');

variableNames = fieldnames(derivedVariables);
assert(~isempty(variableNames), ...
    'saveDerivedFile:NoDerivedVariables', ...
    'At least one derived variable must be supplied.');
assert(~ismember('README', variableNames), ...
    'saveDerivedFile:ReservedVariableName', ...
    'README is reserved for the description supplied separately.');

outputDirectory = fileparts(outputPath);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

payload = derivedVariables;
payload.README = README;
temporaryPath = [tempname(outputDirectory), '.mat'];
temporaryCleanup = onCleanup(@() deleteIfPresent(temporaryPath));
save(temporaryPath, '-struct', 'payload');

savedVariables = string({whos('-file', temporaryPath).name});
expectedVariables = [string(variableNames); "README"];
assert(all(ismember(expectedVariables, savedVariables)), ...
    'saveDerivedFile:IncompleteArtifact', ...
    'The temporary MAT artifact is missing one or more required variables.');

movefile(temporaryPath, outputPath, 'f');
clear temporaryCleanup

end


function deleteIfPresent(filePath)
% Remove an incomplete temporary artifact after an error.

if isfile(filePath)
    delete(filePath);
end

end
