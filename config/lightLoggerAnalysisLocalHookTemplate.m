function lightLoggerAnalysisLocalHook

%  lightLoggerAnalysisLocalHook
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/lightLoggerAnalysisLocalHook.m
%
% Each time you run tbUseProject('lightLoggerAnalysis'), ToolboxToolbox will
% execute your local copy of this file to do setup.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


% Say hello.
projectName = 'lightLoggerAnalysis';

% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

% Get the DropBox path
if ismac
    dbJsonConfigFile = '~/.dropbox/info.json';
    fid = fopen(dbJsonConfigFile);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    dropboxBaseDir = val.business.path;
else
    error('Need to set up DropBox path finding for non-Mac machine')
end

% Set the prefs
setpref(projectName,'dropboxBaseDir',dropboxBaseDir); % main directory path 

% Set the default cal directory to the current project
calLocalData = fullfile(tbLocateProjectSilent(projectName),'cal');
setpref('combiLEDToolbox','CalDataFolder',calLocalData);

% Find the light logger directory 
light_logger_path = "";     

% Search recursively for folders named "lightLogger"
dirMatches = dir(fullfile('~/', '**', 'lightLogger'));

% Filter only directories (sometimes dir returns files too)
dirMatches = dirMatches([dirMatches.isdir]);

% Check if we found any matches
if ~isempty(dirMatches)
    % Take the first match (
        light_logger_path = fullfile(dirMatches(1).folder, dirMatches(1).name);

    % Save the path to a variable or a preference
    fprintf('Found lightLogger folder: %s\n', light_logger_path);
else
    error("ERROR: No lightLogger directory found"); 
end

% Save the LightLogger MATLAB libraries used in LightLogger Analysis
% as a pref
light_logger_libraries_matlab = fullfile(light_logger_path, "libraries_matlab"); 
setpref(projectName, 'light_logger_libraries_matlab', light_logger_libraries_matlab); 

% Save the paths to relevant Python libraries as prefs
Pi_util_path = fullfile(light_logger_path, "raspberry_pi_firmware", "utility", "Pi_uti.py");
setpref(projectName, 'Pi_util_path', Pi_util_path); 
world_util_path = fullfile(light_logger_path, "world", "world_util.py"); 
setpref(projectName, 'world_util_path', world_util_path); 
pupil_util_path = fullfile(light_logger_path, "pupil", "pupil_util.py"); 
setpref(projectName, 'pupil_util_path', pupil_util_path); 
ms_util_path = fullfile(light_logger_path, "ms", "ms_util.py"); 
setpref(projectName, 'ms_util_path', ms_util_path); 