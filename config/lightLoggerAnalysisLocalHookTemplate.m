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

% Set the light logger recordings directory
setpref(projectName,'dataBaseDir',fullfile(dropboxBaseDir,'FLIC_data'));

% Set the default cal directory to this project
calLocalData = fullfile(tbLocateProject(projectName),'cal');
setpref('combiLEDToolbox','CalDataFolder',calLocalData);

% Find the light logger directory 
light_logger_path = tbLocateProject('lightLogger');
combiExperiments_path = tbLocateProject('combiExperiments');
     
% Save the LightLogger MATLAB libraries used in LightLogger Analysis
% as a pref
light_logger_libraries_matlab = fullfile(light_logger_path, "libraries_matlab"); 
setpref(projectName, 'light_logger_libraries_matlab', light_logger_libraries_matlab); 
addpath(light_logger_libraries_matlab);

% Save the path to combiExperiments 
setpref(projectName, 'combiExperiments_path', combiExperiments_path);

% Save the paths to relevant Python libraries as prefs
Pi_util_path = fullfile(light_logger_path, "raspberry_pi_firmware", "utility", "Pi_util.py");
setpref(projectName, 'Pi_util_path', Pi_util_path); 
world_util_path = fullfile(light_logger_path, "world", "world_util.py"); 
setpref(projectName, 'world_util_path', world_util_path); 
pupil_util_path = fullfile(light_logger_path, "pupil", "pupil_util.py"); 
setpref(projectName, 'pupil_util_path', pupil_util_path); 
ms_util_path = fullfile(light_logger_path, "ms", "ms_util.py"); 
setpref(projectName, 'ms_util_path', ms_util_path); 

% Add all functions in the current directory to the path 
light_logger_analysis_path = tbLocateProject('lightLoggerAnalysis');
addpath(genpath(light_logger_analysis_path));


% Next, we will save some important constants that we frequently use 
world_util = import_pyfile(world_util_path);
world_fps = double(world_util.WORLD_CAM_FPS); 
setpref(projectName, "world_fps", world_fps);
clear world_util; 

% For these constants, we measured the modulation depth of the stimulus source as a function of
% contrast using the Klein photometer and the "ChromaSurf" sofwate. The
% values are given here, and used for correction of the observed response
% amplitude below
test_mod_depth_frequencies = [0.25,0.5,1,3,6,12,25,50,100];
test_mod_depth_amplitudes = [1,1,1,1,0.983,0.975,0.939,0.90,0.66];
setpref(projectName, "test_mod_depth_frequencies", test_mod_depth_frequencies); 
setpref(projectName, "test_mod_depth_amplitudes", test_mod_depth_amplitudes);
