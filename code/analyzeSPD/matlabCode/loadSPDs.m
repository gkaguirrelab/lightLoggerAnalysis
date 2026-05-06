function spds = loadSPDS(src_dir, options)
% Load saved SPD result files into a nested MATLAB structure
%
% Syntax:
%   spds = loadSPDS(src_dir)
%   spds = loadSPDS(src_dir, options)
%
% Description:
%   This function is a MATLAB wrapper around the Python helper
%   `load_spds` defined in `spd_util.py`. It searches a directory tree
%   containing previously generated SPD analysis outputs, filters those
%   outputs according to the requested subjects, activities, color modes,
%   and projection types, and then returns the discovered SPD filepaths
%   and associated loaded content as a nested MATLAB structure.
%
%   The wrapper can optionally request that best-fit-line results be
%   included alongside each SPD entry. Internally, the Python helper can
%   serialize the discovered structure to a temporary `.mat` file, which
%   this MATLAB function then loads and returns as a native MATLAB struct.
%
% Inputs:
%   src_dir               - Char/string. Path to the root directory that
%                           contains the saved SPD outputs organized by
%                           color mode, subject, and activity.
%
% Optional key/value pairs:
%   subjects_to_skip      - Cell array / string array. Subject folder
%                           names to exclude from loading.
%   subjects_to_process   - Cell array / string array. Subject folder
%                           names to load. If non-empty, only these
%                           subjects are considered.
%   activities_to_skip    - Cell array / string array. Activity folder
%                           names to exclude from loading.
%   activities_to_process - Cell array / string array. Activity folder
%                           names to load. If non-empty, only these
%                           activities are considered.
%   color_modes_to_skip   - Cell array / string array. Color-mode folder
%                           names to exclude from loading.
%   color_modes_to_process
%                         - Cell array / string array. Color-mode folder
%                           names to load. If non-empty, only these color
%                           modes are considered.
%   projection_types_to_skip
%                         - Cell array / string array. Projection types
%                           to exclude, e.g. {"virtuallyFoveated"}.
%   projection_types_to_process
%                         - Cell array / string array. Projection types
%                           to include. If non-empty, only these
%                           projection types are considered.
%   include_best_fit      - Logical. If true, request that best-fit-line
%                           outputs also be loaded for each SPD entry.
%   output_as_mat         - Char/string. Optional path to a `.mat` file
%                           where the Python helper should write the
%                           loaded SPD structure. If empty, a temporary
%                           `.mat` filepath is created automatically.
%
% Outputs:
%   spds                  - Struct. Nested structure containing the loaded
%                           SPD information. The structure shape is:
%
%                               spds.(color_mode).(subject_id).(activity_name).(projection_type)
%
%                           where each terminal `(projection_type)` field
%                           contains a structure with:
%
%                               .spd
%                               .best_fit
%
%                           In other words, the full output hierarchy is:
%
%                               spds
%                                 -> color_mode
%                                 -> subject_id
%                                 -> activity_name
%                                 -> projection_type
%                                 -> spd / best_fit
%
% Examples:
%{
%   src_dir = "/path/to/combined_spd_outputs";
%
%   spds = loadSPDS(src_dir, ...
%       "subjects_to_process", {"FLIC_2002"}, ...
%       "activities_to_process", {"walkIndoor", "walkOutdoor"}, ...
%       "color_modes_to_process", {"a", "c_lm", "c_s"}, ...
%       "projection_types_to_process", {"virtuallyFoveated", "justProjection"}, ...
%       "include_best_fit", true);
%}
    arguments 
        src_dir
        options.subjects_to_skip = {}
        options.subjects_to_process = {}
        options.activities_to_skip = {}
        options.activities_to_process = {}
        options.color_modes_to_skip = {}
        options.color_modes_to_process = {}
        options.projection_types_to_skip = {}
        options.projection_types_to_process = {}
        options.include_best_fit = false
    end

    % Load the python utility 
    spd_util = import_pyfile(getpref("lightLoggerAnalysis", "spd_util_path")); 

    % Generate the structure with the desired SPDs loaded 
    temp_output_path = sprintf("temp_loadSPDs_%s.mat", char(java.util.UUID.randomUUID)); 

    spd_util.load_spds(src_dir, ...
                       pyargs("subjects_to_skip", options.subjects_to_skip, ...
                              "subjects_to_process", options.subjects_to_process, ...
                              "activities_to_skip", options.activities_to_skip, ...
                              "activities_to_process", options.activities_to_process, ...
                              "color_modes_to_skip", options.color_modes_to_skip, ...
                              "color_modes_to_process", options.color_modes_to_process, ...
                              "projection_types_to_skip", options.projection_types_to_skip, ...
                              "projection_types_to_process", options.projection_types_to_process, ...
                              "include_best_fit", options.include_best_fit, ...
                              "output_as_mat", temp_output_path)...
                        ); 
    py.gc.collect(); 

    % Load the structure in and return
    spds = load(temp_output_path).spds; 

    % Remove the temp output path
    delete(temp_output_path); 

    return 


end 
