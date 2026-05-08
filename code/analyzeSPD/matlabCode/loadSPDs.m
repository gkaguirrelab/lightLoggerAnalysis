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
%   containing previously generated SPD outputs, filters those outputs by
%   subject, activity, color mode, and projection type, and returns a
%   fully loaded MATLAB structure.
%
%   Internally, the Python helper first builds a nested structure of
%   filepaths to the saved `.mat` files. This MATLAB function then
%   recursively walks that filepath structure and replaces each terminal
%   filepath leaf with the loaded `.mat` contents. The returned structure
%   therefore mirrors the on-disk SPD organization while containing loaded
%   MATLAB data rather than path strings.
%
% Inputs:
%   src_dir                 - Char/string. Path to the root directory that
%                             contains the saved SPD outputs organized by
%                             color mode, subject, and activity.
%
% Optional key/value pairs:
%   subjects_to_skip        - Cell array / string array. Subject folder
%                             names to exclude from loading.
%   subjects_to_process     - Cell array / string array. Subject folder
%                             names to include. If non-empty, only these
%                             subjects are considered.
%   activities_to_skip      - Cell array / string array. Activity folder
%                             names to exclude from loading.
%   activities_to_process   - Cell array / string array. Activity folder
%                             names to include. If non-empty, only these
%                             activities are considered.
%   color_modes_to_skip     - Cell array / string array. Color-mode folder
%                             names to exclude from loading.
%   color_modes_to_process  - Cell array / string array. Color-mode folder
%                             names to include. If non-empty, only these
%                             color modes are considered.
%   projection_types_to_skip
%                           - Cell array / string array. Projection types
%                             to exclude, e.g. {"justProjection"}.
%   projection_types_to_process
%                           - Cell array / string array. Projection types
%                             to include. If non-empty, only these
%                             projection types are considered.
%
% Outputs:
%   spds                    - Struct. Nested structure containing the
%                             loaded SPD results. The structure shape is:
%
%                               spds.(color_mode).(subject_id).(activity_name).(projection_type)
%
%                             where each terminal `(projection_type)` node
%                             contains:
%
%                               .spd
%                               .best_fit
%
%                             and each of those fields contains the loaded
%                             contents of the corresponding `.mat` file.
%
%                             So the full output hierarchy is:
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
%       "projection_types_to_process", {"virtuallyFoveated", "justProjection"});
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
    end

    % Load the Python utility module that knows how to discover
    % SPD output files on disk.
    spd_util = import_pyfile(getpref("lightLoggerAnalysis", "spd_util_path")); 

    % Create a unique temporary `.mat` filepath that the Python helper
    % can write into. We use a UUID so concurrent runs do not collide.
    temp_output_path = sprintf("temp_loadSPDs_%s.mat", char(java.util.UUID.randomUUID)); 

    % Ask Python to discover the requested SPD files.
    % `paths_only = true` means the Python side should return a nested
    % structure of filepaths rather than loading the contents itself.
    % `output_as_mat = temp_output_path` tells Python to serialize that
    % nested path structure into a temporary `.mat` file for MATLAB.
    spd_util.load_spds(src_dir, ...
                    pyargs("subjects_to_skip", options.subjects_to_skip, ...
                    "subjects_to_process", options.subjects_to_process, ...
                    "activities_to_skip", options.activities_to_skip, ...
                    "activities_to_process", options.activities_to_process, ...
                    "color_modes_to_skip", options.color_modes_to_skip, ...
                    "color_modes_to_process", options.color_modes_to_process, ...
                    "projection_types_to_skip", options.projection_types_to_skip, ...
                    "projection_types_to_process", options.projection_types_to_process, ...
                    "paths_only", true, ...
                    "output_as_mat", temp_output_path)...
                ); 

    % Run Python garbage collection after the helper call in case the
    % helper created many temporary Python objects.
    py.gc.collect(); 

    % Load the temporary `.mat` file. At this point, `spd_paths` is still
    % only a nested structure of filepaths to the actual SPD result files.
    spd_paths = load(temp_output_path).spds; 

    % Recursively walk the filepath structure and replace each terminal
    % filepath leaf with the loaded `.mat` contents from that filepath.
    spds = iLoadSpdPathStruct(spd_paths); 

    % Remove the temporary `.mat` file now that we have loaded it into
    % native MATLAB memory.
    delete(temp_output_path); 

    return 


end 


function output_struct = iLoadSpdPathStruct(input_struct)
% Recursively walk the nested SPD filepath struct and load terminal
% `.spd` / `.best_fit` mat-file leaves into MATLAB structs.

    % Base case: if the current input is not a struct, there is nothing
    % further to recurse into, so return the value unchanged.
    if (~isstruct(input_struct))
        output_struct = input_struct;
        return;
    end

    % Retrieve the fieldnames at the current node so we can determine
    % whether this is an internal node or a terminal projection leaf.
    input_fields = fieldnames(input_struct);

    % If this node contains `.spd` or `.best_fit`, treat it as a terminal
    % projection-type node whose values are filepaths to `.mat` files.
    if (ismember("spd", input_fields) || ismember("best_fit", input_fields))

        % Initialize the output struct for this terminal node.
        output_struct = struct();

        % If an SPD filepath exists at this node, load that `.mat` file
        % and store its contents under `.spd`.
        if (ismember("spd", input_fields))
            output_struct.spd = iLoadMatFileFromPath(input_struct.spd);
        end

        % If a best-fit filepath exists at this node, load that `.mat`
        % file and store its contents under `.best_fit`.
        if (ismember("best_fit", input_fields))
            output_struct.best_fit = iLoadMatFileFromPath(input_struct.best_fit);
        end

        % Once we have loaded the terminal files, stop recursing.
        return;
    end

    % Otherwise, this is an internal node such as color mode, subject,
    % activity, or projection-type grouping. Create an output struct and
    % recurse into each child field one at a time.
    output_struct = struct();
    for ff = 1:numel(input_fields)

        % Retrieve the current child field name.
        field_name = input_fields{ff};

        % Recurse into that child and store the fully loaded result back
        % into the matching field of the output struct.
        output_struct.(field_name) = iLoadSpdPathStruct(input_struct.(field_name));
    end
end


function loaded_value = iLoadMatFileFromPath(path_value)
% Load a .mat file from a filepath leaf stored in the intermediate SPD
% path structure.

    % Accept MATLAB strings directly.
    if (isstring(path_value))
        mat_path = path_value;

    % Accept character vectors as well, then normalize to string.
    elseif (ischar(path_value))
        mat_path = string(path_value);

    % Anything else means the intermediate path structure was not in the
    % expected filepath format.
    else
        error("Expected terminal SPD path leaf to be char or string, got %s", class(path_value));
    end

    % Load and return the `.mat` file contents at this filepath.
    loaded_value = load(mat_path);
end
