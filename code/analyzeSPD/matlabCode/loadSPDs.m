function spds = loadSPDS(src_dir, options)
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
