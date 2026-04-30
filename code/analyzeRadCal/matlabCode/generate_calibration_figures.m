function generate_calibration_figures(output_path)
    arguments 
        output_path; 
    end     

    % Get the path to the DropBox location where calibration data lives
    dropbox_basedir = getpref('lightLoggerAnalysis', 'dropboxBaseDir'); 

    % Get the path to the radiometricx calibration 
    radiometric_calibration_path = fullfile(dropbox_basedir, 'FLIC_data', 'lightLoggerRadCal', 'W1P1M1'); 
    assert(exist(radiometric_calibration_path, 'dir')); 

    % Next, iterate over the files in this directory
    radiometric_calibration_experiments = dir(radiometric_calibration_path)

    % Iterate over the files 
    for experiment_num = 1:numel(radiometric_calibration_experiments)
        % Get the experiment itself
        experiment = radiometric_calibration_experiments(experiment_num);

        % Skip hidden files
        if(startsWith(experiment.name, "."))
            continue; 
        end 

        % Skip the test file 
        if(experiment.name == "test")
            continue; 
        end 

        % Gather the path to this expeirment 
        experiment_path = fullfile(radiometric_calibration_path, experiment.name); 
        
        % Find the mat file that is in the experiments folder
        experiment_files = dir(experiment_path);
        experiment_results_path = "";
        for file_num = 1:numel(experiment_files)
            file = experiment_files(file_num);

            % Skip hidden files 
            if(startsWith(file.name, '.'))
                continue; 
            end     
            
            if(ismember('.mat', file.name) )
                experiment_results_path = fullfile(experiment_path, file.name); 
                break; 
            end 

        end 
        assert(experiment_results_path ~= "")

        % Load in the calibration data for this experiment
        light_logger_calibration_data = load(experiment_results_path); 
        experiment_fieldnames = fieldnames(light_logger_calibration_data); 
        experiment_fieldname = experiment_fieldnames{1}; 
        light_logger_calibration_data = light_logger_calibration_data.(experiment_fieldname); 

        % Generate the figures 
        analyze_light_logger_calibration_data(light_logger_calibration_data, 'figure_output_path', output_path); 

    end 



end