    function LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder,...
                                                                                apply_digital_gain,...
                                                                                convert_time_units,...
                                                                                convert_to_floats,...
                                                                                use_mean_frame,...
                                                                                mean_axes...
                                                                            )
    arguments 
        experiment_folder {mustBeText}; 
        apply_digital_gain = false;  
        convert_time_units = false; 
        convert_to_floats = false; 
        use_mean_frame {mustBeNumericOrLogical} = false; 
        mean_axes = false; 
    end     

    % If we have not been passed in specific axis to mean, use the default
    % which is to mean each frame
    if(~isstruct(mean_axes))
        mean_axes = struct; 
        mean_axes.W = [1, 2]; 
        mean_axes.P = [1, 2];
        mean_axes.M = [];
    end 

    disp("Conversion | Importing libraries...");
    tbUseProject('lightLoggerAnalysis');

    % Import the Python library used for downloading from Dropbox 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path"));

    % Next, we will load in the recordings and rudimentarily convert the Py.dict to a struct. 
    % Further nested conversion will be needed for each reading method
    disp("Conversion | Parsing recordings...")
    parsed_readings = struct(Pi_util.load_sorted_calibration_files(experiment_folder,...
                                                                   apply_digital_gain,...
                                                                   use_mean_frame,...
                                                                   convert_time_units,...
                                                                   convert_to_floats,...
                                                                   mean_axes...
                                                                  )...
                            ); 

    % Save the path to the current file's directory. we are going to need this when we mess around with tbUses
    path_to_file_dir = fileparts(mfilename("fullpath"));

    % Load in the Calibration Metadata 
    % Include dependencies (without this the rawData structs of the cal field)
    % of the metadata struct will be empty 
    tbUse('combiLEDToolbox'); 
    addpath(path_to_file_dir);
    calibration_metadata = load_calibration_metadata(fullfile(experiment_folder, "calibration_metadata.mat.bytes")); 
    tbUseProject('lightLoggerAnalysis');

    % Convert the subfields to purely MATLAB type
    parsed_readings.ms_linearity = convert_ms_linearity_to_matlab(calibration_metadata.ms_linearity, parsed_readings.ms_linearity);
    parsed_readings.temporal_sensitivity = convert_temporal_sensitivity_to_matlab(calibration_metadata.temporal_sensitivity, parsed_readings.temporal_sensitivity);
    %parsed_readings.phase_fitting = convert_temporal_sensitivity_to_matlab(calibration_metadata.phase_fitting, parsed_readings.phase_fitting);
    %parsed_readings.contrast_gamma = convert_temporal_sensitivity_to_matlab(calibration_metadata.contrast_gamma, parsed_readings.contrast_gamma);

    % Initialize a return struct 
    LightLoggerCalibrationData.metadata = calibration_metadata;
    LightLoggerCalibrationData.readings = parsed_readings; 

    return ; 

end 

% Local function to convert the ms linearity field of the readings dict to pure 
% MATLAB types. Note: the CalibrationData referenced here is the substruct of CalibrationData
% specifically for ms_linearity
function converted_linearity = convert_ms_linearity_to_matlab(linearity_calibration_metadata, ms_linearity_readings)
    % First, let's extract some information about what we intended 
    num_NDF_levels = numel(linearity_calibration_metadata.NDFs);
    num_settings_levels = numel(linearity_calibration_metadata.background_scalars);
    n_measures = linearity_calibration_metadata.n_measures; 
    
    % Now, we will allocate an output array 
    converted_linearity = cell(num_NDF_levels, num_settings_levels, n_measures);

    % Then, we will convert the outer Python list, representing the NDF levels, 
    % to cell array 
    ms_linearity_readings = cell(ms_linearity_readings); 

    % Next, we will convert each inner 1xY py.list to a cell array. 
    % These lists represent the lists of recordings at different settings levels 
    ms_linearity_readings = cellfun(@(x) cell(x), ms_linearity_readings, 'UniformOutput', false);

    % Next, we need to go in all of those recording cells and convert the inner py.list of 
    % measurements at a given settings level to cell 
    for nn = 1:num_NDF_levels
        % Convert all the repeated measurements at this settings level to cell 
        NDF_cell = cellfun(@(x) cell(x), ms_linearity_readings{nn}, 'UniformOutput', false);   

        % Then, convert those 1x1 py.lists representing chunks of recordings for each measurement 
        % to cell 
        for ss = 1:num_settings_levels
            % Retrieve the settings level cell 
            settings_cell = cell(NDF_cell{ss}); 
            
            % Iterate over the measures at this settings level 
            for mm = 1:n_measures
                % Find the index of the setting that was used for this combination 
                % of measurement and settings number 
                settings_idx = linearity_calibration_metadata.background_scalars_orders(nn, mm, ss); 
                setting = linearity_calibration_metadata.background_scalars(settings_idx);

                % Retrieve the measurement cell 
                measurement_cell = cell(settings_cell{mm}); 
                
                % Convert the chunks in the cell to MATLAB struct type 
                measurement_cell = cellfun(@(x) chunk_dict_to_matlab(x), measurement_cell, 'UniformOutput', false); 
                
                % Save the updated measurement cell into the settings cell 
                converted_linearity{nn, settings_idx, mm} = measurement_cell{:}; 
            end 

        end 

    end 

end 


% Local function to convert the temporal sensitivity and 
% phase fitting fields of the readings dict to pure 
% MATLAB types. Note: the CalibrationData referenced here is the substruct of CalibrationData
% specifically for temporal sensitivity
function converted_temporal_sensitivity = convert_temporal_sensitivity_to_matlab(temporal_sensitivity_calibration_metadata, temporal_sensitivity_readings)
    % First, let's retrieve some basic information from the struct 
    num_NDF_levels = numel(temporal_sensitivity_calibration_metadata.NDFs); 
    num_contrast_levels = numel(temporal_sensitivity_calibration_metadata.contrast_levels); 
    num_frequencies = numel(temporal_sensitivity_calibration_metadata.frequencies); 
    n_measures = temporal_sensitivity_calibration_metadata.n_measures; 

    % Allocate a converted output 
    converted_temporal_sensitivity = cell(num_NDF_levels, num_contrast_levels, num_frequencies, n_measures); 

    % Then Let's retrieve the original order of the contrast levels and frequencies 
    contrast_levels = temporal_sensitivity_calibration_metadata.contrast_levels; 
    frequencies = temporal_sensitivity_calibration_metadata.frequencies; 
    
    % Then Let's extract the order in which the contrast and frequencies were exposed 
    contrast_orders = temporal_sensitivity_calibration_metadata.contrast_levels_orders; 
    frequencies_orders = temporal_sensitivity_calibration_metadata.frequencies_orders; 

    % First, .temporal_sensitivity is a 1x1 py.list Convert this outer list to cell.
    % This gives you a cell array with of N py.lists corresponding to the NDF levels 
    temporal_sensitivity_readings = cell(temporal_sensitivity_readings); 

    % Then, convert each NDF level measurement py.list into a cell array as well 
    temporal_sensitivity_readings = cellfun(@(x) cell(x), temporal_sensitivity_readings, 'UniformOutput', false); 

    % Now that we've done some basic conversion, let's do some iteration to fill in our converted cell 
    % First, iterate over the NDF levels 
    for nn = 1:num_NDF_levels
        % Let's retrieve the contrasts measured at this NDF
        contrast_level_readings = cellfun(@(x) cell(x), temporal_sensitivity_readings{nn}, 'UniformOutput', false);
        
        % Next, iterate over the contrast levels 
        for cc = 1:num_contrast_levels
            % Retrieve the cell of frequencies measured at this contrast 
            % Convert those py.list of measures per frequency to cell 
            frequencies_readings = cellfun(@(x) cell(x), contrast_level_readings{cc}, 'UniformOutput', false); 
            
            % Then, iterate over the frequencies 
            for ff = 1:num_frequencies
                % Retrieve the measurements at this frequnecy 
                frequency_readings = cellfun(@(x) cell(x), frequencies_readings{ff}, 'UniformOutput', false);
        
                % Then, iterate over the measure at this frequency 
                for mm = 1:n_measures
                    % Retrieve the cell array of py.dict chunks that compose the measurement
                    measurement_chunks = frequency_readings{mm};

                    % Convert the chunks to MATLAB 
                    measurement_chunks_matlab = cellfun(@(x) chunk_dict_to_matlab(x), measurement_chunks, 'UniformOutput', false);

                    % Retrieve the mmth measurement at this frequency 
                    % convert the chunks of this measurement to MATLAB type, and save in struct 
                    converted_temporal_sensitivity{nn, cc, ff, mm} = measurement_chunks_matlab{:};
                end 
            
            end 

        end 

    end 

    return ; 

end