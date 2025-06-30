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
    %parsed_readings.temporal_sensitivity = convert_temporal_sensitivity_to_matlab(calibration_metadata.temporal_sensitivity, parsed_readings.temporal_sensitivity);
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
function converted_linearity = convert_ms_linearity_to_matlab(calibration_metadata, ms_linearity_readings)
    % First, let's extract some information about what we intended 
    num_NDF_levels = numel(calibration_metadata.NDFs);
    num_settings_levels = numel(calibration_metadata.background_scalars);
    n_measures = calibration_metadata.n_measures; 
    
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
                settings_idx = calibration_metadata.background_scalars_orders(nn, mm, ss); 
                setting = calibration_metadata.background_scalars(settings_idx);

                % Retrieve the measurement cell 
                measurement_cell = cell(settings_cell{mm}); 
                
                % Convert the chunks in the cell to MATLAB struct type 
                measurement_cell = cellfun(@(x) chunk_dict_to_matlab(x), measurement_cell, 'UniformOutput', false); 
                
                %save the updated measurement cell into the settings cell 
                converted_linearity{nn, ss, mm} = measurement_cell{:}; 
            end 

        end 

    end 

end 


% Local function to convert the temporal sensitivity and 
% phase fitting fields of the readings dict to pure 
% MATLAB types. Note: the CalibrationData referenced here is the substruct of CalibrationData
% specifically for temporal sensitivity
function converted_temporal_sensitivity = convert_temporal_sensitivity_to_matlab(CalibrationData, temporal_sensitivity)
    % Let's retrieve the original order of the contrast levels and frequencies 
    contrast_levels = CalibrationData.contrast_levels; 
    frequencies = CalibrationData.frequencies; 
    
    % Let's extract the order in which the contrast and frequencies were exposed 
    contrast_orders = CalibrationData.contrast_levels_orders; 
    frequencies_orders = CalibrationData.frequencies_orders; 
    
    % First, .temporal_sensitivity is a 1x1 py.list Convert this outer list to cell.
    % This gives you a cell array with 1 element, a YxZ py.list
    temporal_sensitivity = cell(temporal_sensitivity); 
    
    % Convert this py.list to a cell array 
    % This give you a cell array like { {YxZ} cell }, 
    % so let's splice out the extra dimension 
    temporal_sensitivity = cellfun(@(x) cell(x), temporal_sensitivity, 'UniformOutput', false); 
                                
    % This cell array is now a cell array where each element is a cell array of 
    % YxZ py.list. The outer cell is the contrast level and the inner cell 
    % is the frequency. First, let's put it back into matrix form. 
    temporal_sensitivity = vertcat(temporal_sensitivity{:});

    % First, let's calculate the actual size of the 
    % recording matrix and ensure it aligns with our expectations
    intended_size = [numel(contrast_levels), numel(frequencies), CalibrationData.n_measures];

    % If we intended to measure nothing, just return empty 
    if(intended_size(1) == 0)
        converted_temporal_sensitivity = {};
        return ; 
    end 

    [present_contrast_levels, present_frequencies] = size(temporal_sensitivity); 

    % Output a warning message here that there were 0 frequencies present 
    % and 0 contrast levels present. This can be intended behavior, 
    % at certain NDF levels where we do not record certain measurements 
    % TODO: Still need to fix this
    if(present_contrast_levels == 0 || present_frequencies == 0)
        converted_temporal_sensitivity = {};
        return ; 
    end 

    present_measures = length(temporal_sensitivity{1});
    actual_size = [present_contrast_levels, present_frequencies, present_measures];

    assert( isequal(intended_size, actual_size) );

    % Now, we get a Contrast x Frequency cell of 1xnMeasures cells. Let's create this into
    % a ContrastxFrequencyxMeasure cell 
    converted_temporal_sensitivity = cell(present_contrast_levels, present_frequencies, present_measures);

    % Now, we will easily iterate over/convert the measurements and store them in this 
    % new converted variable. Note: We will ALSO 
    % do the reordering/unshuffling component here 
    for cc = 1:present_contrast_levels
        for ff = 1:present_frequencies
            % First retrieve the measurement 
            % cell for this contrast and frequency 
            measurement_cell = temporal_sensitivity{cc, ff};
        
            % Every measurement is a py.list of 
            % reading dicts, so first convert this 
            % measurement into cell, then 
            % we will parse the reading dicts to MATLAB 
            % types 
            for nn = 1:present_measures
                video_chunks_as_cell = cell(measurement_cell{nn}); 
                video_chunks_as_cell = cellfun(@(x) chunk_dict_to_matlab(x), video_chunks_as_cell); 

                % Now, let's find out what REAL idx was this combination of measurement number, 
                % contrast number, and frequency number 
                contrast_idx = contrast_orders(nn, cc); 
                frequency_idx = frequencies_orders(nn, cc, ff); 

                % Store the fully converted reading back into the larger cell 
                converted_temporal_sensitivity{contrast_idx, frequency_idx, nn} = video_chunks_as_cell;
                
            end 

        end 

    end 

    return ; 

end