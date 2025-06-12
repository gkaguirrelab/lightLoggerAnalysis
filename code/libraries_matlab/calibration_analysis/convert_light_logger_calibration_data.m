function LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder, NDF, use_mean_frame, mean_axes)
    arguments 
        experiment_folder {mustBeText}; 
        NDF {mustBeInteger}; 
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

    % Include dependencies (without this the rawData struct of the cal field)
    % of the metadata struct will be empty 
    tbUse('combiLEDToolbox'); 

    % Add the path to our misc MATLAB utility functions
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab"));

    % Import the Python library used for downloading from Dropbox 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path"));

    % Load in the Calibration Metadata 
    NDF_folder_path = fullfile(experiment_folder, sprintf("NDF_%d", NDF)); 
    calibration_metadata = load_calibration_metadata(experiment_folder, NDF_folder_path); 
   
    % Next, we will load in the recordings and rudimentarily convert the Py.dict to a struct. 
    % Further nested conversion will be needed for each reading method
    disp("Conversion | Parsing recordings...")
    parsed_readings = struct(Pi_util.load_sorted_calibration_files(NDF_folder_path, true, use_mean_frame, true, true, mean_axes)); 

    % Convert the subfields to purely MATLAB type
    parsed_readings.ms_linearity = convert_ms_linearity_to_matlab(calibration_metadata.ms_linearity, parsed_readings.ms_linearity);
    parsed_readings.temporal_sensitivity = convert_temporal_sensitivity_to_matlab(calibration_metadata.temporal_sensitivity, parsed_readings.temporal_sensitivity);
    parsed_readings.phase_fitting = convert_temporal_sensitivity_to_matlab(calibration_metadata.phase_fitting, parsed_readings.phase_fitting);
    parsed_readings.contrast_gamma = convert_temporal_sensitivity_to_matlab(calibration_metadata.contrast_gamma, parsed_readings.contrast_gamma);

    % Initialize a return struct 
    LightLoggerCalibrationData.metadata = calibration_metadata;
    LightLoggerCalibrationData.readings = parsed_readings; 

    return ; 

end 

% Local function to load in the Calibration struct from the raw bytes 
function calibration_metadata = load_calibration_metadata(experiment_folder, NDF_folder_path)
    % Next, we will construct the path to the CalibrationData struct which contains information about 
    % the experiment
    calibration_metadata_struct_filepath = fullfile(NDF_folder_path, "CalibrationData.mat.bytes"); 
    calibration_metadata_struct_file = fopen(calibration_metadata_struct_filepath, 'rb'); % Open the file whose bytes represent the CalibrationData struct 
    calibration_metadata_struct_bytes = fread(calibration_metadata_struct_file, Inf, '*uint8'); % Load in the bytes from that file
    calibration_metadata = getArrayFromByteStream(calibration_metadata_struct_bytes); % Parse the bytes back into a struct
    fclose(calibration_metadata_struct_file); % Close the file
end 

% Local function to convert the ms linearity field of the readings dict to pure 
% MATLAB types. Note: the CalibrationData referenced here is the substruct of CalibrationData
% specifically for ms_linearity
function converted_linearity = convert_ms_linearity_to_matlab(CalibrationData, ms_linearity)
    % First, let's extract some information about what we intended 
    num_settings_levels = numel(CalibrationData.background_scalars);
    n_measures = CalibrationData.n_measures; 
    
    % .ms_linearity is a YxZ py.list. Convert this outer list to cell.
    % This gives you a cell of py.lists. The outer cell represents the 
    % settings number and the inner one represents the measurements 
    % for this settings level 
    ms_linearity = cell(ms_linearity); 

    % Convert those inner py.list to cell. Now you have a cell of cells 
    ms_linearity = cellfun(@(x) cell(x), ms_linearity, 'UniformOutput', false); 

    % Form this back into a cell matrix of shape settingsLevel x measurementNum 
    ms_linearity = vertcat(ms_linearity{:});

    % Let's do some checking that these dimensions match what we intended 
    intended_size = [num_settings_levels, n_measures]; 
    actual_size = size(ms_linearity); 
    assert(isequal(intended_size, actual_size)); 

    % Convert the inner elements, 1x1 py.list of recording chunks 
    % to cells (thus making them 1x1 py.dict)
    ms_linearity = cellfun(@(x) cell(x), ms_linearity); 

    % Now that we have it in the correct shape, 
    % let's create a cell for the reorganized values 
    converted_linearity = cell(num_settings_levels, n_measures);

    % Iterate over the settings level 
    for ss = 1:num_settings_levels
        % Iterate over the measurements per settings level
        for nn =1:n_measures
            % Find the index of the setting that was used for this combination 
            % of measurement and settings number 
            settings_idx = CalibrationData.background_scalars_orders(nn, ss); 
            
            % Convert this dictionary to MATLAB types and save it into the matrix 
            converted_linearity{settings_idx, nn} = chunk_dict_to_matlab(ms_linearity{ss, nn}); 
        end 
    end 

    return ; 

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