function convert_recording_to_matlab(path_to_recording, output_path,...
                                     apply_digital_gain, use_mean_frame,...
                                     convert_time_units, convert_to_floats,...
                                     mean_axes,...
                                     contains_agc_metadata,...
                                     password...
                                    )
% Convert all of the individual chunks of a video to native MATLAB type and save them 
%
% Syntax:
%   convert_recording_to_matlab(path_to_recording, output_path,...
%                               apply_digital_gain, use_mean_frame,...
%                               convert_time_units, convert_to_floats,...
%                               mean_axes,...
%                               contains_agc_metadata,...
%                               password...
%                              )
%
% Description:
%   Given a path to a recording, read in all of the chunks from all of the 
%   sensors, apply desired transformations, and output the chunks 
%   converted to MATLAB type in output_path. 
%
% Inputs:
%   path_to_recording     - String. The path to the suprafile 
%                           containing all of the chunks 
%
%
%   apply_digital_gain    - Logical. Whether or not to apply 
%                           each frame's associated digital gain
%                           scalar
%
%   use_mean_frame        - Logical. Whether to only use the mean 
%                           pixel value from each frame for each 
%                           of the camera-based sensors.
%
%   convert_time_units    - Logical. Whether to convert the 
%                           timestamps from different sensors 
%                           all into the same units (s).
% 
%   convert_to_floats     - Logical. Whether to convert all 
%                           of the Python np.arrays to float types. 
%                           Can make things easier coming to MATLAB
%
%   mean_axes             - Struct. Represents the axes to mean 
%                           over when use_mean_frame is true. 
%                           Only applies to camera senors. Form
%                           is a struct with fields (WPM) 
%                           where each field is a tuple of axes. 
%                           0-indexed.
%
%   contains_agc_metdata  - Struct. Represents whether a sensor 
%                           has or does not have AGC metadata
%                           for a given recording. Form is 
%                           a struct with fields (WPM) where 
%                           each field is a boolean.
%
%   password               - String. Represents the password 
%                            used to encrypt the data (if encrypted)
%
% Outputs:
%
%   NONE 
%
% Examples:
%{
    path_to_experiment = '/Volumes/EXTERNAL1/test_folder_0';
    output_path = "./output_folder"; 
    apply_digital_gain = true; 
    use_mean_frame = false; 
    convert_time_units = true; 
    convert_to_floats = true; 
    convert_recording_to_matlab(path_to_recording, output_path, apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats, mean_axes, contains_agc_metadata, password)
%}                                    
    arguments 
        path_to_recording {mustBeText}; % Path to the recording file full of chunks 
        output_path {mustBeText}; % Path where the converted recording will be output
        apply_digital_gain {mustBeNumericOrLogical} = false; % Whether or not to apply digital gain values to each frame 
        use_mean_frame {mustBeNumericOrLogical} = false; % Whether to use the mean frame from the camera sensors or the entire frames 
        convert_time_units {mustBeNumericOrLogical} = false; % Whether to convert different time units from the different sensors all to seconds 
        convert_to_floats {mustBeNumericOrLogical} = false; % Whether to convert the Python np.arrays to float types. Can make things easier coming to MATLAB
        mean_axes = false; % The axes per sensor to apply mean over if we want to take some sort of mean. Note: ONLY for camera sensors  
        contains_agc_metadata = false; % Flags for each sensor if its metadata matrix contains AGC data
        password {mustBeText} = "1234"; % The password needed to decrypt encrypted + compressed files (.blosc files)
    end 
    
    % Begin timing the function 
    tic; 

    % Apply the default time ranges to splice out of the video (the entire video)
    % Default is a 1x2 None tuple. This signifies the entire video
    time_ranges = struct; 
    
    % Splice the entire video from all sensors 
    sensor_names = {'W', 'P', 'M'}; 
    for ss = 1:numel(sensor_names)  
        time_ranges.(sensor_names{ss}) = py.tuple({py.None, py.None}); 
    end 

    % Apply default mean axes if is not selected
    if(~isstruct(mean_axes))
        % Default is mean frame for each camera sensor 
        % Nothing for the MS
        mean_axes = struct; 

        % Define the default mean axes per sensor 
        sensor_names = {'W', 'P'}; 
        for ss = 1:numel(sensor_names)  
            mean_axes.(sensor_names{ss}) = py.tuple({1, 2}); 
        end 
        mean_axes.M = py.tuple();

    end 


    % Apply default contains AGC metadata if not supplied
    if(~isstruct(contains_agc_metadata))
        % Default is mean frame for each camera sensor 
        % Nothing for the MS
        contains_agc_metadata = struct; 

        % Define the default mean axes per sensor 
        contains_agc_metadata.W = true; 
        sensor_names = {'P', 'M'}; 
        for ss = 1:numel(sensor_names)  
            contains_agc_metadata.(sensor_names{ss}) = false; 
        end 
        
    end 



    % Because parsing chunks in MATLAB, for even a small number, is time consuming and memory 
    % intensive, and impossible for large videos, we will simply each single chunk 
    % at a time, then output this as a MATLAB compatible file so that one can simply read 
    % them in later 

    % Import the Pi_util python utility library 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    

    % First, let's find out the max number of chunks in the video for each sensor.
    % To do this, we will gather the chunk files per sensor 
    chunk_files_per_sensor = struct(Pi_util.group_sensors_files(path_to_recording));
    num_chunks_per_sensor = struct; 
    sensor_names = fieldnames(chunk_files_per_sensor); 
    for ss = 1:numel(sensor_names)
        % Retrieve the sensor name 
        sensor_name = sensor_names{ss};

        % Save how many chunks this sensor has
        num_chunks_per_sensor.(sensor_name) = length(chunk_files_per_sensor.(sensor_name)); 

    end     

    % Then, we will iterate over the max number of chunks
    max_num_chunks = max(cell2mat(struct2cell(num_chunks_per_sensor)));
    for ch = 1:max_num_chunks
        fprintf("Main | Processing chunk %d/%d...\n", ch, max_num_chunks); 

        % Initialize a struct to denote the chunk numbers to 
        % extract per sensor 
        chunk_ranges = struct; 
        
        % Set the range of chunks that we will extract 
        for ss = 1:numel(sensor_names)
            % Retrieve the sensor name 
            sensor_name = sensor_names{ss}; 

            % Construct the chunk ranges to splice 
            % (which will just be this current chunk)
            % for this sensor (if it has that many chunks)
            if(ch > num_chunks_per_sensor.(sensor_name))
                chunk_range = {py.None, py.None};
            else
                chunk_range = {ch-1, ch}; % Have to ensure it is exclusive Python list indexing
            end 

            % Save the chunk range for this sensor 
            chunk_ranges.(sensor_name) = py.tuple(chunk_range); 
        end 
        
        % Let's construct the argumetns to parse chunks 
        chunks_cell = parse_chunks(path_to_recording, ...
                                   apply_digital_gain, use_mean_frame,...
                                   convert_time_units, convert_to_floats,...
                                   time_ranges, chunk_ranges,...
                                   mean_axes,...
                                   contains_agc_metadata,...
                                   false,...
                                   password...
                                  ); 
        chunk = chunks_cell{ch}; % Fine to use MATLAB index here because we switch back 
        
        % Construct the output path for this chunk 
        output_filepath = fullfile(output_path, sprintf("chunk_%d.mat", ch-1)); 
        
        % Save the output chunk 
        save(output_filepath, "chunk", '-v7.3'); 

    end 
    
    % Calculate the time needed for conversion 
    elapsed_seconds = toc; 
    fprintf("Main | Conversion took %.3f seconds\n", elapsed_seconds); 

end