function convert_recording_to_matlab(path_to_recording, output_path,...
                                     apply_digital_gain, use_mean_frame,...
                                     convert_time_units, convert_to_floats,...
                                     mean_axes,...
                                     contains_agc_metadata,...
                                     password...
                                    )
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
    
    % Because parsing chunks in MATLAB, for even a small number, is time consuming and memory 
    % intensive, and impossible for large videos, we will simply each single chunk 
    % at a time, then output this as a MATLAB compatible file so that one can simply read 
    % them in later 
    
    % First, let's find out the max number of chunks in the video for each sensor
    num_chunks_per_sensor.W = 3; 
    num_chunks_per_sensor.P = 2; 
    num_chunks_per_sensor.M = 1;  % Dummy data 
    

    % Then, we will iterate over the max number of chunks
    max_num_chunks = max([num_chunks_per_sensor.W, num_chunks_per_sensor.P, num_chunks_per_sensor.M]);  
    for ch = 1:max_num_chunks
        % Construct the chunk ranges to splice 
        chunk_ranges.W = 1; 
        chunk_ranges.P = 1; 
        chunk_ranges.M = 1; 
        
        % Let's construct the argumetns to parse chunks 
        chunk_cells = parse_chunks(path_to_recording, ...
                                   apply_digital_gain, use_mean_frame,...
                                   convert_time_units, convert_to_floats,...
                                   false, chunk_ranges,...
                                   mean_axes,...
                                   contains_agc_metadata,...
                                   password...
                                  ); 
        chunk = chunks_cell{1};      
        
        % Construct the output path for this chunk 
        output_filepath = fullfile(output_path, sprintf("chunk_%d.mat", ch)); 
        
        % Save the output chunk 
        save(output_filepath, "chunk"); 

    end 
    




end