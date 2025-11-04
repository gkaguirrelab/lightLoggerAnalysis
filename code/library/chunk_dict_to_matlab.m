function converted = chunk_dict_to_matlab(py_chunk_dict)
% Convert a py.dict representing a single chunk of recording data into purely MATLAB type
%
% Syntax:
%   chunk_dict_to_matlab(py_chunk_dict)
%
% Description:
%   Given a py.dict representing a single chunk of recording from 
%   the light logger, with fields W, P, M, S, corresponding to 
%   each of the sensors, convert all of the data from the sensors
%   to purely MATLAB type
%
% Inputs:
%   py_chunk_dict    - py.dict. Python dictionary representing 
%                      all the sensor data for a single chunk 
% Outputs:
%
%   converted        - Struct. Struct equivalent of the input 
%                      dictionary with all data converted to MATLAB 
%                      native types. 
%
% Examples:
%{
    % Load in the Pi util hepler file 
    Pi_util = import_pyfile(Pi_util_path);

    % Parse the experiment using Python and return as py.list
    chunks_as_py = Pi_util.parse_chunks(path_to_experiment,...
                                        apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats,...
                                        time_ranges, chunk_ranges, mean_axes, contains_agc_metadata,...
                                        password...
                                       );
    
    % Convert outer Python list to cell 
    chunks = cell(chunks_as_py);  
    
    % Convert the first chunk to MATLAB native type 
    MATLAB_chunk = chunk_dict_to_matlab(chunks{1})
%}


    % Convert the py.dict containing dicts of all sensors' data 
    % to struct
    chunk_as_struct =  struct(py_chunk_dict); 

    % Iterate over the sensors' dicts and convert them to MATLAB type 
    field_names = fieldnames(chunk_as_struct);
    for ss = 1:numel(field_names)
        % Retrieve this sensor's data dict as a struct
        sensor_struct = struct(chunk_as_struct.(field_names{ss})); 

        % Convert the temporal and value fields to MATLAB double type 
        % Treat the MS slightly differently as it has subfields per sensor 
        if(field_names{ss} == 'M')
            % Retrieve the t struct 
            t_struct = struct(sensor_struct.('t')); 

            % Next, we will iterate over the sensor names and convert their t vectors 
            sub_field_names = fieldnames(t_struct);
            for sf = 1:numel(sub_field_names)
                t_struct.(sub_field_names{sf}) = double(t_struct.(sub_field_names{sf})); 
            end

            % Save the convert t struct 
            sensor_struct.('t') = t_struct; 

        else
            % Simply convert the t numpy array to double  
            sensor_struct.('t') = double(sensor_struct.('t'));
        end 
        
        % Treat the MS field slightly differently as its values are further 
        % broken down into another dictionary per sensor
        if(field_names{ss} == 'M')
            % Retrieve the values dictionary and convert to 
            % struct right away
            data_frames_struct = struct(sensor_struct.('v'));

            % Next, we will go over and convert all the numpy array to MATLAB array
            sub_field_names = fieldnames(data_frames_struct);
            for sf = 1:numel(sub_field_names)
                data_frames_struct.(sub_field_names{sf}) = double(data_frames_struct.(sub_field_names{sf})); 
            end
            
            % Save the converted subdict 
            sensor_struct.('v') = data_frames_struct;
        
        % Otherwise, simply convert the value field to double 
        else 
            sensor_struct.('v') = double(sensor_struct.('v'));  
        end
        
        % If we are converting the world camera, we must also convert its sensor field 
        if(field_names{ss} == "W" || field_names{ss} == "P") 
            % Retrieve the settings field of the struct 
            settings_struct = struct(sensor_struct.('settings')); 

            % Now iterate over the settings fields and convert them to doubles 
            sub_field_names = fieldnames(settings_struct); 
            for sf = 1:numel(sub_field_names)
                settings_struct.(sub_field_names{sf}) = double(settings_struct.(sub_field_names{sf})); 
            end 
            
            % Save the converted subdict 
            sensor_struct.('settings') = settings_struct; 

        end 

        % Replace the Python data type with the MATLAB data type 
        chunk_as_struct.(field_names{ss}) = sensor_struct; 

    end
    
    % Finalize the convertion 
    converted = chunk_as_struct;

end 