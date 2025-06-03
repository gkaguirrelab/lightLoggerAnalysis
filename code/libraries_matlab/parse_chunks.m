function chunks = parse_chunks(path_to_experiment, ...
                               apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats,...
                               time_ranges, chunk_ranges, mean_axes, contains_agc_metadata,...
                               password,...
                               matlab_analysis_libraries_path,...
                               Pi_util_path...
                              )
% Parse the chunks of a recording from the light logger into a cell of chunk structs
%
% Syntax:
% function chunks = parse_chunks(path_to_experiment, ...
%                                apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats,...
%                                time_ranges, chunk_ranges, mean_axes, contains_agc_metadata,...
%                                password,...
%                                matlab_libraries_path,...
%                                Pi_util_path...
%                               )
%
% Description:
%   Readings in the desired chunks for a recording using a Python helper 
%   function, applying the desired transformations selected, then
%   returns the chunks to MATLAB in a cell array of chunk structs, 
%   with all data converted to native MATLAB types. 
%
% Inputs:
%   path_to_experiment    - String. The path to the suprafile 
%                           containing all of the chunks.
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
%   time_ranges           - Struct. Represents the selected 
%                           range of time to parse. Form is 
%                           a struct with fields (WPM)
%                           where each field is a tuple 
%                           of a range [start, end). 
%
%   chunk_ranges          - Struct. Represents the selected 
%                           range of chunks to parse. Form is 
%                           a struct with fields (WPM)
%                           where each field is a tuple 
%                           of a range [start, end). Note:
%                           the use of this and time_ranges is
%                           exclusively OR. 0-indexed. 
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
%   matlab_analysis_libraries_path  - String. Path to the utilized MATLAB 
%                                     helper functions (chunk_dict_to_matlab)
%                                       
%
%   Pi_util_path           - String. Path to the Pi_util.py helper file
%
%
% Outputs:
%
%   chunks                - Cell. A cell containing structs 
%                           of all of the sensors' data 
%                           for a given chunk.  
%
% Examples:
%{
    path_to_experiment = '/Volumes/EXTERNAL1/test_folder_0';
    chunks = parse_chunks(path_to_experiment); 
%}

    % Parse and validate the input arguments
    arguments 
        path_to_experiment {mustBeText}; % The path to the suprafolder for this experiment
        apply_digital_gain {mustBeNumericOrLogical} = false; % Whether or not to apply digital gain values to each frame 
        use_mean_frame {mustBeNumericOrLogical} = false; % Whether to use the mean frame from the camera sensors or the entire frames 
        convert_time_units {mustBeNumericOrLogical} = false; % Whether to convert different time units from the different sensors all to seconds 
        convert_to_floats {mustBeNumericOrLogical} = false; % Whether to convert the Python np.arrays to float types. Can make things easier coming to MATLAB
        time_ranges = false; % The timestamps to splice out of a video in the form [start, end) (relative to start of the video)
        chunk_ranges = false; % The chunk numbers to splice out of a video in the form [start, end] (0-indexed)
        mean_axes = false; % The axes per sensor to apply mean over if we want to take some sort of mean. Note: ONLY for camera sensors  
        contains_agc_metadata = false; % Flags for each sensor if its metadata matrix contains AGC data
        password {mustBeText} = "1234"; % The password needed to decrypt encrypted + compressed files (.blosc files)
        matlab_analysis_libraries_path {mustBeText} = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "libraries_matlab"); % Path to the utilized MATLAB helper functions (chunk_dict_to_matlab)
        Pi_util_path {mustBeText} = fullfile(fileparts(mfilename("fullpath")), 'Pi_util');  % Path to the Pi_util.py helper file
    end 

    % Append path to the individual chunk parsing helper function 
    matlab_libraries_path = matlab_analysis_libraries_path;
    addpath(matlab_libraries_path); 

    % Load in the Pi util hepler file 
    Pi_util = import_pyfile(Pi_util_path);

    % Apply the default time ranges to splice out of the video if not supplied 
    if(~isstruct(time_ranges))
        % Default is a 1x2 None tuple. This signifies the entire video
        time_ranges = struct; 
        
        % Splice the entire video from all sensors 
        sensor_names = {'W', 'P', 'M'}; 
        for ss = 1:numel(sensor_names)  
            time_ranges.(sensor_names{ss}) = py.tuple({py.None, py.None}); 
        end 
    end 

    % Apply the default chunk ranges to splice out of the video if not supplied.
    % Note: this XOR timestamps can be used 
    if(~isstruct(chunk_ranges))
        % Default is 1x2 tuple (0, None). This signifies the entire video
        chunk_ranges = struct; 

        % Splice out the entire video from all sensors 
        sensor_names = {'W', 'P', 'M'}; 
        for ss = 1:numel(sensor_names)  
            chunk_ranges.(sensor_names{ss}) = py.tuple({0, py.None}); 
        end 

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

    % Parse the experiment using Python and return as py.list
    chunks_as_py = Pi_util.parse_chunks(path_to_experiment,...
                                        apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats,...
                                        time_ranges, chunk_ranges, mean_axes, contains_agc_metadata,...
                                        password...
                                       );

    % Convert outer Python list to cell 
    chunks = cell(chunks_as_py);

    % Iterate over the chunks and convert them to fully MATLAB types 
    for cc = 1:numel(chunks)
        fprintf("Converting chunk: %d/%d\n", cc, numel(chunks));

        % Replace the Python data type with the MATLAB data type. 
        chunks{cc} = chunk_dict_to_matlab(chunks{cc}); 
    end 
    
end
