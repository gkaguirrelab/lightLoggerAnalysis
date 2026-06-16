function chunks = parse_chunks(path_to_experiment,...
                               apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats,...
                               apply_phase_correction, apply_RGB_correction, apply_fielding_function,...
                               time_ranges, chunk_ranges, mean_axes, contains_agc_metadata,...
                               verbose, password, differentiate_color...
                              )
% Parse the chunks of a recording from the light logger into a cell of chunk structs
%
% Syntax:
%   chunks = parse_chunks(path_to_experiment)
%   chunks = parse_chunks(path_to_experiment, apply_digital_gain, use_mean_frame, ...)
%
% Description:
%   Reads the desired chunks for a recording using a Python helper
%   function, applying the desired transformations selected, then
%   returns the chunks to MATLAB in a cell array of chunk structs,
%   with all data converted to native MATLAB types.
%
%   When convert_to_floats is false (default), value data preserves
%   hardware-native types: uint8 for camera sensors (W, P) and uint16
%   for the minispect (M). When true, all values are returned as double.
%   Timestamps and AGC settings are always double regardless.
%
% Inputs:
%   path_to_experiment      - String. The path to the suprafile
%                             containing all of the chunks.
%
%   apply_digital_gain      - Logical. Whether or not to apply
%                             each frame's associated digital gain
%                             scalar.
%
%   use_mean_frame          - Logical. Whether to only use the mean
%                             pixel value from each frame for each
%                             of the camera-based sensors.
%
%   convert_time_units      - Logical. Whether to convert the
%                             timestamps from different sensors
%                             all into the same units (s).
%
%   convert_to_floats       - Logical. Whether to convert all sensor
%                             value arrays to double. When false,
%                             camera data is uint8 and MS data is
%                             uint16.
%
%   apply_phase_correction  - Logical. Whether or not to apply
%                             our calculated offsets between sensors
%                             to their time vectors.
%
%   apply_RGB_correction    - Logical. Whether or not to apply
%                             our calculated RGB scalars to world.
%
%   apply_fielding_function - Logical. Whether or not to apply
%                             the fielding function to world.
%
%   time_ranges             - Struct. Represents the selected
%                             range of time to parse. Form is
%                             a struct with fields (WPM)
%                             where each field is a tuple
%                             of a range [start, end).
%
%   chunk_ranges            - Struct. Represents the selected
%                             range of chunks to parse. Form is
%                             a struct with fields (WPM)
%                             where each field is a tuple
%                             of a range [start, end). Note:
%                             the use of this and time_ranges is
%                             exclusively OR. 0-indexed.
%
%   mean_axes               - Struct. Represents the axes to mean
%                             over when use_mean_frame is true.
%                             Only applies to camera sensors. Form
%                             is a struct with fields (WPM)
%                             where each field is a tuple of axes.
%                             0-indexed.
%
%   contains_agc_metadata   - Struct. Represents whether a sensor
%                             has or does not have AGC metadata
%                             for a given recording. Form is
%                             a struct with fields (WPM) where
%                             each field is a boolean.
%
%   verbose                 - Logical. Whether or not to print progress
%                             to the terminal.
%
%   password                - String. Represents the password
%                             used to encrypt the data (if encrypted).
%
%   differentiate_color     - Logical. When true (requires
%                             use_mean_frame=true), separates world
%                             camera Bayer data into R, G, B, and frame
%                             mean channels before spatial averaging.
%                             The final dimension of W.v indexes
%                             [R_mean, G_mean, B_mean, frame_mean],
%                             while the leading dimensions depend on
%                             which axes of the original
%                             [n_frames, n_rows, n_cols] data were
%                             averaged by mean_axes.W. For example,
%                             mean_axes.W=[1,2] yields [n_frames, 4]
%                             and mean_axes.W=[0,1,2] yields [1, 4].
%                             Forces W.v to double. Only applies to
%                             the world camera (W); the pupil camera
%                             (P) is unaffected.
%
% Outputs:
%   chunks                  - Cell. A cell containing structs
%                             of all of the sensors' data
%                             for a given chunk.
%
% Examples:
%{
    % Basic usage with float conversion
    path_to_experiment = '/Volumes/EXTERNAL1/test_folder_0';
    chunks = parse_chunks(path_to_experiment, true, false, true, true, ...
                          true, true, true, false);
%}
%{
    % With color-differentiated world camera means
    % W.v ends in a 4-channel dimension [R, G, B, frame_mean]
    chunks = parse_chunks('/Volumes/EXTERNAL1/test_folder_0', ...
                          false, true, true, false, ...  % dgain=off, mean=on, time=on, float=off
                          false, false, false, ...       % phase=off, rgb=off, field=off
                          false, false, false, false, ...% ranges/axes/agc defaults
                          true, '1234', true);           % verbose, password, differentiate_color
    class(chunks{1}.W.v)  % double (forced by differentiate_color)
    size(chunks{1}.W.v)   % e.g. [n_frames, 4] or [1, 4]
    class(chunks{1}.P.v)  % uint8 (pupil unaffected)
%}

    % Parse and validate the input arguments
    arguments 
        path_to_experiment {mustBeText}; % The path to the suprafolder for this experiment
        apply_digital_gain {mustBeNumericOrLogical} = false; % Whether or not to apply digital gain values to each frame 
        use_mean_frame {mustBeNumericOrLogical} = false; % Whether to use the mean frame from the camera sensors or the entire frames 
        convert_time_units {mustBeNumericOrLogical} = false; % Whether to convert different time units from the different sensors all to seconds 
        convert_to_floats {mustBeNumericOrLogical} = false; % Whether to convert the Python np.arrays to float types. Can make things easier coming to MATLAB
        apply_phase_correction {mustBeNumericOrLogical} = false; % Whether or not to apply our calculated phase offsets to the sensor's t vectors
        apply_RGB_correction {mustBeNumericOrLogical} = false; % Whether or not to apply our RGB scalars to world
        apply_fielding_function {mustBeNumericOrLogical} = false; % Whether or not to apply the fielding function to world
        time_ranges = false; % The timestamps to splice out of a video in the form [start, end) (relative to start of the video)
        chunk_ranges = false; % The chunk numbers to splice out of a video in the form [start, end] (0-indexed)
        mean_axes = false; % The axes per sensor to apply mean over if we want to take some sort of mean. Note: ONLY for camera sensors  
        contains_agc_metadata = false; % Flags for each sensor if its metadata matrix contains AGC data
        verbose {mustBeNumericOrLogical} = true; % If you want progress output to the terminal
        password {mustBeText} = "1234"; % The password needed to decrypt encrypted + compressed files (.blosc files)
        differentiate_color {mustBeNumericOrLogical} = false; % Whether to separate world camera Bayer data into R, G, B, and frame mean channels when use_mean_frame is true
    end

    % differentiate_color requires use_mean_frame
    assert(~differentiate_color || use_mean_frame, ...
        "differentiate_color requires use_mean_frame to be true");

    % Append path to the individual chunk parsing helper function 
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab")); 

    % Load in the Pi util helper file 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path"));

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
    thing.W = {0, py.None};
    thing.M = {0, py.None};
    chunks_as_py = Pi_util.parse_chunks(path_to_experiment,...
                                        apply_digital_gain, use_mean_frame, convert_time_units, convert_to_floats, apply_phase_correction,...
                                        apply_RGB_correction, apply_fielding_function, time_ranges, thing, mean_axes, contains_agc_metadata,...
                                        password, pyargs('differentiate_color', differentiate_color)...
                                       );

    % Convert outer Python list to cell 
    chunks = cell(chunks_as_py);

    % Iterate over the chunks and convert them to fully MATLAB types 
    for cc = 1:numel(chunks)
        if(verbose)
            fprintf("Converting chunk: %d/%d\n", cc, numel(chunks));
        end 

        % Replace the Python data type with the MATLAB data type. 
        chunks{cc} = chunk_dict_to_matlab_lla(chunks{cc}, 'convert_to_floats', convert_to_floats, 'differentiate_color', differentiate_color);
    end 
    
end
