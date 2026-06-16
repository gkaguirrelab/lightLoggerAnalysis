    function LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder, options)
    % Load, parse, and convert light logger calibration data from Python to MATLAB
    %
    % Syntax:
    %   LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder)
    %   LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder, 'apply_digital_gain', true)
    %   LightLoggerCalibrationData = convert_light_logger_calibration_data(experiment_folder, Name, Value)
    %
    % Description:
    %   Given a path to a calibration experiment folder (containing NDF
    %   subfolders and a calibration_metadata.mat.bytes file), this function
    %   retrieves the sorted folder paths for each calibration operation via
    %   a Python helper (calibration_util.py), then parses and converts each
    %   individual measurement one at a time in MATLAB. This approach keeps
    %   memory usage low by only holding one measurement's raw data in memory
    %   at a time, converting it to MATLAB types, and freeing the Python
    %   memory before moving to the next measurement.
    %
    %   Each calibration operation's readings are stored in a multi-dimensional
    %   cell array indexed by NDF level, stimulus parameter(s), and measurement
    %   number. Individual cells contain chunk structs with fields for each
    %   sensor (W, P, M, S), where each sensor struct has t (timestamps),
    %   v (values), and optionally settings (AGC metadata).
    %
    %   When convert_to_floats is false (default), value data preserves
    %   hardware-native types: uint8 for camera sensors (W, P) and uint16
    %   for the minispect (M). When true, all values are returned as double.
    %   Timestamps and AGC settings are always returned as double regardless.
    %
    % Inputs:
    %   experiment_folder           - String. Path to the top-level calibration
    %                                 experiment directory. Must contain one or
    %                                 more NDF_* subfolders with recording data
    %                                 and a calibration_metadata.mat.bytes file.
    %
    % Name-Value Arguments:
    %   'apply_digital_gain'        - Logical (default: false). Whether to
    %                                 multiply each camera frame by its
    %                                 associated digital gain scalar from the
    %                                 AGC metadata.
    %
    %   'apply_radiometric_correction' - Logical (default: false). Whether to
    %                                 apply per-channel RGB radiometric scalars
    %                                 to the world camera Bayer data.
    %
    %   'convert_time_units'        - Logical (default: false). Whether to
    %                                 convert timestamps from each sensor's
    %                                 native units into seconds. World camera
    %                                 timestamps are in nanoseconds; MS and
    %                                 pupil timestamps are already in seconds.
    %
    %   'convert_to_floats'         - Logical (default: false). Whether to
    %                                 cast all sensor value arrays to double.
    %                                 When false, camera data is uint8 and MS
    %                                 data is uint16.
    %
    %   'use_mean_frame'            - Logical (default: false). Whether to
    %                                 reduce camera frame buffers by averaging
    %                                 over the spatial axes specified by
    %                                 mean_axes.
    %
    %   'mean_axes'                 - Struct or false (default: false). Axes
    %                                 over which to average when use_mean_frame
    %                                 is true. Struct with fields W, P, M where
    %                                 each value is a vector of 0-indexed axes.
    %                                 If false, defaults to W=[1,2], P=[1,2],
    %                                 M=[].
    %
    %   'differentiate_color'       - Logical (default: false). When true
    %                                 (requires use_mean_frame=true), the world
    %                                 camera Bayer data is separated into R, G,
    %                                 B, and overall frame mean channels before
    %                                 spatial averaging. The final dimension of
    %                                 W.v indexes [R_mean, G_mean, B_mean,
    %                                 frame_mean], while the leading dimensions
    %                                 depend on which axes in the original
    %                                 [n_frames, n_rows, n_cols] data were
    %                                 averaged by mean_axes.W. For example,
    %                                 mean_axes.W=[1,2] yields [n_frames, 4],
    %                                 while mean_axes.W=[0,1,2] yields [1, 4].
    %                                 Forces W.v to double regardless of
    %                                 convert_to_floats. Only applies to the
    %                                 world camera (W); the pupil camera (P) is
    %                                 unaffected.
    %
    %   'verbose'                   - Logical (default: true). Whether to print
    %                                 progress messages during parsing.
    %
    % Outputs:
    %   LightLoggerCalibrationData  - Struct with two fields:
    %       .metadata               - Struct. The calibration metadata loaded
    %                                 from calibration_metadata.mat.bytes,
    %                                 containing stimulus parameters, NDF
    %                                 levels, randomization orders, etc. for
    %                                 each calibration operation.
    %       .readings               - Struct with fields:
    %           .ms_linearity       - Cell array (NDF x settings x measure).
    %                                 Minispect linearity measurements at
    %                                 varying background scalars.
    %           .temporal_sensitivity - Cell array (NDF x contrast x freq x
    %                                 measure). Temporal sensitivity
    %                                 measurements de-randomized to canonical
    %                                 contrast and frequency order.
    %           .phase_fitting      - Cell array (NDF x contrast x freq x
    %                                 measure). Phase fitting measurements,
    %                                 same structure as temporal_sensitivity.
    %           .contrast_gamma     - Cell array (NDF x contrast x freq x
    %                                 measure). Contrast gamma measurements,
    %                                 same structure as temporal_sensitivity.
    %           .world_linearity    - Cell array (NDF x contrast_target x
    %                                 settings x measure). World camera
    %                                 linearity measurements.
    %
    % Examples:
    %{
        % Basic usage: load calibration data with default settings
        experiment_folder = '/Volumes/T7 Shield/5cameraLinearity';
        CalData = convert_light_logger_calibration_data(experiment_folder);

        % Access the first ms_linearity measurement at NDF=1, settings=1
        chunk = CalData.readings.ms_linearity{1, 1, 1};
        disp(class(chunk.M.v.AS));  % uint16

        % Access the first world_linearity measurement
        chunk = CalData.readings.world_linearity{1, 1, 1, 1};
        disp(class(chunk.W.v));  % uint8
    %}
    %{
        % Load with digital gain correction and float conversion
        CalData = convert_light_logger_calibration_data(...
            '/Volumes/T7 Shield/5cameraLinearity', ...
            'apply_digital_gain', true, ...
            'convert_to_floats', true, ...
            'convert_time_units', true ...
        );

        % All values are now double
        chunk = CalData.readings.temporal_sensitivity{1, 1, 1, 1};
        disp(class(chunk.W.v));  % double
    %}
    %{
        % Load with mean frame reduction and custom axes
        mean_axes = struct('W', [1, 2], 'P', [1, 2], 'M', []);
        CalData = convert_light_logger_calibration_data(...
            '/Volumes/T7 Shield/5cameraLinearity', ...
            'use_mean_frame', true, ...
            'mean_axes', mean_axes, ...
            'apply_radiometric_correction', true, ...
            'verbose', false ...
        );
    %}
    %{
        % Load with Bayer color channel separation
        % W.v ends in a 4-channel dimension [R, G, B, frame_mean]
        CalData = convert_light_logger_calibration_data(...
            '/Volumes/T7 Shield/5cameraLinearity', ...
            'use_mean_frame', true, ...
            'differentiate_color', true ...
        );

        chunk = CalData.readings.temporal_sensitivity{1, 1, 1, 1};
        size(chunk.W.v)       % e.g. [n_frames, 4] or [1, 4]
        disp(class(chunk.W.v))  % double (forced by differentiate_color)
        disp(class(chunk.P.v))  % uint8 (pupil unaffected)
    %}

    arguments
        experiment_folder {mustBeText};
        options.apply_digital_gain {mustBeNumericOrLogical} = false;
        options.apply_radiometric_correction {mustBeNumericOrLogical} = false;
        options.convert_time_units {mustBeNumericOrLogical} = false;
        options.convert_to_floats {mustBeNumericOrLogical} = false;
        options.use_mean_frame {mustBeNumericOrLogical} = false;
        options.mean_axes = false;
        options.differentiate_color {mustBeNumericOrLogical} = false;
        options.verbose {mustBeNumericOrLogical} = true;
    end

    % differentiate_color requires use_mean_frame
    assert(~options.differentiate_color || options.use_mean_frame, ...
        "differentiate_color requires use_mean_frame to be true");

    mean_axes = options.mean_axes;

    % If we have not been passed in specific axis to mean, use the default
    % which is to mean each frame
    if(~isstruct(mean_axes))
        mean_axes = struct;
        mean_axes.W = [1, 2];
        mean_axes.P = [1, 2];
        mean_axes.M = [];
    end

    % Import Python libraries
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab"));
    calibration_util = import_pyfile(getpref("lightLoggerAnalysis", "calibration_util_path"));
    chunk_io = import_pyfile(getpref("lightLoggerAnalysis", "chunk_io_path"));

    % Retrieve sorted folder paths only (no parsing in Python)
    if(options.verbose) 
        disp(experiment_folder)
    end

    sorted_paths = struct(calibration_util.load_sorted_calibration_files(experiment_folder, pyargs('parse_files', false)));

    % Save the path to the current file's directory. we are going to need this when we mess around with tbUses
    path_to_file_dir = fileparts(mfilename("fullpath"));

    % Load in the Calibration Metadata
    % Include dependencies (without this the rawData structs of the cal field)
    % of the metadata struct will be empty
    tbUse('combiLEDToolbox');
    addpath(path_to_file_dir);

    calibration_metadata_path = fullfile(experiment_folder, "calibration_metadata.mat.bytes"); 
    if(options.verbose)
        fprintf("Calbibration metadata path: %s\n", calibration_metadata_path)
    end

    calibration_metadata = load_calibration_metadata(calibration_metadata_path);
    tbUseProject('lightLoggerAnalysis');

    % Build the parse options struct to pass to subfunctions
    parse_opts.apply_digital_gain = options.apply_digital_gain;
    parse_opts.apply_radiometric_correction = options.apply_radiometric_correction;
    parse_opts.convert_time_units = options.convert_time_units;
    parse_opts.convert_to_floats = options.convert_to_floats;
    parse_opts.use_mean_frame = options.use_mean_frame;
    parse_opts.mean_axes = mean_axes;
    parse_opts.differentiate_color = options.differentiate_color;
    parse_opts.verbose = options.verbose;

    % Parse and convert each calibration operation one measurement at a time
    parsed_readings.ms_linearity = convert_ms_linearity_to_matlab(calibration_metadata.ms_linearity, sorted_paths.ms_linearity, chunk_io, parse_opts);
    parsed_readings.temporal_sensitivity = convert_temporal_sensitivity_to_matlab(calibration_metadata.temporal_sensitivity, sorted_paths.temporal_sensitivity, chunk_io, parse_opts);
    parsed_readings.phase_fitting = convert_temporal_sensitivity_to_matlab(calibration_metadata.phase_fitting, sorted_paths.phase_fitting, chunk_io, parse_opts);
    parsed_readings.contrast_gamma = convert_temporal_sensitivity_to_matlab(calibration_metadata.contrast_gamma, sorted_paths.contrast_gamma, chunk_io, parse_opts);
    parsed_readings.world_linearity = convert_camera_linearity_to_matlab(calibration_metadata.world_linearity, sorted_paths.world_linearity, chunk_io, parse_opts);

    % Initialize a return struct
    LightLoggerCalibrationData.metadata = calibration_metadata;
    LightLoggerCalibrationData.readings = parsed_readings;

    % Clear any python memory
    py.gc.collect();

    return ;

end


function converted_linearity = convert_ms_linearity_to_matlab(linearity_calibration_metadata, ms_linearity_paths, chunk_io, opts)

    num_NDF_levels = numel(linearity_calibration_metadata.NDFs);
    num_settings_levels = numel(linearity_calibration_metadata.background_scalars);
    n_measures = linearity_calibration_metadata.n_measures;

    converted_linearity = cell(num_NDF_levels, num_settings_levels, n_measures);

    if(opts.verbose); fprintf("[MS Linearity]\n"); end

    % Convert the outer Python list to cell array of NDF levels
    NDF_cells = cell(ms_linearity_paths);

    for nn = 1:num_NDF_levels
        if(opts.verbose); fprintf("NDF %d/%d\n", nn, num_NDF_levels); end

        % Convert this NDF's settings-level list to cell
        settings_cells = cell(NDF_cells{nn});

        for ss = 1:num_settings_levels
            if(opts.verbose); fprintf("  Settings level %d/%d\n", ss, num_settings_levels); end

            % Convert the measurements list at this settings level to cell
            measurement_cells = cell(settings_cells{ss});

            for mm = 1:n_measures
                if(opts.verbose); fprintf("    Measurement %d/%d\n", mm, n_measures); end

                % Get the path for this measurement
                measurement_path = string(measurement_cells{mm});

                % Parse this single measurement in Python
                chunks_py = chunk_io.parse_chunks(measurement_path, ...
                    pyargs( ...
                        'apply_digital_gain', opts.apply_digital_gain, ...
                        'use_mean_frame', opts.use_mean_frame, ...
                        'convert_time_units', opts.convert_time_units, ...
                        'convert_to_float', opts.convert_to_floats, ...
                        'apply_RGB_correction', opts.apply_radiometric_correction, ...
                        'mean_axes', opts.mean_axes, ...
                        'differentiate_color', opts.differentiate_color ...
                    ) ...
                );

                % Convert to MATLAB types one chunk at a time
                chunks_cell = cell(chunks_py);
                clear chunks_py;
                for ci = 1:numel(chunks_cell)
                    chunks_cell{ci} = chunk_dict_to_matlab_lla(chunks_cell{ci}, ...
                        'convert_to_floats', opts.convert_to_floats, ...
                        'differentiate_color', opts.differentiate_color);
                end

                converted_linearity{nn, ss, mm} = chunks_cell{:};

                % Free Python memory immediately
                clear chunks_cell;
                py.gc.collect();
            end
        end
    end

end


function converted_temporal_sensitivity = convert_temporal_sensitivity_to_matlab(temporal_sensitivity_calibration_metadata, temporal_sensitivity_paths, chunk_io, opts)

    num_NDF_levels = numel(temporal_sensitivity_calibration_metadata.NDFs);
    num_contrast_levels = numel(temporal_sensitivity_calibration_metadata.contrast_levels);
    num_frequencies = numel(temporal_sensitivity_calibration_metadata.frequencies);
    n_measures = temporal_sensitivity_calibration_metadata.n_measures;

    converted_temporal_sensitivity = cell(num_NDF_levels, num_contrast_levels, num_frequencies, n_measures);

    contrast_orders = temporal_sensitivity_calibration_metadata.contrast_levels_orders;
    frequencies_orders = temporal_sensitivity_calibration_metadata.frequencies_orders;

    if(opts.verbose); fprintf("[Temporal Sensitivity / Phase Fitting / Contrast Gamma]\n"); end

    % Convert the outer Python list to cell array of NDF levels
    NDF_cells = cell(temporal_sensitivity_paths);

    for nn = 1:num_NDF_levels
        if(opts.verbose); fprintf("NDF %d/%d\n", nn, num_NDF_levels); end

        % Convert this NDF's contrast-level list to cell
        contrast_cells = cell(NDF_cells{nn});

        for cc = 1:num_contrast_levels
            if(opts.verbose); fprintf("  Contrast level %d/%d\n", cc, num_contrast_levels); end

            % Convert the frequencies list at this contrast to cell
            frequency_cells = cell(contrast_cells{cc});

            for ff = 1:num_frequencies
                if(opts.verbose); fprintf("    Frequency %d/%d\n", ff, num_frequencies); end

                % Convert the measurements list at this frequency to cell
                measurement_cells = cell(frequency_cells{ff});

                for mm = 1:n_measures
                    if(opts.verbose); fprintf("      Measurement %d/%d\n", mm, n_measures); end

                    % Get the path for this measurement
                    measurement_path = string(measurement_cells{mm});

                    % Parse this single measurement in Python
                    chunks_py = chunk_io.parse_chunks(measurement_path, ...
                        pyargs( ...
                            'apply_digital_gain', opts.apply_digital_gain, ...
                            'use_mean_frame', opts.use_mean_frame, ...
                            'convert_time_units', opts.convert_time_units, ...
                            'convert_to_float', opts.convert_to_floats, ...
                            'apply_RGB_correction', opts.apply_radiometric_correction, ...
                            'mean_axes', opts.mean_axes, ...
                            'contains_agc_metadata_dict', py.dict(pyargs('W', true, 'P', false, 'M', false)), ...
                            'differentiate_color', opts.differentiate_color ...
                        ) ...
                    );

                    % Convert to MATLAB types one chunk at a time
                    chunks_cell = cell(chunks_py);
                    clear chunks_py;
                    for ci = 1:numel(chunks_cell)
                        chunks_cell{ci} = chunk_dict_to_matlab_lla(chunks_cell{ci}, ...
                            'convert_to_floats', opts.convert_to_floats, ...
                            'differentiate_color', opts.differentiate_color);
                    end

                    % De-randomize: find the canonical index for this contrast and frequency
                    contrast_idx = contrast_orders(nn, mm, cc);
                    frequency_idx = frequencies_orders(nn, mm, cc, ff);

                    converted_temporal_sensitivity{nn, contrast_idx, frequency_idx, mm} = chunks_cell{:};

                    % Free Python memory immediately
                    clear chunks_cell;
                    py.gc.collect();
                end
            end
        end
    end

end


function converted_linearity = convert_camera_linearity_to_matlab(linearity_calibration_metadata, world_linearity_paths, chunk_io, opts)

    num_NDF_levels = numel(linearity_calibration_metadata.NDFs);
    num_contrast_targets = numel(linearity_calibration_metadata.contrast_agc_targets);
    num_settings_levels = numel(linearity_calibration_metadata.background_scalars);
    n_measures = linearity_calibration_metadata.n_measures;

    converted_linearity = cell(num_NDF_levels, num_contrast_targets, num_settings_levels, n_measures);

    if(opts.verbose); fprintf("[World Camera Linearity]\n"); end

    % Convert the outer Python list to cell array of NDF levels
    NDF_cells = cell(world_linearity_paths);

    for nn = 1:num_NDF_levels
        if(opts.verbose); fprintf("NDF %d/%d\n", nn, num_NDF_levels); end

        % Convert this NDF's contrast-target list to cell
        contrast_cells = cell(NDF_cells{nn});

        for cc = 1:num_contrast_targets
            if(opts.verbose); fprintf("  Contrast target %d/%d\n", cc, num_contrast_targets); end

            % Convert the settings list at this contrast target to cell
            settings_cells = cell(contrast_cells{cc});

            for ss = 1:num_settings_levels
                if(opts.verbose); fprintf("    Settings level %d/%d\n", ss, num_settings_levels); end

                % Convert the measurements list at this settings level to cell
                measurement_cells = cell(settings_cells{ss});

                for mm = 1:n_measures
                    if(opts.verbose); fprintf("      Measurement %d/%d\n", mm, n_measures); end

                    % Get the path for this measurement
                    measurement_path = string(measurement_cells{mm});

                    % Parse this single measurement in Python
                    chunks_py = chunk_io.parse_chunks(measurement_path, ...
                        pyargs( ...
                            'apply_digital_gain', opts.apply_digital_gain, ...
                            'use_mean_frame', opts.use_mean_frame, ...
                            'convert_time_units', opts.convert_time_units, ...
                            'convert_to_float', opts.convert_to_floats, ...
                            'apply_RGB_correction', opts.apply_radiometric_correction, ...
                            'mean_axes', opts.mean_axes, ...
                            'differentiate_color', opts.differentiate_color ...
                        ) ...
                    );

                    % Convert to MATLAB types one chunk at a time
                    chunks_cell = cell(chunks_py);
                    clear chunks_py;
                    for ci = 1:numel(chunks_cell)
                        chunks_cell{ci} = chunk_dict_to_matlab_lla(chunks_cell{ci}, ...
                            'convert_to_floats', opts.convert_to_floats, ...
                            'differentiate_color', opts.differentiate_color);

                        if opts.verbose && isfield(chunks_cell{ci}, 'W') && isfield(chunks_cell{ci}.W, 'v')
                            world_values = chunks_cell{ci}.W.v;
                            world_size = sprintf('%dx', size(world_values));
                            world_size = world_size(1:end-1);
                            fprintf("        Chunk %d/%d W.v size: [%s], class: %s\n", ...
                                ci, numel(chunks_cell), world_size, class(world_values));
                        end
                    end

                    converted_linearity{nn, cc, ss, mm} = chunks_cell{:};

                    % Free Python memory immediately
                    clear chunks_cell;
                    py.gc.collect();
                end
            end
        end
    end

end
