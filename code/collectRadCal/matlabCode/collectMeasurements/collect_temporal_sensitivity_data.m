function [success, temporal_sensitivity_calibration_metadata] = collect_temporal_sensitivity_data(device_num, ...
                                                                                                  temporal_sensitivity_calibration_metadata,...
                                                                                                  bluetooth_central,...
                                                                                                  bluetooth_client,...
                                                                                                  label,...
                                                                                                  dropbox_savedir,...
                                                                                                  local_savedir,...
                                                                                                  cooldown_callback...
                                                                                                ) 
% Collect calibration data for temporal sensitivity of 
% all of the sensors in the light logger
%
% Syntax:
%   [success, temporal_sensitivity_calibration_metadata] = collect_temporal_sensitivity_data(device_num, temporal_sensitivity_calibration_metadata, bluetooth_central, bluetooth_client, label, dropbox_savedir, local_savedir, cooldown_callback)
%
% Description:
%   This function coordinates the temporal-sensitivity calibration by
%   iterating over NDF levels, randomized contrast orders, randomized
%   frequency orders within each contrast, and repeated measurements. For
%   each condition it configures the CombiLED modulation, starts the
%   flicker, triggers a light-logger calibration recording over Bluetooth,
%   waits until the peripheral returns to its idle state, stops the
%   modulation, and records completion in the metadata so interrupted runs
%   can later resume from the remaining conditions.
%
% Inputs:
%   device_num            - Scalar. Identifier of the target light logger
%                           peripheral.
%   temporal_sensitivity_calibration_metadata
%                         - Struct. Calibration plan and progress record
%                           for the temporal-sensitivity stage.
%   bluetooth_central     - Python module exposing the Bluetooth control
%                           wrappers.
%   bluetooth_client      - Persistent BLE client handle returned by the
%                           connection helper.
%   label                 - String. Optional filename prefix for the
%                           collected recordings.
%   dropbox_savedir       - String. Optional cloud output directory.
%   local_savedir         - String. Optional local output directory.
%   cooldown_callback     - Function handle invoked between NDF blocks
%                           when more unfinished work remains.
%
% Outputs:
%
%   success               - Logical scalar. Returns true on success and
%                           false on failure.
%   temporal_sensitivity_calibration_metadata
%                         - Updated metadata struct with completion flags,
%                           stored modulations, and any saved error
%                           message.
%
% Examples:
%{
    [success, temporal_meta] = collect_temporal_sensitivity_data( ...
        device_num, temporal_meta, bluetooth_central, bluetooth_client, ...
        "", cloud_output_dir, local_output_dir, ...
        @wait_for_combiLED_cooldown_between_NDFs);
%}

    % Validate the arguments 
    arguments
        device_num;
        temporal_sensitivity_calibration_metadata;
        bluetooth_central;
        bluetooth_client;
        label = "";
        dropbox_savedir = "";
        local_savedir = "";
        cooldown_callback = [];
    end 

    % First, extract some information from the calibration struct 
    temporal_sensitivity_calibration_metadata.last_error_message = "";
    NDFs = temporal_sensitivity_calibration_metadata.NDFs; 
    background = temporal_sensitivity_calibration_metadata.background; 
    contrast_levels = temporal_sensitivity_calibration_metadata.contrast_levels; 
    contrast_levels_orders = temporal_sensitivity_calibration_metadata.contrast_levels_orders; 
    frequencies = temporal_sensitivity_calibration_metadata.frequencies; 
    frequencies_orders = temporal_sensitivity_calibration_metadata.frequencies_orders; 
    n_seconds = temporal_sensitivity_calibration_metadata.recording_seconds; 
    n_measures = temporal_sensitivity_calibration_metadata.n_measures;  

    % First, we will iterate across the NDF levels 
    for nn = 1:numel(NDFs)
        % Retrieve the current NDF 
        NDF = NDFs(nn); 

        % Determine if there are any measurements to be collected 
        % If there are not, we can simply skip this NDF 
        completed_measurements_NDF = temporal_sensitivity_calibration_metadata.completed_measurements(nn, :, :, :);
        if(all(completed_measurements_NDF(:) == true))
            continue 
        end 

        % Initialize the CombiLED for this NDF level 
        % Ensure we have constants from lightlogger
        tbUseProject('lightLogger'); 
        [CombiLED, cal] = initialize_combiLED_lightLogger(getpref("lightLogger", "combiExperiments_path"), round(NDF)); 

        % Save the cal file for this NDF 
        cal_files = temporal_sensitivity_calibration_metadata.cal_files;
        cal_files{nn} = cal; 
        temporal_sensitivity_calibration_metadata.cal_files = cal_files; 

        % Construct the modulation to dispaly
        observerAgeInYears = 30; % Build general information about the modulation on the CombiLED 
        pupilDiameterMm = 3; 

        photoreceptors = photoreceptorDictionaryHuman('observerAgeInYears',observerAgeInYears,'pupilDiameterMm',pupilDiameterMm);
        modResult = designModulation('LightFlux', photoreceptors, cal); % Construct the modulation 
        modResult.settingsHigh = (background)'; % Set the color profile of the modulation 
        modResult.settingsBackground = (background./2)';
        modResult.settingsLow = (background.*0)';

        % Save the modulation for this NDF
        modulations = temporal_sensitivity_calibration_metadata.modulations;
        modulations{nn} = modResult; 
        temporal_sensitivity_calibration_metadata.modulations = modulations; 

        % Prepare the CombiLED for the modulation 
        CombiLED.setSettings(modResult);
        pause(0.5); 
        CombiLED.setWaveformIndex(1);
        pause(0.5);

        % Retrieve the sensors and settings for this NDF level 
        sensors = temporal_sensitivity_calibration_metadata.sensors_and_settings{nn};

        % Pause to ensure the NDF is attached 
        fprintf("Preparing to capture. Ensure NDF: %f is attached. Press any key to continue: \n", NDF); 
        pause(); 
    
        % Then, we will perform 3 iterations of measurements 
        for mm = 1:n_measures
            % Here, we will select the order of the contrasts 
            % that will be shown 
            contrast_order = squeeze(contrast_levels_orders(nn, mm, :));

            % Put the contrasts in this order 
            shuffled_contrasts = contrast_levels(contrast_order); 

            assert( numel(shuffled_contrasts) == numel(contrast_levels) ); 

            % Next, iterate over the shuffled contrast levels
            for cc = 1:numel(shuffled_contrasts)
                % Retrieve the contrast level that we will expose 
                contrast = shuffled_contrasts(cc); 
                contrast_idx = contrast_order(cc);
                
                % For each contrast level, 
                % randomize the frequency order  
                frequencies_order = squeeze(frequencies_orders(nn, mm, contrast_idx, :));

                % Put the frequencies in this order 
                shuffled_frequencies = frequencies(frequencies_order);    

                assert(numel(shuffled_frequencies) == numel(frequencies) ) ; 

                % Iterate over the shuffled frequencies 
                for ff = 1:numel(shuffled_frequencies)
                    % Retrieve the frequency we will expose 
                    freq = shuffled_frequencies(ff); 
                    freq_idx = frequencies_order(ff);

                    % Note that the CombiLED has a roll-off in modulation depth
                    % with temporal frequency. We do not correct for this here,
                    % but instead apply the correction at the data analysis
                    % stage.

                    % Set the contrast level on the CombiLED 
                    CombiLED.setContrast(contrast);
                    pause(0.5);

                    % Prepare the CombiLED to flicker at this frequency 
                    CombiLED.setFrequency(freq);
                    pause(0.5);

                    fprintf("Temporal Sensitivity | NDF #: (%d/%d) Contrast #: (%d/%d) / Freq #: (%d/%d) C: %.3f F: %.3f | Measurement #(%d/%d)\n", ...
                        nn, numel(NDFs), cc, numel(contrast_levels), ff, numel(frequencies), contrast, freq, mm, n_measures);
                    
                    % If we have already performed this combination of contrast frequency 
                    % and measure, skip it 
                    if(temporal_sensitivity_calibration_metadata.completed_measurements(nn, contrast_idx, freq_idx, mm))
                        pause(0.5);
                        continue ; 
                    end 

                    % Start the modulation on the CombiLED 
                    CombiLED.startModulation();
                    pause(0.5);

                    try
                        % Generate a message to send to the RPi
                        update_message = py_call_module_attr(bluetooth_central, "initialize_update_message"); % Initialize the message dict 
                        
                        filename = label + sprintf("TemporalSensitivity_%dcontrastIdx_%dfreqIdx_%dmeasurementIdx", contrast_idx, freq_idx, mm); % Define the filename for the recording (replace. in floats with x)

                        % Compose the target output dir for this NDF level
                        cloud_output_dir = ""; 
                        if(dropbox_savedir ~= "")
                            cloud_output_dir = fullfile(dropbox_savedir, sprintf("NDF%f", NDF));
                        end 
                        
                        local_output_dir = "";
                        if(local_savedir ~= "")
                            local_output_dir = fullfile(local_savedir, sprintf("NDF%f", NDF));
                        end 


                        % If the desired recording seconds is not long enough to see a full modulation, 
                        % we will edit it to be the max single chunk length 
                        recording_seconds = n_seconds; 
                        if(freq * recording_seconds < 1)
                            warning(sprintf("Recording seconds: %d not enough to see a single modulation at %.2f, defaulting to max length (29)", recording_seconds, freq));

                            recording_seconds = 29; 
                        end 

                        % Then, let's check if this is enough to see 1 modulation. If not, we error 
                        if(freq * recording_seconds < 1)
                            fprintf("ERROR: Recording seconds: %d is not enough to see a single modulation at %.2f\n", recording_seconds, freq);

                            success = 0; 
                            return; 
                        end 

                        py_call_module_attr(bluetooth_central, "generate_calibration_state", update_message, py.str(filename),... % Put it all together in the dict 
                                            cloud_output_dir, local_output_dir,...
                                            true, py.int(recording_seconds),... 
                                            py.int(30), sensors);

                        % Send a message to the RPi to make a recording 
                        py_call_module_attr(bluetooth_central, "message_peripheral_matlab_wrapper", device_num, update_message, bluetooth_client);

                        % Read the state from the light logger, and pause until it is changed 
                        % from calibration to wait 
                        % if it changes to error, then throw an error
                        while(true)
                            % Read the current state from the light logger.
                            % Transient CoreBluetooth read errors (e.g. CBATT
                            % "Unlikely error", code 14) are retried rather than
                            % aborting the whole calibration, since the recording
                            % has already been triggered on the light logger.
                            lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client);
                            state_name = string(char(lightlogger_state.state));     

                            % If the state is error, raise an error on this machine 
                            % as well 
                            if(state_name == "error") 
                                temporal_sensitivity_calibration_metadata.last_error_message = describe_lightlogger_error(lightlogger_state);
                                success = 0 ; 
                                return ; 
                            end 

                            % If the light logger is done recording + uploading 
                            % break from the loop
                            if(state_name == "wait")
                                break 
                            end
                            
                            % Pause for some time between reads
                            pause(0.5); 
                        end
                    catch ME
                        try
                            CombiLED.stopModulation();
                        catch
                        end
                        report_bluetooth_failure("Temporal sensitivity", NDF, contrast_idx, freq_idx, mm, ME);
                        temporal_sensitivity_calibration_metadata.last_error_message = ...
                            sprintf("Temporal sensitivity bluetooth failure at NDF %.3f contrast %d frequency %d measurement %d.\n%s", ...
                                    NDF, contrast_idx, freq_idx, mm, getReport(ME, "extended", "hyperlinks", "off"));
                        success = 0;
                        return;
                    end

                    % End the modulation 
                    CombiLED.stopModulation(); 
                    pause(0.5);

                    % Otherwise, mark this measurement as completed
                    temporal_sensitivity_calibration_metadata.completed_measurements(nn, contrast_idx, freq_idx, mm) = true; 

                end % Frequencies 

            end % Contrasts

        end % Measures 

        if(should_cooldown_before_next_NDF(temporal_sensitivity_calibration_metadata.completed_measurements, nn) && ~isempty(cooldown_callback))
            cooldown_callback(CombiLED, NDF);
        end

        % Close the CombiLED for this NDF 
        CombiLED.serialClose(); 

    end % NDF

    % Exit with success if nothing went awry
    success = 1; 

    return ; 

end 

function tf = should_cooldown_before_next_NDF(completed_measurements, current_NDF_idx)
% Decide whether another unfinished NDF remains after the current one.
%
% Syntax:
%   tf = should_cooldown_before_next_NDF(completed_measurements, current_NDF_idx)
%
% Description:
%   This helper scans the completion tensor beyond the current NDF index
%   and reports whether any future measurement remains unfinished so the
%   caller can decide whether to enforce CombiLED cooldown before
%   continuing.
% Inputs:
%   completed_measurements   - Logical completion tensor over NDF,
%                              contrast, frequency, and repeat
%                              dimensions.
%   current_NDF_idx          - Scalar index of the NDF that just finished.
%
% Outputs:
%   tf                       - Logical scalar indicating whether later NDFs
%                              still have unfinished measurements.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    tf = false;
    for idx = (current_NDF_idx + 1):size(completed_measurements, 1)
        completed_measurements_NDF = completed_measurements(idx, :, :, :);
        if(any(completed_measurements_NDF(:) == false))
            tf = true;
            return;
        end
    end
end

function value = py_module_attr(module, attr_name)
% Fetch an attribute from a Python module or object.
%
% Syntax:
%   value = py_module_attr(module, attr_name)
%
% Description:
%   This helper centralizes `py.getattr` access so the surrounding
%   calibration code can call Python-side functionality without repeating
%   lookup boilerplate.
% Inputs:
%   module                   - Python module or object to inspect.
%   attr_name                - Name of the attribute to retrieve.
%
% Outputs:
%   value                    - Retrieved Python attribute.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    value = py.getattr(module, attr_name);
end

function value = py_call_module_attr(module, attr_name, varargin)
% Fetch and immediately call a Python attribute.
%
% Syntax:
%   value = py_call_module_attr(module, attr_name, varargin)
%
% Description:
%   This helper looks up a named Python callable and invokes it with the
%   supplied arguments, keeping the MATLAB-side Bluetooth code concise.
% Inputs:
%   module                   - Python module or object to inspect.
%   attr_name                - Name of the callable attribute to invoke.
%   varargin                 - Arguments forwarded to the Python callable.
%
% Outputs:
%   value                    - Return value from the Python callable.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    callable = py_module_attr(module, attr_name);
    value = callable(varargin{:});
end

function lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client)
% Read light-logger state with retry handling for transient BLE failures.
%
% Syntax:
%   lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client)
%
% Description:
%   Bluetooth reads can occasionally fail even though the device has
%   already accepted the recording command. This helper retries a limited
%   number of times before rethrowing the final exception.
% Inputs:
%   bluetooth_central        - Python Bluetooth helper module.
%   device_num               - Peripheral identifier used by the wrapper.
%   bluetooth_client         - Persistent BLE client handle.
%
% Outputs:
%   lightlogger_state        - MATLAB struct copy of the peripheral state.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    max_attempts = 5;
    retry_pause_seconds = 1.0;

    for attempt = 1:max_attempts
        try
            lightlogger_state = struct(py_call_module_attr(bluetooth_central, "read_peripheral_matlab_wrapper", device_num, bluetooth_client));
            return;
        catch ME
            if(attempt == max_attempts)
                rethrow(ME);
            end

            fprintf(2, "Transient bluetooth read failure (attempt %d/%d): %s\nRetrying in %.1f s...\n", ...
                    attempt, max_attempts, ME.message, retry_pause_seconds);
            pause(retry_pause_seconds);
        end
    end
end

function report_bluetooth_failure(measurement_name, NDF, contrast_idx, freq_idx, measurement_num, ME)
% Print a structured summary of a Bluetooth failure to stderr.
%
% Syntax:
%   report_bluetooth_failure(measurement_name, NDF, contrast_idx, freq_idx, measurement_num, ME)
%
% Description:
%   This helper formats the active measurement context together with the
%   full MATLAB exception report so failures can be diagnosed from logs.
% Inputs:
%   measurement_name         - Label for the active calibration stage.
%   NDF                      - NDF value active when the failure occurred.
%   contrast_idx             - Contrast index active at failure.
%   freq_idx                 - Frequency index active at failure.
%   measurement_num          - Repeat index active at failure.
%   ME                       - MATLAB exception that was thrown.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    fprintf(2, "%s bluetooth failure at NDF %.3f contrast %d frequency %d measurement %d.\n", measurement_name, NDF, contrast_idx, freq_idx, measurement_num);
    fprintf(2, "%s\n", getReport(ME, "extended", "hyperlinks", "off"));
end

function message = describe_lightlogger_error(lightlogger_state)
% Build a readable error report from a light-logger error state payload.
%
% Syntax:
%   message = describe_lightlogger_error(lightlogger_state)
%
% Description:
%   This helper extracts high-level state fields, nested write-error
%   metadata, and a raw struct dump into one newline-delimited diagnostic
%   string that can be saved inside the calibration metadata.
% Inputs:
%   lightlogger_state        - Peripheral state struct returned by the
%                              Bluetooth wrapper.
%
% Outputs:
%   message                  - Multi-line diagnostic string summarizing the
%                              error state.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    message_parts = strings(0, 1);
    message_parts(end + 1) = "Light logger entered error state.";
    if(isfield(lightlogger_state, "state"))
        message_parts(end + 1) = "state: " + string(lightlogger_state.state);
    end

    if(isfield(lightlogger_state, "error_message"))
        error_message = string(lightlogger_state.error_message);
        if(strlength(strtrim(error_message)) > 0)
            message_parts(end + 1) = "error_message: " + error_message;
        end
    end

    if(isfield(lightlogger_state, "info"))
        info_struct = to_plain_struct(lightlogger_state.info);
        if(isstruct(info_struct) && isfield(info_struct, "write_error"))
            write_error = to_plain_struct(info_struct.write_error);
            if(isstruct(write_error))
                if(isfield(write_error, "message"))
                    value = string(write_error.message);
                    if(strlength(strtrim(value)) > 0)
                        message_parts(end + 1) = "write_error.message: " + value;
                    end
                end
                if(isfield(write_error, "exception_type"))
                    value = string(write_error.exception_type);
                    if(strlength(strtrim(value)) > 0)
                        message_parts(end + 1) = "write_error.exception_type: " + value;
                    end
                end
                if(isfield(write_error, "exception_message"))
                    value = string(write_error.exception_message);
                    if(strlength(strtrim(value)) > 0)
                        message_parts(end + 1) = "write_error.exception_message: " + value;
                    end
                end
                if(isfield(write_error, "traceback"))
                    value = string(write_error.traceback);
                    if(strlength(strtrim(value)) > 0)
                        message_parts(end + 1) = "write_error.traceback: " + value;
                    end
                end
                if(isfield(write_error, "header"))
                    header_struct = to_plain_struct(write_error.header);
                    header_lines = flatten_struct_fields(header_struct, "write_error.header");
                    if(~isempty(header_lines))
                        message_parts = [message_parts; header_lines(:)];
                    end
                end
            end
        end
    end

    raw_state_dump = strtrim(string(evalc("disp(lightlogger_state)")));
    if(strlength(raw_state_dump) > 0)
        message_parts(end + 1) = "raw_state_dump:";
        message_parts(end + 1) = raw_state_dump;
    end

    message = strjoin(message_parts, newline);
end

function value = to_plain_struct(value)
% Normalize Python dict-like values into plain MATLAB structs when possible.
%
% Syntax:
%   value = to_plain_struct(value)
%
% Description:
%   This helper converts `py.dict` objects into MATLAB structs and treats
%   `py.None` as an empty struct so downstream formatting code can inspect
%   nested error payloads with ordinary MATLAB struct logic.
% Inputs:
%   value                    - Struct-like, dict-like, or `py.None` value
%                              to normalize.
%
% Outputs:
%   value                    - Normalized MATLAB struct or original value.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    if(isstruct(value))
        return;
    end

    if(isa(value, "py.dict"))
        value = struct(value);
        return;
    end

    if(isa(value, "py.NoneType"))
        value = struct;
    end
end

function lines = flatten_struct_fields(input_struct, prefix)
% Flatten nested struct contents into prefixed text lines.
%
% Syntax:
%   lines = flatten_struct_fields(input_struct, prefix)
%
% Description:
%   This helper recursively walks a struct and formats each leaf value as a
%   `prefix.field = value` style string.
% Inputs:
%   input_struct             - Struct whose contents should be flattened.
%   prefix                   - Prefix to prepend to each flattened field
%                              name.
%
% Outputs:
%   lines                    - String array of flattened `field=value`
%                              diagnostics.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    lines = strings(0, 1);
    if(~isstruct(input_struct))
        return;
    end

    field_names = fieldnames(input_struct);
    for idx = 1:numel(field_names)
        field_name = field_names{idx};
        field_value = input_struct.(field_name);
        field_prefix = prefix + "." + string(field_name);

        if(isstruct(field_value))
            nested_lines = flatten_struct_fields(field_value, field_prefix);
            if(~isempty(nested_lines))
                lines = [lines; nested_lines(:)];
            end
        else
            lines(end + 1) = field_prefix + ": " + string(local_value_to_text(field_value));
        end
    end
end

function text = local_value_to_text(value)
% Convert one diagnostic value into a printable string.
%
% Syntax:
%   text = local_value_to_text(value)
%
% Description:
%   This helper converts nested diagnostic values into a concise string
%   representation so error payloads can be serialized into readable text
%   reports.
% Inputs:
%   value                    - Value to convert into text.
%
% Outputs:
%   text                     - String scalar representation of the value.
%
% Examples:
%{
    % See collect_temporal_sensitivity_data.m for usage context.
%}

    if(isstring(value) || ischar(value))
        text = string(value);
        return;
    end

    if(islogical(value) || isnumeric(value))
        if(isscalar(value))
            text = string(value);
        else
            text = strtrim(string(mat2str(value)));
        end
        return;
    end

    text = strtrim(string(evalc("disp(value)")));
end
