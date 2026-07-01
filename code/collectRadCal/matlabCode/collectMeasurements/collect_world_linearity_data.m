function [success, world_linearity_calibration_metadata] = collect_world_linearity_data(device_num,...
                                                                                         world_linearity_calibration_metadata,...
                                                                                         bluetooth_central,...
                                                                                         bluetooth_client,...
                                                                                         label,...
                                                                                         dropbox_savedir,...
                                                                                         local_savedir,...
                                                                                         cooldown_callback...
                                                                                        )
% Collect calibration data for the linearity of the world camera.
%
% Syntax:
%   [success, world_linearity_calibration_metadata] = collect_world_linearity_data(device_num, world_linearity_calibration_metadata, bluetooth_central, bluetooth_client, label, dropbox_savedir, local_savedir, cooldown_callback)
%
% Description:
%   This function coordinates the world-camera linearity calibration by
%   iterating over NDF levels, fixed AGC-target-derived camera settings,
%   randomized CombiLED intensity scalars, and repeated measurements. For
%   each requested condition it initializes the CombiLED for the active
%   NDF, applies the scaled background spectrum, triggers a light-logger
%   calibration recording over Bluetooth, waits until the peripheral
%   returns to its idle state, and marks the condition as completed in the
%   metadata struct so interrupted runs can be resumed.
%
% Inputs:
%   device_num                           - Int. The device id of the RPi
%                                          that will run the calibration.
%
%   world_linearity_calibration_metadata - Struct. Metadata and progress
%                                          state for the world camera
%                                          linearity calibration.
%
%   bluetooth_central                    - Object. Represents a Python
%                                          module (bluetooth_control.py)
%                                          with functions for bluetooth
%                                          communication.
%
%   label                                - String. Optional notes to
%                                          attach to the front of the
%                                          saved filenames.
%
%   dropbox_savedir                      - String. The path to the
%                                          directory on DropBox where the
%                                          files will be output.
%
%   local_savedir                        - String. The path to the local
%                                          directory where files will be
%                                          saved.
%
% Outputs:
%
%   success                              - Logical scalar. Returns true on
%                                          success and false on failure.
%   world_linearity_calibration_metadata - Updated metadata struct with
%                                          completion flags, calibration
%                                          files, and any saved error
%                                          message.
%

    % Validate the arguments
    arguments
        device_num;
        world_linearity_calibration_metadata;
        bluetooth_central;
        bluetooth_client;
        label = "";
        dropbox_savedir = "";
        local_savedir = "";
        cooldown_callback = [];
    end

    % Extract some information from the world_linearity calibration struct
    world_linearity_calibration_metadata.last_error_message = "";
    NDFs = world_linearity_calibration_metadata.NDFs;
    contrast_agc_targets = world_linearity_calibration_metadata.contrast_agc_targets; 
    background = world_linearity_calibration_metadata.background;
    settings_scalars = world_linearity_calibration_metadata.background_scalars;
    settings_scalars_orders = world_linearity_calibration_metadata.background_scalars_orders;
    n_seconds = world_linearity_calibration_metadata.recording_seconds;
    n_measures = world_linearity_calibration_metadata.n_measures;

    % Iterate over the NDF levels
    for nn = 1:numel(NDFs)
        % Retrieve the current NDF level
        NDF = NDFs(nn);

        % Determine if there are any measurements to be collected.
        % If there are not, we can simply skip this NDF.
        completed_measurements_NDF = world_linearity_calibration_metadata.completed_measurements(nn, :, :, :);
        if(all(completed_measurements_NDF(:) == true))
            continue
        end

        % Initialize the CombiLED for this NDF level
        tbUseProject('lightLogger');
        [CombiLED, cal] = initialize_combiLED_lightLogger(getpref("lightLogger", "combiExperiments_path"), round(NDF));

        % Save the cal file for this NDF
        cal_files = world_linearity_calibration_metadata.cal_files;
        cal_files{nn} = cal;
        world_linearity_calibration_metadata.cal_files = cal_files;

        % Pause to ensure the NDF is attached
        fprintf("Preparing to capture world camera linearity. Ensure NDF: %f is attached. Press any key to continue: \n", NDF)
        pause();

        % Iterate over the contrast AGC targets for this NDF levels
        for cc = 1:numel(contrast_agc_targets)
            % Retrieve the current contrast AGC target
            contrast_agc_target = contrast_agc_targets(cc); 
            contrast_target_fieldname = "x" + strrep(string(contrast_agc_target), ".", "x"); 
            contrast_target_sensors = world_linearity_calibration_metadata.sensors_and_settings.(contrast_target_fieldname); 

            % Retrieve the sensors and settings for this NDF level
            sensors = contrast_target_sensors{nn};

            % Iterate over the number of measurements to make
            for mm = 1:n_measures
                % Retrieve the randomized order of the settings for this measure
                settings_scalars_order = squeeze(settings_scalars_orders(nn, mm, :));

                % Shuffle the array
                settings_scalars_shuffled = settings_scalars(settings_scalars_order);

                % Iterate over the setting scalars for the CombiLED
                for setting_order_idx = 1:numel(settings_scalars_shuffled)
                    settings_scalar_idx = settings_scalars_order(setting_order_idx);
                    fprintf("World Camera Linearity | NDF (%d/%d) Contrast AGC target (%d/%d) Setting (%d/%d) Scalar: %.3f M: (%d/%d)\n", ...
                            nn, numel(NDFs), cc, numel(contrast_agc_targets),...
                            setting_order_idx, numel(settings_scalars_shuffled), settings_scalars_shuffled(setting_order_idx),...
                            mm, n_measures...
                        );

                    % Determine if we already visited this combination of
                    % settings and measurement. If we did, skip it.
                    if(world_linearity_calibration_metadata.completed_measurements(nn, cc, settings_scalar_idx, mm))
                        continue;
                    end

                    % Retrieve the scalar for this settings level
                    settings_scalar = settings_scalars_shuffled(setting_order_idx);

                    % Calculate the settings by applying the scalar to the
                    % background state
                    settings = background * settings_scalar;

                    % Apply the settings to the CombiLED
                    CombiLED.setPrimaries(settings);

                    try
                        % Generate a message to send to the RPi
                        update_message = py_call_module_attr(bluetooth_central, "initialize_update_message");

                        filename = label + sprintf("WorldCameraLinearity_%dcontrastTargetIdx_%dsettingsIdx_%dmeasurementIdx", ...
                                                   cc, settings_scalar_idx, mm);

                        % Compose the target output dir for this NDF level
                        cloud_output_dir = "";
                        if(dropbox_savedir ~= "")
                            cloud_output_dir = fullfile(dropbox_savedir, sprintf("NDF%f", NDF));
                        end

                        local_output_dir = "";
                        if(local_savedir ~= "")
                            local_output_dir = fullfile(local_savedir, sprintf("NDF%f", NDF));
                        end

                        py_call_module_attr(bluetooth_central, "generate_calibration_state", update_message, py.str(filename),...
                                            cloud_output_dir, local_output_dir,...
                                            true, py.int(n_seconds),...
                                            py.int(30), sensors...
                                            );

                        % Send a message to the RPi to make a recording
                        py_call_module_attr(bluetooth_central, "message_peripheral_matlab_wrapper", device_num, update_message, bluetooth_client);

                        % Read the state from the light logger, and pause until it
                        % is changed from calibration to wait. If it changes to
                        % error, then throw an error.
                        while(true)
                            % Read the current state from the light logger.
                            % Transient CoreBluetooth read errors (e.g. CBATT
                            % "Unlikely error", code 14) are retried rather than
                            % aborting the whole calibration, since the recording
                            % has already been triggered on the light logger.
                            lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client);
                            state_name = string(char(lightlogger_state.state));

                            % If the state is error, raise an error on this
                            % machine as well
                            if(state_name == "error")
                                world_linearity_calibration_metadata.last_error_message = describe_lightlogger_error(lightlogger_state);
                                success = 0;
                                return;
                            end

                            % If the light logger is done recording + uploading,
                            % break from the loop
                            if(state_name == "wait")
                                break
                            end

                            % Pause for some time between reads
                            pause(0.5);
                        end
                    catch ME
                        report_bluetooth_failure("World linearity", NDF, cc, mm, ME);
                        world_linearity_calibration_metadata.last_error_message = ...
                            sprintf("World linearity bluetooth failure at NDF %.3f contrast %d measurement %d.\n%s", ...
                                    NDF, cc, mm, getReport(ME, "extended", "hyperlinks", "off"));
                        success = 0;
                        return;
                    end

                    % Otherwise, we successfully completed and we will mark
                    % this measurement as done
                    world_linearity_calibration_metadata.completed_measurements(nn, cc, settings_scalar_idx, mm) = true;

                end % Settings

            end % Measures
        end % Contrast AGC targets

        if(should_cooldown_before_next_NDF(world_linearity_calibration_metadata.completed_measurements, nn) && ~isempty(cooldown_callback))
            cooldown_callback(CombiLED, NDF);
        end

        % Close the CombiLED for this NDF level
        CombiLED.serialClose();

    end % NDFs

    % Set the success flag and return
    success = 1;

end

function tf = should_cooldown_before_next_NDF(completed_measurements, current_NDF_idx)
% Decide whether another unfinished NDF remains after the current one.
%
% Syntax:
%   tf = should_cooldown_before_next_NDF(completed_measurements, current_NDF_idx)
%
% Description:
%   This helper scans the completion tensor beyond the current NDF index
%   and reports whether any future measurement remains unfinished. The
%   caller uses this to decide whether to invoke the CombiLED cooldown
%   routine before advancing.
% Inputs:
%   completed_measurements   - Logical completion tensor over NDF,
%                              contrast target, setting, and repeat
%                              dimensions.
%   current_NDF_idx          - Scalar index of the NDF that just finished.
%
% Outputs:
%   tf                       - Logical scalar indicating whether later NDFs
%                              still have unfinished measurements.
%
% Examples:
%{
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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

function report_bluetooth_failure(measurement_name, NDF, contrast_idx, measurement_num, ME)
% Print a structured summary of a Bluetooth failure to stderr.
%
% Syntax:
%   report_bluetooth_failure(measurement_name, NDF, contrast_idx, measurement_num, ME)
%
% Description:
%   This helper formats the active measurement context together with the
%   full MATLAB exception report so failures can be diagnosed from logs.
% Inputs:
%   measurement_name         - Label for the active calibration stage.
%   NDF                      - NDF value active when the failure occurred.
%   contrast_idx             - Contrast-target index active at failure.
%   measurement_num          - Repeat index active at failure.
%   ME                       - MATLAB exception that was thrown.
%
% Outputs:
%   None.
%
% Examples:
%{
    % See collect_world_linearity_data.m for usage context.
%}

    fprintf(2, "%s bluetooth failure at NDF %.3f contrast %d measurement %d.\n", measurement_name, NDF, contrast_idx, measurement_num);
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
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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
    % See collect_world_linearity_data.m for usage context.
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
