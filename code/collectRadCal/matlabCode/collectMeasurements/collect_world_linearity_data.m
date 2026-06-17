function [success, world_linearity_calibration_metadata] = collect_world_linearity_data(device_num,...
                                                                                         world_linearity_calibration_metadata,...
                                                                                         bluetooth_central,...
                                                                                         bluetooth_client,...
                                                                                         label,...
                                                                                         dropbox_savedir,...
                                                                                         local_savedir...
                                                                                        )
% Collect calibration data for the linearity of the world camera.
%
% Syntax:
%   success = collect_world_linearity_data(device_num, world_linearity_calibration_metadata, bluetooth_central, label, dropbox_savedir, local_savedir)
%
% Description:
%   Collect measurements on the linearity of the world camera by exposing
%   it, connected via the light logger's Bluetooth functionality, to
%   linearly increasing static settings on the CombiLED and recording
%   short videos at each step. These videos are then uploaded to DropBox.
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
%   success                              - Int. Returns 1 on success, 0 on
%                                          failure.
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
        [CombiLED, cal] = initialize_combiLED(getpref("lightLogger", "combiExperiments_path"), NDF);

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

        % Close the CombiLED for this NDF level
        CombiLED.serialClose();

    end % NDFs

    % Set the success flag and return
    success = 1;

end

function value = py_module_attr(module, attr_name)
    value = py.getattr(module, attr_name);
end

function value = py_call_module_attr(module, attr_name, varargin)
    callable = py_module_attr(module, attr_name);
    value = callable(varargin{:});
end

function lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client)
    % Read the peripheral state, retrying on transient bluetooth failures.
    %
    % CoreBluetooth/Bleak occasionally throws transient read errors (for
    % example "Unlikely error", CBATTErrorDomain code 14). These are not
    % indicative of a real failure of the calibration, so we retry the read a
    % few times before giving up and letting the error propagate.
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
    fprintf(2, "%s bluetooth failure at NDF %.3f contrast %d measurement %d.\n", measurement_name, NDF, contrast_idx, measurement_num);
    fprintf(2, "%s\n", getReport(ME, "extended", "hyperlinks", "off"));
end

function message = describe_lightlogger_error(lightlogger_state)
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
