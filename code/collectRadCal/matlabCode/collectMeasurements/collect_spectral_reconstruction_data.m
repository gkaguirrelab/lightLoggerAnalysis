function [success, spectral_reconstruction_calibration_metadata] = collect_spectral_reconstruction_data(device_num,...
                                                                                                         spectral_reconstruction_calibration_metadata,...
                                                                                                         bluetooth_central,...
                                                                                                         bluetooth_client,...
                                                                                                         label,...
                                                                                                         dropbox_savedir,...
                                                                                                         local_savedir,...
                                                                                                         cooldown_callback...
                                                                                                        )
% Collect paired minispect and world-camera spectral reconstruction data.
%
% Syntax:
%   [success, spectral_reconstruction_calibration_metadata] = collect_spectral_reconstruction_data(device_num, spectral_reconstruction_calibration_metadata, bluetooth_central, bluetooth_client, label, dropbox_savedir, local_savedir, cooldown_callback)
%
% Description:
%   This function coordinates spectral-reconstruction calibration by
%   iterating over NDF levels and randomized CombiLED primary vectors.
%   Ten independently randomized eight-primary vectors are used at each
%   NDF, with three recordings made at each vector. All setting/repeat
%   pairs are exposed in a saved randomized order. Each recording includes
%   the minispect and world camera, with fixed world-camera settings and
%   AGC disabled. Completion is recorded in the metadata so interrupted
%   runs can resume without repeating finished measurements.
%
% Inputs:
%   device_num            - Scalar. Identifier of the target light logger
%                           peripheral.
%   spectral_reconstruction_calibration_metadata
%                         - Struct. Calibration plan and progress record
%                           for the spectral-reconstruction stage.
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
%   success               - Logical scalar. Returns true on success and
%                           false on failure.
%   spectral_reconstruction_calibration_metadata
%                         - Updated metadata struct with completion flags,
%                           calibration files, and any saved error
%                           message.

    arguments
        device_num;
        spectral_reconstruction_calibration_metadata;
        bluetooth_central;
        bluetooth_client;
        label = "";
        dropbox_savedir = "";
        local_savedir = "";
        cooldown_callback = [];
    end

    % Extract information from the spectral reconstruction calibration struct.
    spectral_reconstruction_calibration_metadata.last_error_message = "";
    NDFs = spectral_reconstruction_calibration_metadata.NDFs;
    combiLED_settings = spectral_reconstruction_calibration_metadata.combiLED_settings;
    measurement_orders = spectral_reconstruction_calibration_metadata.measurement_orders;
    n_seconds = spectral_reconstruction_calibration_metadata.recording_seconds;
    n_measures = spectral_reconstruction_calibration_metadata.n_measures;

    % Iterate over the NDF levels.
    for nn = 1:numel(NDFs)
        NDF = NDFs(nn);

        % Skip this NDF when all of its measurements are already complete.
        completed_measurements_NDF = spectral_reconstruction_calibration_metadata.completed_measurements(nn, :, :);
        if(all(completed_measurements_NDF(:) == true))
            continue
        end

        % Initialize the CombiLED for this NDF level.
        tbUseProject('lightLogger');
        [CombiLED, cal] = initialize_combiLED_lightLogger(getpref("lightLogger", "combiExperiments_path"), round(NDF));

        % Save the calibration file for this NDF.
        cal_files = spectral_reconstruction_calibration_metadata.cal_files;
        cal_files{nn} = cal;
        spectral_reconstruction_calibration_metadata.cal_files = cal_files;

        % Retrieve the minispect and fixed world-camera settings for this NDF.
        sensors = spectral_reconstruction_calibration_metadata.sensors_and_settings{nn};

        fprintf("Preparing to capture spectral reconstruction data. Ensure NDF: %f is attached. Press any key to continue: \n", NDF)
        pause();

        % Make three passes through the settings. Each pass contains every
        % settings vector exactly once, in a newly randomized order.
        for measurement_idx = 1:n_measures
            settings_order = squeeze(measurement_orders(nn, measurement_idx, :));

            for settings_order_idx = 1:numel(settings_order)
                settings_idx = settings_order(settings_order_idx);

                fprintf("Spectral Reconstruction | NDF (%d/%d) Setting (%d/%d) M: (%d/%d)\n", ...
                        nn, numel(NDFs), ...
                        settings_order_idx, numel(settings_order), ...
                        measurement_idx, n_measures ...
                       );

                if(spectral_reconstruction_calibration_metadata.completed_measurements(nn, settings_idx, measurement_idx))
                    continue;
                end

                settings = squeeze(combiLED_settings(nn, settings_idx, :))';
                fprintf("Spectral Reconstruction | Settings index: %d | CombiLED primaries: %s\n", ...
                        settings_idx, mat2str(settings, 6));
                CombiLED.setPrimaries(settings);

                try
                    update_message = py_call_module_attr(bluetooth_central, "initialize_update_message");

                    filename = label + sprintf("SpectralReconstruction_%dNDFIdx_%0.3fNDF_%dsettingsIdx_%dmeasurementIdx", ...
                                               nn, NDF, settings_idx, measurement_idx);

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

                    py_call_module_attr(bluetooth_central, "message_peripheral_matlab_wrapper", device_num, update_message, bluetooth_client);

                    while(true)
                        lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client);
                        state_name = string(char(lightlogger_state.state));

                        if(state_name == "error")
                            spectral_reconstruction_calibration_metadata.last_error_message = describe_lightlogger_error(lightlogger_state);
                            success = 0;
                            return;
                        end

                        if(state_name == "wait")
                            break
                        end

                        pause(0.5);
                    end
                catch ME
                    report_bluetooth_failure("Spectral reconstruction", NDF, settings_idx, measurement_idx, ME);
                    spectral_reconstruction_calibration_metadata.last_error_message = ...
                        sprintf("Spectral reconstruction bluetooth failure at NDF %.3f setting %d measurement %d.\n%s", ...
                                NDF, settings_idx, measurement_idx, getReport(ME, "extended", "hyperlinks", "off"));
                    success = 0;
                    return;
                end

                spectral_reconstruction_calibration_metadata.completed_measurements(nn, settings_idx, measurement_idx) = true;
            end
        end

        if(should_cooldown_before_next_NDF(spectral_reconstruction_calibration_metadata.completed_measurements, nn) && ~isempty(cooldown_callback))
            cooldown_callback(CombiLED, NDF);
        end

        CombiLED.serialClose();
    end

    success = 1;
end

function tf = should_cooldown_before_next_NDF(completed_measurements, current_NDF_idx)
% Decide whether another unfinished NDF remains after the current one.

    tf = false;
    for idx = (current_NDF_idx + 1):size(completed_measurements, 1)
        completed_measurements_NDF = completed_measurements(idx, :, :);
        if(any(completed_measurements_NDF(:) == false))
            tf = true;
            return;
        end
    end
end

function value = py_module_attr(module, attr_name)
% Fetch an attribute from a Python module or object.

    value = py.getattr(module, attr_name);
end

function value = py_call_module_attr(module, attr_name, varargin)
% Fetch and immediately call a Python attribute.

    callable = py_module_attr(module, attr_name);
    value = callable(varargin{:});
end

function lightlogger_state = read_peripheral_state_with_retry(bluetooth_central, device_num, bluetooth_client)
% Read light-logger state with retry handling for transient BLE failures.

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

function report_bluetooth_failure(measurement_name, NDF, settings_idx, measurement_idx, ME)
% Print a structured summary of a Bluetooth failure to stderr.

    fprintf(2, "%s bluetooth failure at NDF %.3f setting %d measurement %d.\n", ...
            measurement_name, NDF, settings_idx, measurement_idx);
    fprintf(2, "%s\n", getReport(ME, "extended", "hyperlinks", "off"));
end

function message = describe_lightlogger_error(lightlogger_state)
% Build a readable error report from a light-logger error-state payload.

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
                write_error_fields = {"message", "exception_type", "exception_message", "traceback"};
                for field_idx = 1:numel(write_error_fields)
                    field_name = write_error_fields{field_idx};
                    if(isfield(write_error, field_name))
                        value = string(write_error.(field_name));
                        if(strlength(strtrim(value)) > 0)
                            message_parts(end + 1) = "write_error." + field_name + ": " + value;
                        end
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

    if(isa(value, "py.NoneType"))
        value = struct;
        return;
    end

    if(isa(value, "py.dict"))
        value = struct(value);
    end
end
