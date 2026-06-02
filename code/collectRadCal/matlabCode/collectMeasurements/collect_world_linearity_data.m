function [success, world_linearity_calibration_metadata] = collect_world_linearity_data(device_num,...
                                                                                         world_linearity_calibration_metadata,...
                                                                                         bluetooth_central,...
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
        label = "";
        dropbox_savedir = "";
        local_savedir = "";
    end

    % Extract some information from the world_linearity calibration struct
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

                    % Generate a message to send to the RPi
                    update_message = bluetooth_central.initialize_update_message();

                    filename = label + sprintf("WorldCameraLinearity_%0.2fcontrastTarget_%dsettingsIdx_%dmeasurementIdx", ...
                                               contrast_agc_target, settings_scalar_idx, mm);

                    % Compose the target output dir for this NDF level
                    cloud_output_dir = "";
                    if(dropbox_savedir ~= "")
                        cloud_output_dir = fullfile(dropbox_savedir, sprintf("NDF%f", NDF));
                    end

                    local_output_dir = "";
                    if(local_savedir ~= "")
                        local_output_dir = fullfile(local_savedir, sprintf("NDF%f", NDF));
                    end

                    bluetooth_central.generate_calibration_state(update_message, py.str(filename),...
                                                                cloud_output_dir, local_output_dir,...
                                                                true, py.int(n_seconds),...
                                                                py.int(30), sensors...
                                                                );

                    % Send a message to the RPi to make a recording
                    bluetooth_central.message_peripheral_matlab_wrapper(device_num, update_message);

                    % Read the state from the light logger, and pause until it
                    % is changed from calibration to wait. If it changes to
                    % error, then throw an error.
                    while(true)
                        % Read the current state from the light logger
                        lightlogger_state = struct(bluetooth_central.read_peripheral_matlab_wrapper(device_num));
                        state_name = string(char(lightlogger_state.state));

                        % If the state is error, raise an error on this
                        % machine as well
                        if(state_name == "error")
                            success = 0;
                            return;
                        end

                        % If the light logger is done recording + uploading,
                        % break from the loop
                        if(state_name == "wait")
                            break
                        end

                        % Pause for some time between reads
                        pause(2);
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
