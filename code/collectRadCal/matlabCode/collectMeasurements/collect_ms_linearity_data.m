function [success, ms_linearity_calibration_metadata] = collect_ms_linearity_data(device_num,...
                                                                             ms_linearity_calibration_metadata,... 
                                                                             bluetooth_central,...
                                                                             label,...
                                                                             dropbox_savedir,...
                                                                             local_savedir...
                                                                            )        
% Collect calibration data for the linearity of the light-sensing 
% chips of the light logger (AS/TS) 
%
% Syntax:
%   success = collect_ms_linearity_data(label, dropbox_savedir, CombiLED, background, settings_scalars, n_seconds, bluetooth_central)        
%
% Description:
%   Collect measurements on the linearity of the sensors in the MS 
%   by exposing it, connected via the light logger's Bluetooth 
%   functionality, to linearly increasing static settings on the CombiLED, 
%   and recording short videos at each step. These videos are then 
%   uploaded to DropBox. 
%
% Inputs:
%   label                 - String. Optional notes to attach
%                           to the front of the saved filenames
%
%   dropbox_savedir       - String. The path to the directory 
%                           on DropBox where the files will be 
%                           output 
%
%   CombiLED              - Object. Serial controller for the CombiLED 
%
%   background            - Vector. Represents the initial settings 
%                           of the CombiLED that we will then scale. 
%
%   settings_scalars      - Vector. Represents the linearly spaced 
%                           scalars in order that we will loop over 
%                           and apply to the background. 
%
%   n_seconds             - Int. The duration of the recording at each
%                           settings level
%
%   bluetooth_central     - Object. Represents a Python module (bluetooth_control.py) 
%                           with functions for bluetooth communication. 
%
% Outputs:
%
%   success               - Int. Returns 1 on success, 0 on failure.                 
%
% Examples:
%{
    label = "";
    dropbox_save_dir = "/FLIC_data/LightLogger_RadCal"; 
    CombiLED = combiLEDControl(); % You need to do further initialization on this (selecting cal file), this is just for the sake of example 
    background = ones(1, 8);
    settings_scalars = linspace(0, 0.95, 9);
    n_seconds = 12;
    bluetooth_central = import_pyfile("bluetooth_central.py") % You will need to give the full path to this file
    success = collect_ms_linearity_data(label, dropbox_savedir, CombiLED, background, settings_scalars, n_seconds, bluetooth_central)  
%}

    % Validate the arguments 
    arguments 
        device_num; 
        ms_linearity_calibration_metadata; 
        bluetooth_central;  
        label = ""; 
        dropbox_savedir = ""; 
        local_savedir = ""; 
    end 

    % Extract some information from the MS_linearity calibration struct 
    NDFs = ms_linearity_calibration_metadata.NDFs; 
    background = ms_linearity_calibration_metadata.background; 
    settings_scalars = ms_linearity_calibration_metadata.background_scalars; 
    settings_scalars_orders = ms_linearity_calibration_metadata.background_scalars_orders;
    n_seconds = ms_linearity_calibration_metadata.recording_seconds;
    n_measures = ms_linearity_calibration_metadata.n_measures; 

    % Iterate ovre the NDF levels
    for nn = 1:numel(NDFs); 
        % Retrieve the current NDF level 
        NDF = NDFs(nn); 

        % Determine if there are any measurements to be collected 
        % If there are not, we can simply skip this NDF 
        completed_measurements_NDF = ms_linearity_calibration_metadata.completed_measurements(nn, :, :);
        if(all(completed_measurements_NDF(:) == true))
            continue 
        end 

        % Initialize the CombiLED for this NDF level 
        tbUseProject('lightLogger');
        [CombiLED, cal] = initialize_combiLED(getpref("lightLogger", "combiExperiments_path"), NDF); 

        % Save the cal file for this NDF 
        cal_files = ms_linearity_calibration_metadata.cal_files;
        cal_files{nn} = cal; 
        ms_linearity_calibration_metadata.cal_files = cal_files; 

        % Retrieve the sensors and settings for this NDF level 
        sensors = ms_linearity_calibration_metadata.sensors_and_settings{nn};

        % Pause to ensure the NDF is attached 
        fprintf("Preparing to capture. Ensure NDF: %f is attached. Press any key to continue: \n", NDF)
        pause(); 

        % Iterate over the number of measurements to make 
        for mm = 1:n_measures
            % Retrieve the randomized order of the settings 
            % for this measure 
            settings_scalars_order = settings_scalars_orders(nn, mm, :); 

            % Shuffle the array 
            settings_scalars_shuffled = settings_scalars(settings_scalars_order); 

            % Iterate over the setting scalars for the CombiLED 
            for setting_idx = 1:numel(settings_scalars_shuffled)
                fprintf("MS Linearity | NDF (%d/%d) Setting (%d/%d) Scalar: %.3f M: (%d/%d)\n", ...
                        nn, numel(NDFs),...
                        setting_idx, numel(settings_scalars_shuffled), settings_scalars_shuffled(setting_idx),... 
                        mm, n_measures...
                    );
                
                % Determine if we already visited this combination of settings and measurement. If we did 
                % skip it
                if(ms_linearity_calibration_metadata.completed_measurements(nn, setting_idx, mm))
                    continue; 
                end 

                % Retrieve the scalar for this settings level 
                settings_scalar = settings_scalars_shuffled(setting_idx);
                
                % Calculate the settings by applying the scalar to 
                % the background state
                settings = background * settings_scalar; 
                
                % Apply the settings to the CombiLED 
                CombiLED.setPrimaries(settings);
                
                % Generate a message to send to the RPi
                update_message = bluetooth_central.initialize_update_message();

                filename = label + sprintf("MSlinearity_%dsettingsIdx_%dmeasurementIdx", setting_idx, mm);
                
                % Compose the target output dir for this NDF level
                cloud_output_dir = ""; 
                if(dropbox_savedir ~= "")
                    cloud_output_dir = fullfile(dropbox_savedir, sprintf("NDF%f", NDF));
                end 
                
                local_output_dir = "";
                if(local_savedir ~= "")
                    local_output_dir = fullfile(local_savedir, sprintf("NDF%f", NDF));
                end 
                
                
                bluetooth_central.generate_calibration_state(update_message, py.str(filename),... % Put it all together in the dict 
                                                             cloud_output_dir, local_output_dir,...
                                                             true, py.int(n_seconds),... 
                                                             py.int(30), sensors...
                                                            );

                % Send a message to the RPi to make a recording 
                bluetooth_central.message_peripheral_matlab_wrapper(device_num, update_message);

                % Read the state from the light logger, and pause until it is changed 
                % from calibration to wait 
                % if it changes to error, then throw an error
                while(true)
                    % Read the current state from the light logger 
                    lightlogger_state = struct(bluetooth_central.read_peripheral_matlab_wrapper(device_num));
                    state_name = string(char(lightlogger_state.state));     

                    % If the state is error, raise an error on this machine 
                    % as well 
                    if(state_name == "error") 
                        success = 0 ; 
                        return ; 
                    end 

                    % If the light logger is done recording + uploading 
                    % break from the loop
                    if(state_name == "wait")
                        break 
                    end
                    
                    % Pause for some time between reads
                    pause(2); 
                end 

                % Otherwise, we succesfully completed and we will mark this 
                % measurement as done 
                ms_linearity_calibration_metadata.completed_measurements(nn, setting_idx, mm) = true; 
            
            end % Settings  

        end % Measures 
    
        % Close the CombiLED for this NDF level 
        CombiLED.serialClose(); 

    end % NDFs 

    % Set the success flag and return 
    success = 1; 

end 