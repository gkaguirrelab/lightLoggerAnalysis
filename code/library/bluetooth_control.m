function bluetooth_control()

    % First, get the current directory the user is in so that we can 
    % return to it after importing the Python modules 
    current_dir = pwd();

    % Retrieve the path to the phone firmware directory
    [lightlogger_dir, ~, ~] = fileparts(fileparts(mfilename("fullpath"))); 
    phone_firmware_dir = fullfile(lightlogger_dir, "phone_firmware"); 

    % Import the Python utility module from that path 
    cd(phone_firmware_dir); 
    bluetooth_central = py.importlib.import_module('bluetooth_central'); 

    % After the module has been imported, we can return to the directory 
    % where we started 
    cd(current_dir);

    % Initialize a message that we will send to the RPi for bluetooth control 
    update_message = bluetooth_central.initialize_update_message();

    % Prompt the user to select a state to generate
    states = {"wait", "science", "calibration", "upload", "snapshot"};
    disp("Enter a state to generate")
    for ss = 1:numel(states)
        fprintf("%d | %s\n", ss, states{ss}); 
    end

    state_name = input("Enter a state to generate: ");

    % Generate a state name based on user choice 
    switch(state_name)
        case 1 
            % No prompt is needed for generating the wait state, as it homes the device 
            disp("----GENERATING WAIT STATE----")
            bluetooth_central.generate_wait_state(update_message);

        case 2
            disp("----GENERATING SCIENCE STATE----")
            % Retrieve a name for the recording 
            recording_name = input("Enter a name for the recording (myCalibration): ", "s"); 

            % Determine a chunk duration for the recording 
            chunk_duration = input("Enter a length in seconds for chunk/buffer size (1-30): "); 

            % Select which sensors to use
            sensor_names = {"W", "P", "M"};
            
            disp("Select which sensors to use: ");
            for ss = 1:numel(sensor_names)
                fprintf("%d | %s\n", ss, sensor_names{ss});
            end
            sensors_to_use = num2cell(input("Enter sensors to use (123): ", "s")); 
            
            % Build a sensor struct based on the sensors to use 
            sensors_struct = struct; 
            for cc = 1:numel(sensors_to_use)
                % Retrieve the name of the sensor to initialize parameters for 
                sensor_num = str2double(sensors_to_use{cc});
                sensor_name = sensor_names{sensor_num}; 

                sensors_struct.(sensor_name) = struct; 
                
                % Right now, only the world camera takes parameters
                if(sensor_name ~= "W")
                    continue; 
                end 

                % Otherwise, prompt for world information
                gain = input("Sensor: W | Select an initial gain to use for recording (1.0-10.0): ");
                exposure = input("Sensor: W | Select an initial exposure value to use for recording (37-5000): ");
                use_agc = input("Sensor: W | Select whether or not to use our custom AGC (0-1): ");
                save_agc_metadata = input("Sensor: W | Select whether or not to use our custom AGC (0-1): ");

                sensors_struct.(sensor_name).gain = gain; 
                sensors_struct.(sensor_name).exposure = exposure;
                sensors_struct.(sensor_name).use_agc = use_agc;
                sensors_struct.(sensor_name).save_agc_metadata = save_agc_metadata;
            end 

            % Generate the state from the user input 
            bluetooth_central.generate_science_state(update_message, recording_name, chunk_duration, sensors_struct);

        case 3
            disp("----GENERATING CALIBRATION STATE----")
            
            % Retrieve a name for the recording 
            recording_name = input("Enter a name for the recording (myCalibration): ", "s"); 

            % Retrieve a length in seconds for recording 
            total_duration = input("Enter a length in seconds for the total duration of recording: ");

            % Determine a chunk duration for the recording 
            chunk_duration = input("Enter a length in seconds for chunk/buffer size (1-30): "); 

            % Determine whether or not to delete videos after uploading 
            delete_recording = input("Delete the recording after upload? (0/1)"); 

            % Select which sensors to use
            sensor_names = {"W", "P", "M"};
            
            disp("Select which sensors to use: ");
            for ss = 1:numel(sensor_names)
                fprintf("%d | %s\n", ss, sensor_names{ss});
            end
            sensors_to_use = num2cell(input("Enter sensors to use (123): ", "s")); 
            
            % Build a sensor struct based on the sensors to use 
            sensors_struct = struct; 
            for cc = 1:numel(sensors_to_use)
                % Retrieve the name of the sensor to initialize parameters for 
                sensor_num = str2double(sensors_to_use{cc});
                sensor_name = sensor_names{sensor_num}; 

                sensors_struct.(sensor_name) = struct; 
                
                % Right now, only the world camera takes parameters
                if(sensor_name ~= "W")
                    continue; 
                end 

                % Otherwise, prompt for world information
                gain = input("Sensor: W | Select an initial gain to use for recording (1.0-10.0): ");
                exposure = input("Sensor: W | Select an initial exposure value to use for recording (37-5000): ");
                use_agc = input("Sensor: W | Select whether or not to use our custom AGC (0-1): ");
                save_agc_metadata = input("Sensor: W | Select whether or not to use our custom AGC (0-1): ");

                sensors_struct.(sensor_name).gain = gain; 
                sensors_struct.(sensor_name).exposure = exposure;
                sensors_struct.(sensor_name).use_agc = use_agc;
                sensors_struct.(sensor_name).save_agc_metadata = save_agc_metadata;
            end 
            
            % Generate the state from the user input 
            bluetooth_central.generate_science_state(update_message, recording_name, chunk_duration, sensors_struct);
        
        case 4
            disp("----GENERATING UPLOAD STATE----")

            % Retrieve a name for the recording to upload 
            recording_name = input("Enter a name for the recording (myUpload): ", "s"); 

            % Determine whether to delete the video after upload 
            delete_recording = input("Delete the recording after upload? (0/1)"); 
            
            % Generate the state from the user input 
            bluetooth_central.generate_upload_state(update_message, recording_name, delete_recording);

        case 5
            disp("----GENERATING SNAPSHOT STATE----")

            % Select which sensors to use
            sensor_names = {"W", "P"};
        
            disp("Select which sensors to use: ");
            for ss = 1:numel(sensor_names)
                fprintf("%d | %s\n", ss, sensor_names{ss});
            end
            sensors_to_use = num2cell(input("Enter sensors to use (12): ", "s")); 
            
            % Build a sensor struct based on the sensors to use 
            sensors_struct = struct; 
            for cc = 1:numel(sensors_to_use)
                % Retrieve the name of the sensor to initialize parameters for 
                sensor_num = str2double(sensors_to_use{cc});
                sensor_name = sensor_names{sensor_num}; 

                sensors_struct.(sensor_name) = struct; 
                
                % Right now, only the world camera takes parameters
                if(sensor_name ~= "W")
                    continue; 
                end 

                % Otherwise, prompt for world information
                gain = input("Sensor: W | Select an initial gain to use for recording (1.0-10.0): ");
                exposure = input("Sensor: W | Select an initial exposure value to use for recording (37-5000): ");
                use_agc = input("Sensor: W | Select whether or not to use our custom AGC (0-1): ");

                sensors_struct.(sensor_name).gain = gain; 
                sensors_struct.(sensor_name).exposure = exposure;
                sensors_struct.(sensor_name).use_agc = use_agc;
            end 
            
            % Generate the state from the user input 
            bluetooth_central.generate_snapshot_state(update_message, sensors_struct);

    end 

    disp("You have generated the following state: "); 
    disp(update_message);
    disp("Press any key to upload to device."); 
    pause() 

    % Send a message to the RPi 
    bluetooth_central.message_peripheral_matlab_wrapper(update_message);

    % Visualize the snapshots on the local machine if we sent a snapshot state message
    if(state_name == 5)
        bluetooth_central.visualize_snapshots(py.None)
    end 

end 