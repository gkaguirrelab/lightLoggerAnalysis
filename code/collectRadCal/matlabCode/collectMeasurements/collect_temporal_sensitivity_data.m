function [success, temporal_sensitivity_calibration_metadata] = collect_temporal_sensitivity_data(device_num, ...
                                                                                                  temporal_sensitivity_calibration_metadata,...
                                                                                                  bluetooth_central,...
                                                                                                  label,...
                                                                                                  dropbox_savedir,...
                                                                                                  local_savedir...
                                                                                                ) 
% Collect calibration data for temporal sensitivity of 
% all of the sensors in the light logger
%
% Syntax:
%   success = collect_temporal_sensitivity_data(label, dropbox_savedir, CombiLED, modResult, contrast_levels, frequencies, n_seconds, bluetooth_central, sensors)
%
% Description:
%   Collect measurements on the tempporal sensitivity all of the sensors in
%   the light logger by exposing it, connected via the light logger's
%   Bluetooth functionality, to modulations at various frequencies and
%   contrast levels and recording short videos for each. These videos are
%   then uploaded to DropBox.
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
%   modResult             - Struct. Contains information about the
%                           modulation that will be displayed  
%  
%   contrast_levels       - Vector. Represents the list of contrast levels
%                           at which we will modulate 
%
%   frequencies           - Vector. Represents the list of frequencies 
%                           of modulation we will epose to the light logger 
%
%   n_seconds             - Int. The duration of the recording at each
%                           settings level
%
%   bluetooth_central     - Object. Represents a Python module (bluetooth_control.py) 
%                           with functions for bluetooth communication. 
%
%   sensors               - Struct. Represents the sensors to record with 
%                           and their settings (gain, exposure, agc)
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
    modResult = designModulation('LightFlux', photoreceptors, cal); % You will need to load in the cal object and construct the photorecptors object. This is just for the sake of example
    contrast_levels = [0.5, 0,6];
    frequencies = [1, 2, 3, 4, 5];
    n_seconds = 12;
    bluetooth_central = import_pyfile("bluetooth_central.py") % You will need to give the full path to this file
    sensors.M = struct; % You can also initialize more with the first initial of the sensor. Then in the substruct, specificy settings
    sensors.W.agc = false; 
    success = collect_temporal_sensitivity_data(label, dropbox_savedir, CombiLED, modResult, contrast_levels, frequencies, n_seconds, bluetooth_central, sensors)
%}

    % Validate the arguments 
    arguments 
        device_num; 
        temporal_sensitivity_calibration_metadata; 
        bluetooth_central; 
        label = ""; 
        dropbox_savedir = ""; 
        local_savedir = ""; 
    end 

    % First, extract some information from the calibration struct 
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
        completed_measurements_NDF = temporal_sensitivity_calibration_metadata.completed_measurements(nn, :, :);
        if(all(completed_measurements_NDF(:) == true))
            continue 
        end 

        % Initialize the CombiLED for this NDF level 
        % Ensure we have constants from lightlogger
        tbUseProject('lightLogger'); 
        [CombiLED, cal] = initialize_combiLED(getpref("lightLogger", "combiExperiments_path"), NDF); 

        % Save the cal file for this NDF 
        cal_files = temporal_sensitivity_calibration_metadata.cal_files;
        cal_files{nn} = cal; 
        temporal_sensitivity_calibration_metadata.cal_files = cal_files; 

        % Construct the modulation to dispaly
        recording_seconds = 12; % The duration to record for
        observerAgeInYears = 30; % Build general information about the modulation on the CombiLED 
        pupilDiameterMm = 3; 

        photoreceptors = photoreceptorDictionaryHuman('observerAgeInYears',observerAgeInYears,'pupilDiameterMm',pupilDiameterMm);
        modResult = designModulation('LightFlux', photoreceptors, cal); % Construct the modulation 
        modResult.settingsHigh = (background)'; % Set the color profile of the modulation 
        modResult.settingsBackground = (background./2)';
        modResult.settingsLow = (background.*0)';

        % Save the modulation for this NDF
        modulations = temporal_sensitivity_calibration_metadata.modulations;
        modulations{nn} = cal; 
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
                
                % For each contrast level, 
                % randomize the frequency order  
                frequencies_order = squeeze(frequencies_orders(nn, mm, cc, :));

                % Put the frequencies in this order 
                shuffled_frequencies = frequencies(frequencies_order);    

                assert(numel(shuffled_frequencies) == numel(frequencies) ) ; 

                % Iterate over the shuffled frequencies 
                for ff = 1:numel(shuffled_frequencies)
                    % Retrieve the frequency we will expose 
                    freq = shuffled_frequencies(ff); 

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
                    if(temporal_sensitivity_calibration_metadata.completed_measurements(nn, cc, ff, mm))
                        pause(0.5);
                        continue ; 
                    end 

                    % Start the modulation on the CombiLED 
                    CombiLED.startModulation();
                    pause(0.5);

                    % Generate a message to send to the RPi
                    update_message = bluetooth_central.initialize_update_message(); % Initialize the message dict 
                    
                    filename = label + sprintf("TemporalSensitivity_%dcontrastIdx_%dfreqIdx_%dmeasurementIdx", cc, ff, mm); % Define the filename for the recording (replace. in floats with x)

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

                        sucess = 0; 
                        return; 
                    end 

                    bluetooth_central.generate_calibration_state(update_message, py.str(filename),... % Put it all together in the dict 
                                                                 cloud_output_dir, local_output_dir,...
                                                                 true, py.int(recording_seconds),... 
                                                                 py.int(30), sensors);

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

                    % End the modulation 
                    CombiLED.stopModulation(); 
                    pause(0.5);

                    % Otherwise, mark this measurement as completed
                    temporal_sensitivity_calibration_metadata.completed_measurements(nn, cc, ff, mm) = true; 

                end % Frequencies 

            end % Contrasts

        end % Measures 

        % Close the CombiLED for this NDF 
        CombiLED.serialClose(); 

    end % NDF

    % Exit with success if nothing went awry
    success = 1; 

    return ; 

end 