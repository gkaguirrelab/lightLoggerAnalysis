function collect_light_logger_calibration_data(experiment_name, device_num, sensor_ids, NDFs, calibration_metadata, upload_video_data)
% Collect calibration data for all of the sensors given their ids on a specific light logger. 
%
% Syntax:
%   collect_light_logger_calibration_data(device_num, sensor_ids, calibration_metadata, upload_video_data)
%
% Description:
%   Perform serveral calibration operations between the light logger,
%   connected remotely via Bluetooth, as well as the CombiLED, 
%   connected over the serial port. They are 1. measure the 
%   linearity of the light sensing chips of the MS, 
%   2. measure the linearity of the world camera,
%   3. measure the temporal sensitivity of all sensors,
%   4. measure the offset between sensors in phase/time of all sensors, and
%   5. measure the contrast gamma of the sensors of the world camera. 
%   Upload these either to DropBox or an external SSD. 
%
% Inputs:
%   device_num            - Int. The deviceID of the RPI 
%                           that will run the calibration.
%                           We need this to target the right 
%                           one to connect to over Bluetooth.  
%   
%   sensor_ids            - Vector. The numerical id of 
%                           the sensors in the form 
%                           [W_id, P_id, M_id]
%
%   NDF                   - Numeric. The NDF that is present for
%                           this calibration.
%   
%   calibration_metadata  - Struct. Primarily optional, 
%                           as a new one will be generated 
%                           and saved to DropBox, when this is 
%                           false. Otherwise, simply take the 
%                           incompleted struct and pass it back 
%                           into the function 
%
%   upload_video_data    -  Bool. Whether or not to upload 
%                           videos onto DropBox, as opposed 
%                           to storing them on an external SSD. 
%
% Outputs:
%
%   NONE                
%
% Examples:
%{
    experiment_name = "test"; 
    device_num = 1; 
    sensor_ids = [1, 1, 1]; 
    NDFs.ms_linearity = [0, 1, 2, 3, 4, 5, 6]; 
    NDFs.world_linearity = [];
    NDFs.temporal_sensitivity = [];
    NDFs.phase_fitting = [];
    NDFs.contrast_gamma = [];
    calibration_metadata = false; 
    upload_video_data = false; 
    collect_light_logger_calibration_data(experiment_name, device_num, sensor_ids, NDFs, calibration_metadata, upload_video_data);
%}  


    % Parse and validate the arguments 
    arguments
        experiment_name {mustBeText} % The name of the calibration that will determine the folder it is saved in, e.g. MSOnly 
        device_num {mustBeInteger}; % The device id of the RPI that will run the calibration
        sensor_ids {mustBeVector}; % The unique ids for each of the sensors
        NDFs  % The NDFs to measure for given sensor(s) in struct form
        calibration_metadata = false; % A  calibration metadata data struct as it will be defined later. 
                                      % in this function. This is optional 
                                      % and useful primarily for picking up when a recording 
                                      % has failed due to non-deterministic networking or bluetooth 
                                      % errors.  
        upload_video_data = false; 
        
    end     

    % Terminate pyenv to get it out of bad states and force a restart 
    terminate(pyenv);

    % Ensure the sensor ids is a vector with at most 3 elements 
    % [W_id, P_id, M_id]
    sensor_ids = sensor_ids(:); % Flatten the id vector 
    assert(numel(sensor_ids) > 0 && numel(sensor_ids) <= 3); % Ensure it has the 3 elements

    % Build the experiment name from the sensor ids by interveaving 
    % the sensor id number after the sensor initial
    sensor_labels = {'W', 'P', 'M'};
    sensor_id_parts = strings(1, numel(sensor_ids));
    for ii = 1:numel(sensor_ids)
        sensor_id_parts(ii) = sensor_labels{ii} + string(sensor_ids(ii));
    end
    sensor_ids = strjoin(sensor_id_parts, "");

    % This is a hacky workaround because the project dependencies are not
    % properly sorted out and tbUseProject can drop this folder from path.
    current_dir = fileparts(mfilename('fullpath'));

    disp("Calibration | Importing libraries...");
    tbUseProject('lightLogger');
    addpath(current_dir);

    % Import the Python library used to communicate 
    % via Bluetooth with the peripheral 
    bluetooth_central = import_pyfile_light_logger(getpref("lightLogger", "bluetooth_central_path"));

    % Import the Python library used for uploading to Dropbox 
    upload_mode = import_pyfile_light_logger(getpref("lightLogger", "upload_mode_path"));

    % Import the utility files of the world and pupil cameras 
    tbUseProject('lightLoggerAnalysis'); 
    disp("Importing utility files"); 
    world_util_path = getpref("lightLoggerAnalysis", "world_util_path"); 
    pupil_util_path = getpref("lightLoggerAnalysis", "pupil_util_path"); 
    Pi_util_path = getpref("lightLogger", "Pi_util_path"); 
    fprintf("\t world_util from path %s\n", world_util_path); 
    fprintf("\t pupil_util from path %s\n", pupil_util_path); 
    fprintf("\t Pi_util form path %s\n",  Pi_util_path); 

    world_util = import_pyfile(world_util_path);
    pupil_util = import_pyfile(pupil_util_path);
    Pi_util = import_pyfile(Pi_util_path); 

    fprintf("\t world_util actual path %s\n", char(py.str(py.getattr(world_util, '__file__'))));
    fprintf("\t pupil_util actual path %s\n", char(py.str(py.getattr(pupil_util, '__file__'))));
    fprintf("\t Pi_util actual path %s\n", char(py.str(py.getattr(Pi_util, '__file__'))));



    % Verify that the devices can be found
    disp("Calibration | Initializing devices...");

    % Establish a persistent BLE connection to the device so we do not
    % re-discover on every message
    disp("Calibration | Connecting to device over BLE...");
    bluetooth_client = py_call_module_attr(bluetooth_central, "connect_peripheral_matlab_wrapper", device_num);
    ble_cleanup = onCleanup(@() py_call_module_attr(bluetooth_central, "disconnect_peripheral_matlab_wrapper", bluetooth_client));
    disp("Calibration | BLE connection established.");

    % Set up and initialize the CombiLED
    % Save the combiExperiments path and add CombiLED Toolbox to the path
    % Return to the lightLogger project
    tbUseProject('lightLogger');
    addpath(current_dir);

    % Define the background we will use for the CombiLED 
    CombiLED_background = ones([1, 8]); 

    % Construct the calibration metadata based on the sensors and NDFs provided  
    % if a metadata struct is not already provided 
    if(~isstruct(calibration_metadata))
        calibration_metadata = initialize_calibration_data(calibration_metadata, ...
                                                           NDFs,...
                                                           CombiLED_background,...
                                                           world_util,...
                                                           pupil_util...
                                                          ); 
    else
        calibration_metadata = rebuild_calibration_metadata(calibration_metadata, world_util);
    end 

    disp("THIS IS CALIBRATION METADATA")
    disp(calibration_metadata)

    % Step 2: Initialize output directories 
    dropbox_savedir = fullfile(getpref('lightLogger', 'dropbox_calibration_dir'), sensor_ids, experiment_name); 
    cloud_output_dir = dropbox_savedir; 
    local_savedir = "";
    if(upload_video_data == false)
        cloud_output_dir = "";
    end 
        
    % Gather the target local savedir 
    default_external_output_dir = fullfile("/media/rpiControl/EXTERNAL_1/", experiment_name);
    external_output_dir = input(sprintf("Enter a location to save the local data [Enter for default=%s" + "]: ", default_external_output_dir), 's');
    if(strtrim(external_output_dir) == "")
        external_output_dir = default_external_output_dir; 
    end 

    % Step 3: Collect the calibration measurements

    % A. Calcualte the MS linearity at all NDF levels
    disp("Calibration | Collecting MS linearity with struct:");
    disp(calibration_metadata.ms_linearity)

    % Attempt to collect the calibration measurements 
    [success, ms_linearity_calibration_metadata] = collect_ms_linearity_data(device_num,...
                                                                             calibration_metadata.ms_linearity,...
                                                                             bluetooth_central,...
                                                                             bluetooth_client,...
                                                                             "",...
                                                                             cloud_output_dir,...
                                                                             external_output_dir,...
                                                                             @wait_for_combiLED_cooldown_between_NDFs...
                                                                            ); 
    calibration_metadata.ms_linearity = ms_linearity_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured, error cleanly
    if(~success)
        raise_calibration_failure("MS linearity calibration quit early on the light logger.", calibration_metadata.ms_linearity);
    end 
    

    tbUseProject('lightLogger');
    addpath(current_dir);

    % B. Calculate the world camera linearity at all NDF levels
    disp("Calibration | Collecting world camera linearity with struct:");
    disp(calibration_metadata.world_linearity)

    % Attempt to collect this data
    [success, world_linearity_calibration_metadata] = collect_world_linearity_data(device_num,...
                                                                                   calibration_metadata.world_linearity,...
                                                                                   bluetooth_central,...
                                                                                   bluetooth_client,...
                                                                                   "",...
                                                                                   cloud_output_dir,...
                                                                                   external_output_dir,...
                                                                                   @wait_for_combiLED_cooldown_between_NDFs...
                                                                                  );
    calibration_metadata.world_linearity = world_linearity_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    if(~success)
        raise_calibration_failure("World camera linearity calibration quit early on the light logger.", calibration_metadata.world_linearity);
    end

    tbUseProject('lightLogger');
    addpath(current_dir);

    % C. Calibrate the temporal sensitivity of the different sensors 
    %    Only collect this measurement at NDF [0, 5) 
    disp("Calibration | Collecting Temporal Sensitivity with struct: ");
    disp(calibration_metadata.temporal_sensitivity) ;
    
    % Attempt to collect this data
    [success, temporal_sensitivity_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                             calibration_metadata.temporal_sensitivity,...
                                                                                             bluetooth_central,...
                                                                                             bluetooth_client,...
                                                                                             "",...
                                                                                             cloud_output_dir,...
                                                                                             external_output_dir,...
                                                                                             @wait_for_combiLED_cooldown_between_NDFs...
                                                                                            ); 
    calibration_metadata.temporal_sensitivity = temporal_sensitivity_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured error cleanly                                                                                         
    if(~success)       
        raise_calibration_failure("Temporal sensitivity calibration quit early on the light logger.", calibration_metadata.temporal_sensitivity);
    end

    tbUseProject('lightLogger');
    addpath(current_dir);

    % D. Characterize the phase offset between the pairs of sensors 
    %    Only collect this measurement at the 1 NDF level 
    disp("Calibration | Collecting Phase Fitting with struct: ");
    disp(calibration_metadata.phase_fitting);

    % Attempt to collect this measurement 
    [success, phase_fitting_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                      calibration_metadata.phase_fitting,...
                                                                                      bluetooth_central,...
                                                                                      bluetooth_client,...
                                                                                      "PhaseFitting_",...
                                                                                      cloud_output_dir,...
                                                                                      external_output_dir,...
                                                                                      @wait_for_combiLED_cooldown_between_NDFs...
                                                                                     ); 
    calibration_metadata.phase_fitting = phase_fitting_calibration_metadata;    
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured error cleanly                                                                                   
    if(~success) 
        raise_calibration_failure("Phase fitting calibration quit early on the light logger.", calibration_metadata.phase_fitting);
    end 

    % E. Calculate the contrast gamma function
    % Only collect this measurement at the 1 NDF level
    disp("Calibration | Collecting Contrast Gamma with struct: ");
    disp(calibration_metadata.contrast_gamma);
    
    % Attempt to collect this measurement
    [success, contrast_gamma_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                       calibration_metadata.contrast_gamma,...
                                                                                       bluetooth_central,...
                                                                                       bluetooth_client,...
                                                                                       "ContrastGamma_",...
                                                                                       cloud_output_dir,...
                                                                                       external_output_dir,...
                                                                                       @wait_for_combiLED_cooldown_between_NDFs...
                                                                                      ); 
    calibration_metadata.contrast_gamma = contrast_gamma_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we've gotten in the recording 

    % If an error occured error cleanly                                                                                      
    if(~success)                                                                                    
        raise_calibration_failure("Contrast gamma calibration quit early on the light logger.", calibration_metadata.contrast_gamma);
    end

    % Step 4: Save object with all of the parameters to the folder on DropBox
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode)

end 

function wait_for_combiLED_cooldown_between_NDFs(CombiLED, NDF)
    cooldown_seconds = 10 * 60;
    zero_primaries = zeros(1, 8);

    CombiLED.setPrimaries(zero_primaries);
    fprintf("\nCalibration | Completed NDF %.3f. CombiLED primaries set to all zeros.\n", NDF);
    fprintf("Calibration | The CombiLED must cool down for 10 minutes before the next NDF.\n");

    for seconds_remaining = cooldown_seconds:-1:1
        minutes_remaining = floor(seconds_remaining / 60);
        seconds_in_minute = mod(seconds_remaining, 60);
        fprintf("\rCalibration | Cooldown remaining: %02d:%02d", minutes_remaining, seconds_in_minute);
        pause(1);
    end

    fprintf("\rCalibration | Cooldown remaining: 00:00\n");
end

% Local function to set up a Calibration Struct if it has not been passed in as an 
% argument 
function CalibrationData = initialize_calibration_data(CalibrationData,...
                                                       NDFs,...
                                                       background,...
                                                       world_util,...
                                                       pupil_util...
                                                      )

    % Initialize an object to save the parameters and metadata of this calibration data 
    CalibrationData = struct; 

    % Extract the world camera sensor mode 
    sensor_mode_idx = 0; 


    % A. Set up the substructs we will use for each of the Calibration measures 
    CalibrationData.ms_linearity = struct; 
    CalibrationData.world_linearity = struct;
    CalibrationData.temporal_sensitivity = struct; 
    CalibrationData.phase_fitting = struct; 
    CalibrationData.contrast_gamma = struct;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % { MS LINEARITY } - We will collect this between NDFs 0-6
    
    k_settings_levels = 10; % Define the number of settings levels to record
    n_measures = 3; % The number of measurements to make at a given settings level 
    settings_scalars = linspace(0.05, 0.95, k_settings_levels);  % Define the settings values we will explore 
    settings_levels_orders = randomize_settings_orders(numel(NDFs.ms_linearity), settings_scalars, n_measures); % The order we will set the settingsi n in, randomized per measure 
    recording_seconds = 12; % Define how long a given recording will be per setting 
    completed_measurements = false(numel(NDFs.ms_linearity), k_settings_levels, n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure
    
    % Build the sensors and settings per NDF
    sensors_and_settings = {}; 
    for ii = 1:numel(NDFs.ms_linearity)
        % Retrieve the given NDF 
        NDF = NDFs.ms_linearity(ii); 

        % Initialize the sensors we will use for this measurement 
        sensors.M.use_LED = false; 
        
        % Save this struct 
        sensors_and_settings{ii} = sensors; 
    end 
    clear NDF; 
    clear sensors; 

    % Save the metadata for this calibration measurement
    CalibrationData.ms_linearity.NDFs = NDFs.ms_linearity; 
    CalibrationData.ms_linearity.cal_files = cell(numel(NDFs.ms_linearity), 1);
    CalibrationData.ms_linearity.sensors_and_settings = sensors_and_settings; 
    CalibrationData.ms_linearity.background = background; 
    CalibrationData.ms_linearity.background_scalars = settings_scalars; 
    CalibrationData.ms_linearity.background_scalars_orders = settings_levels_orders; 
    CalibrationData.ms_linearity.n_measures = n_measures;
    CalibrationData.ms_linearity.recording_seconds = recording_seconds;  
    CalibrationData.ms_linearity.completed_measurements = completed_measurements; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % { WORLD CAMERA LINEARITY } - We will collect this between NDFs 0-4
    contrast_agc_targets = [0.1]; 
    k_settings_levels = 10; % Define the number of settings levels to record
    n_measures = 3; % The number of measurements to make at a given settings level
    settings_scalars = linspace(0.05, 0.95, k_settings_levels); % Define the settings values we will explore
    settings_levels_orders = randomize_settings_orders(numel(NDFs.world_linearity), settings_scalars, n_measures); % The order we will set the settings in, randomized per measure
    recording_seconds = 15; % Define how long a given recording will be per setting
    completed_measurements = false(numel(NDFs.world_linearity), numel(contrast_agc_targets), k_settings_levels, n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure

    % Build the sensors and settings per NDF
    clear sensors_and_settings;
    sensors_and_settings = struct();

    % For each of the AGC contrast targets, we will buold the settings 
    for tt = 1:numel(contrast_agc_targets)
        contrast_target = contrast_agc_targets(tt); 

        % Let's get the NDF settings for this contrast level    
        settings_contrast_dict = py_call_module_attr(world_util, "get_world_contrast_level_ndf_settings", contrast_target); 

        % Initialize the container for the settings per NDF for this contrast target 
        contrast_target_fieldname = "x" + strrep(string(contrast_target), ".", "x"); 
        agc_target_sensors_and_settings = {}; 

        for ii = 1:numel(NDFs.world_linearity)
            % Retrieve the given NDF
            NDF = NDFs.world_linearity(ii);
            
            % Build the world camera settings for this NDF level
            fixed_settings = double(py_getitem(settings_contrast_dict, py.float(NDF)));
            sensors.W.Again = fixed_settings(1);
            sensors.W.Dgain = fixed_settings(2);
            sensors.W.exposure = int32(fixed_settings(3));
            sensors.W.agc = false;
            sensors.W.save_agc_metadata = true;
            sensors.W.sensor_mode_idx = sensor_mode_idx;
            sensors.W.awb = false;
            sensors.W.noise_mode = false;

            % Save this struct
            agc_target_sensors_and_settings{ii} = sensors;
        end
        sensors_and_settings.(contrast_target_fieldname) = agc_target_sensors_and_settings; 
    
        clear NDF;
        clear sensors;
    end 

    % Save the metadata for this calibration measurement
    CalibrationData.world_linearity.NDFs = NDFs.world_linearity;
    CalibrationData.world_linearity.contrast_agc_targets = contrast_agc_targets; 
    CalibrationData.world_linearity.cal_files = cell(numel(NDFs.world_linearity), 1);
    
    CalibrationData.world_linearity.sensors_and_settings = sensors_and_settings;
    
    
    CalibrationData.world_linearity.background = background;
    CalibrationData.world_linearity.background_scalars = settings_scalars;
    CalibrationData.world_linearity.background_scalars_orders = settings_levels_orders;
    CalibrationData.world_linearity.n_measures = n_measures;
    CalibrationData.world_linearity.recording_seconds = recording_seconds;
    CalibrationData.world_linearity.completed_measurements = completed_measurements;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % { TEMPORAL SENSITIVITY } - We will collect this for NDFs [0, 5)
    recording_seconds = 12; % Define the length in seconds of the recordings
    contrast_levels = [0.5]; % Define the list of contrast levels we will explore
    frequencies = [0.25, 0.5, 1, 3, 6, 12, 25, 50, 100 ]; % Define the list of frequencies we will modulate at (1/29 is a test for the standard 30 second chunks size)
    n_measures = 3;  % The number of measurements to make at each combination of contrast + frequency 
    contrast_levels_orders = randomize_contrast_orders(numel(NDFs.temporal_sensitivity), numel(contrast_levels), n_measures); % The order we will expose contrasts in, randomized per measure 
    frequencies_orders = randomize_frequency_orders(numel(NDFs.temporal_sensitivity), n_measures, numel(contrast_levels), numel(frequencies)); % The orders in which we will expose the frequencies at a given contrast, randomized per contrast
    completed_measurements = false(numel(NDFs.temporal_sensitivity), numel(contrast_levels), numel(frequencies), n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure

    % Initialize the sensors we will use for this measurement 
    clear sensors_and_settings; 

    % Build the sensors and settings per NDF
    sensors_and_settings = {}; 
    for ii = 1:numel(NDFs.temporal_sensitivity)
        % Retrieve the given NDF 
        NDF = NDFs.temporal_sensitivity(ii); 

        disp("NDF")

        % Create some mapping for the gain and exposure of the world 
        % per NDF level so it behaves properly at the start of the recording 
        % (no warmup needed)
        world_ndf_level_settings = py_module_attr(world_util, "WORLD_NDF_LEVEL_SETTINGS");
        fixed_settings = double(py_getitem(world_ndf_level_settings, py.float(NDF)));
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2);
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  % TODO: Turn this back on 
        sensors.W.save_agc_metadata = true; 
        sensors.W.sensor_mode_idx = sensor_mode_idx;
        sensors.W.awb = false; 
        sensors.W.noise_mode = false; 

        disp(sensors.W)
        
        % Save this struct 
        sensors_and_settings{ii} = sensors; 
    end 

    clear NDF; 
    clear sensors; 

    % Save the metadata for this calibration measurement
    CalibrationData.temporal_sensitivity.NDFs = NDFs.temporal_sensitivity; 
    CalibrationData.temporal_sensitivity.cal_files = cell(numel(NDFs.temporal_sensitivity), 1);
    CalibrationData.temporal_sensitivity.background = background; 
    CalibrationData.temporal_sensitivity.sensors_and_settings = sensors_and_settings; 
    CalibrationData.temporal_sensitivity.modulations = cell(numel(NDFs.temporal_sensitivity), 1); 
    CalibrationData.temporal_sensitivity.recording_seconds = recording_seconds;
    CalibrationData.temporal_sensitivity.completed_measurements = completed_measurements; 
    CalibrationData.temporal_sensitivity.contrast_levels = contrast_levels; 
    CalibrationData.temporal_sensitivity.contrast_levels_orders = contrast_levels_orders; 
    CalibrationData.temporal_sensitivity.frequencies = frequencies; 
    CalibrationData.temporal_sensitivity.frequencies_orders = frequencies_orders; 
    CalibrationData.temporal_sensitivity.n_measures = n_measures; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % { PHASE FITTING } - We will collect this for only NDF 1 
    contrast_levels = [0.5]; % Define the list of contrast levels we will explore
    frequencies = [0.1]; % Define the list of frequencies we will modulate at 
    n_measures = 3; % The number of measurements to make at each combination of contrast + frequency 
    contrast_levels_orders = randomize_contrast_orders(numel(NDFs.phase_fitting), numel(contrast_levels), n_measures); % The order we will expose contrasts in, randomized per measure 
    frequencies_orders = randomize_frequency_orders(numel(NDFs.phase_fitting), n_measures, numel(contrast_levels), numel(frequencies)); % The orders in which we will expose the frequencies at a given contrast, randomized per contrast
    recording_seconds = 12; % The duration to record for after warming up 
    completed_measurements = false(numel(NDFs.phase_fitting), numel(contrast_levels), numel(frequencies), n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure

    % Initialize the sensors for this measurement. Crucially, turn the AGC OFF
    clear sensors_and_settings; 
    
    % Build the sensors and settings per NDF
    sensors_and_settings = {}; 
    for ii = 1:numel(NDFs.phase_fitting)
        % Retrieve the given NDF 
        NDF = NDFs.phase_fitting(ii); 

        % Build the settings for the M at this NDF level
        sensors.M.use_LED = false; % Turn off indicator light for calibration (so the light doesn't illuminate the sphere)

        % Build the settings for the W at this NDF level 
        world_ndf_level_settings = py_module_attr(world_util, "WORLD_NDF_LEVEL_SETTINGS");
        fixed_settings = double(py_getitem(world_ndf_level_settings, py.float(NDF)));
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2); 
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  
        sensors.W.save_agc_metadata = true;
        sensors.W.sensor_mode_idx = sensor_mode_idx;
        sensors.W.awb = false; 
        sensors.W.noise_mode = false; 
        
        % Save this struct 
        sensors_and_settings{ii} = sensors; 
    end 
    clear NDF; 
    clear sensors; 

    % Save the metadata for this calibration measurement
    CalibrationData.phase_fitting.NDFs = NDFs.phase_fitting;  
    CalibrationData.phase_fitting.cal_files = cell(numel(NDFs.phase_fitting), 1); 
    CalibrationData.phase_fitting.background = background; 
    CalibrationData.phase_fitting.sensors_and_settings = sensors_and_settings; 
    CalibrationData.phase_fitting.modulations = cell(numel(NDFs.phase_fitting), 1); 
    CalibrationData.phase_fitting.recording_seconds = recording_seconds;
    CalibrationData.phase_fitting.completed_measurements = completed_measurements;
    CalibrationData.phase_fitting.contrast_levels = contrast_levels;  
    CalibrationData.phase_fitting.contrast_levels_orders = contrast_levels_orders;  
    CalibrationData.phase_fitting.frequencies = frequencies; 
    CalibrationData.phase_fitting.frequencies_orders = frequencies_orders; 
    CalibrationData.phase_fitting.n_measures = n_measures; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % { CONTRAST GAMMA } - We will collect this only at NDF1 
    k_contrast_levels = 10; % Define the number of contrast levels to measure
    contrast_levels = linspace(0.05, 0.95, k_contrast_levels); % Define the contrast levels we will measure at 
    frequencies = [10]; % Define the frequencies per contrast level 
    n_measures = 3; % Define the number of measurements to make at each combination of contrast + frequency 
    contrast_levels_orders = randomize_contrast_orders(numel(NDFs.contrast_gamma), numel(contrast_levels), n_measures); % The order we will expose contrasts in, randomized per measure 
    frequencies_orders = randomize_frequency_orders(numel(NDFs.contrast_gamma), n_measures, numel(contrast_levels), numel(frequencies)); % The orders in which we will expose the frequencies at a given contrast, randomized per contrast
    completed_measurements = false(numel(NDFs.contrast_gamma), numel(contrast_levels), numel(frequencies), n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure

    % Initailize the sensors for this measurement. Use only world and turn agc OFF 
    clear sensors_and_settings; 

    % Build the sensors and settings per NDF
    sensors_and_settings = {}; 
    for ii = 1:numel(NDFs.contrast_gamma)
        % Retrieve the given NDF 
        NDF = NDFs.contrast_gamma(ii); 

        % Build the settings for this NDF level 
        world_ndf_level_settings = py_module_attr(world_util, "WORLD_NDF_LEVEL_SETTINGS");
        fixed_settings = double(py_getitem(world_ndf_level_settings, py.float(NDF)));
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2);
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  
        sensors.W.save_agc_metadata = true; 
        sensors.W.sensor_mode_idx = sensor_mode_idx;
        sensors.W.awb = false; 
        sensors.W.noise_mode = false; 
        
        % Save this struct 
        sensors_and_settings{ii} = sensors; 
    end 
    clear NDF; 
    clear sensors; 

    % Save the metadata for this calibration measurement 
    CalibrationData.contrast_gamma.NDFs = NDFs.contrast_gamma; 
    CalibrationData.contrast_gamma.cal_files = cell(numel(NDFs.contrast_gamma), 1); 
    CalibrationData.contrast_gamma.background = background; 
    CalibrationData.contrast_gamma.sensors_and_settings = sensors_and_settings; 
    CalibrationData.contrast_gamma.modulations = cell(numel(NDFs.contrast_gamma), 1); ; 
    CalibrationData.contrast_gamma.recording_seconds = recording_seconds;
    CalibrationData.contrast_gamma.completed_measurements = completed_measurements; 
    CalibrationData.contrast_gamma.contrast_levels = contrast_levels;  
    CalibrationData.contrast_gamma.contrast_levels_orders = contrast_levels_orders;  
    CalibrationData.contrast_gamma.frequencies = frequencies; 
    CalibrationData.contrast_gamma.frequencies_orders = frequencies_orders; 
    CalibrationData.contrast_gamma.n_measures = n_measures; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    disp("THIS IS CALIBRATION DATA BEFORE I RETURN")
    disp(CalibrationData)

    return ; 

end 

function calibration_metadata = rebuild_calibration_metadata(loaded_calibration_metadata, world_util)
    % Rebuild a fresh calibration metadata struct while preserving prior
    % progress and outputs. The main thing we need to refresh is the
    % world-camera sensor settings payload so resumed runs don't carry stale
    % sensor_mode structs or embedded Python-derived state.
    calibration_metadata = loaded_calibration_metadata;

    calibration_fields = {"world_linearity", "temporal_sensitivity", "phase_fitting", "contrast_gamma"};
    for field_idx = 1:numel(calibration_fields)
        calibration_field = calibration_fields{field_idx};
        if(~isfield(calibration_metadata, calibration_field))
            continue;
        end

        if(~isfield(calibration_metadata.(calibration_field), "sensors_and_settings"))
            continue;
        end

        calibration_metadata.(calibration_field).sensors_and_settings = ...
            rebuild_sensor_settings_collection( ...
                loaded_calibration_metadata.(calibration_field).sensors_and_settings, ...
                world_util ...
            );
    end
end

function rebuilt_collection = rebuild_sensor_settings_collection(loaded_collection, world_util)
    if(iscell(loaded_collection))
        rebuilt_collection = cell(size(loaded_collection));
        for idx = 1:numel(loaded_collection)
            rebuilt_collection{idx} = rebuild_single_sensor_settings_entry(loaded_collection{idx}, world_util);
        end
        return;
    end

    if(isstruct(loaded_collection))
        rebuilt_collection = struct;
        field_names = fieldnames(loaded_collection);
        for idx = 1:numel(field_names)
            field_name = field_names{idx};
            rebuilt_collection.(field_name) = rebuild_sensor_settings_collection(loaded_collection.(field_name), world_util);
        end
        return;
    end

    rebuilt_collection = loaded_collection;
end

function sensors = rebuild_single_sensor_settings_entry(loaded_sensors, world_util)
    sensors = loaded_sensors;
    if(~isstruct(sensors))
        return;
    end

    if(isfield(sensors, "W"))
        if(~isstruct(sensors.W))
            sensors.W = struct;
        end
        sensors.W = rebuild_world_sensor_settings(sensors.W, world_util);
    end
end

function world_sensor_settings = rebuild_world_sensor_settings(world_sensor_settings, world_util)
    default_sensor_mode_idx = 0;
    resolved_sensor_mode_idx = default_sensor_mode_idx;

    if(isfield(world_sensor_settings, "sensor_mode_idx"))
        candidate_idx = world_sensor_settings.sensor_mode_idx;
        if(isnumeric(candidate_idx) && isscalar(candidate_idx) && ~isnan(candidate_idx))
            resolved_sensor_mode_idx = double(candidate_idx);
        end
    elseif(isfield(world_sensor_settings, "sensor_mode"))
        resolved_sensor_mode_idx = resolve_world_sensor_mode_idx(world_sensor_settings.sensor_mode, world_util, default_sensor_mode_idx);
    end

    if(isfield(world_sensor_settings, "sensor_mode"))
        world_sensor_settings = rmfield(world_sensor_settings, "sensor_mode");
    end
    world_sensor_settings.sensor_mode_idx = resolved_sensor_mode_idx;
end

function sensor_mode_idx = resolve_world_sensor_mode_idx(sensor_mode, world_util, fallback_idx)
    sensor_mode_idx = fallback_idx;

    if(~isstruct(sensor_mode))
        return;
    end

    world_camera_custom_modes = py_module_attr(world_util, "WORLD_CAMERA_CUSTOM_MODES");
    n_modes = int64(py.len(world_camera_custom_modes));
    for mode_idx = 0:(n_modes - 1)
        candidate_mode = py_getitem(world_camera_custom_modes, int32(mode_idx));
        if(world_sensor_modes_match(sensor_mode, candidate_mode))
            sensor_mode_idx = double(mode_idx);
            return;
        end
    end
end

function tf = world_sensor_modes_match(matlab_sensor_mode, py_sensor_mode)
    tf = strings_match(fetch_struct_string(matlab_sensor_mode, "format"), fetch_py_string(py_sensor_mode, "format")) && ...
         strings_match(fetch_struct_string(matlab_sensor_mode, "unpacked"), fetch_py_string(py_sensor_mode, "unpacked")) && ...
         scalars_match(fetch_struct_scalar(matlab_sensor_mode, "bit_depth"), fetch_py_scalar(py_sensor_mode, "bit_depth")) && ...
         scalars_match(fetch_struct_scalar(matlab_sensor_mode, "fps"), fetch_py_scalar(py_sensor_mode, "fps")) && ...
         vectors_match(fetch_struct_vector(matlab_sensor_mode, "size"), fetch_py_vector(py_sensor_mode, "size")) && ...
         vectors_match(fetch_struct_vector(matlab_sensor_mode, "crop_limits"), fetch_py_vector(py_sensor_mode, "crop_limits"));
end

function value = fetch_struct_string(input_struct, field_name)
    value = "";
    if(isfield(input_struct, field_name))
        value = string(input_struct.(field_name));
    end
end

function value = fetch_py_string(py_mapping, key_name)
    value = string(py.str(py_getitem(py_mapping, key_name)));
end

function value = fetch_struct_scalar(input_struct, field_name)
    value = NaN;
    if(isfield(input_struct, field_name))
        field_value = input_struct.(field_name);
        if(isnumeric(field_value) && isscalar(field_value))
            value = double(field_value);
        end
    end
end

function value = fetch_py_scalar(py_mapping, key_name)
    value = double(py_getitem(py_mapping, key_name));
end

function value = fetch_struct_vector(input_struct, field_name)
    value = [];
    if(isfield(input_struct, field_name))
        field_value = input_struct.(field_name);
        if(isnumeric(field_value))
            value = double(field_value(:).');
        end
    end
end

function value = fetch_py_vector(py_mapping, key_name)
    py_vector = py_getitem(py_mapping, key_name);
    n_vals = int64(py.len(py_vector));
    value = zeros(1, n_vals);
    for idx = 1:n_vals
        value(idx) = double(py_getitem(py_vector, int32(idx - 1)));
    end
end

function tf = strings_match(lhs, rhs)
    tf = strlength(lhs) > 0 && strlength(rhs) > 0 && lhs == rhs;
end

function tf = scalars_match(lhs, rhs)
    tf = ~isnan(lhs) && ~isnan(rhs) && abs(lhs - rhs) < 1e-9;
end

function tf = vectors_match(lhs, rhs)
    tf = ~isempty(lhs) && isequal(size(lhs), size(rhs)) && all(abs(lhs - rhs) < 1e-9);
end

% Local function to generate a matrix of randomized orders of settings scalars 
function order_mat = randomize_settings_orders(n_NDFs, settings_scalars, n_measures)
    % First, let's generate a mat that is the number of measurements (the rows)
    % by number of settings scalars (the cols)
    order_mat = zeros([n_NDFs, n_measures, numel(settings_scalars)]); 

   % Next, we will iterate over the NDFs/ measurements and randomize 
    % the order per measurement 
    for nn = 1:n_NDFs
        for mm = 1:n_measures
            % Calculate a random permutation of the indices 
            random_order = randperm(numel(settings_scalars));

            % Save this order in to the order mat 
            order_mat(nn, mm, :) = random_order;
        end 
    end 

    return ; 

end 

% Local function to generate a matrix of randomized 
% contrast orders
function order_mat = randomize_contrast_orders(n_NDFs, n_contrast_levels, n_measures)
    % First, let's generate a mat that is the number of measurements (the rows)
    % by the number of contrast levels (the cols)
    order_mat = zeros([n_NDFs, n_measures, n_contrast_levels]); 

    % Next, we will iterate over the measurements and randomize 
    % the order per measurement 
    for nn = 1:n_NDFs
        for mm = 1:n_measures
            % Calculate a random permutation of the indices 
            random_order = randperm(n_contrast_levels);

            % Save this order in to the order mat 
            order_mat(nn, mm, :) = random_order;
        end 
    end 

    return ; 
end 

% Local function to generate a matrix of randomized frequency orders 
% For each of the contrast levels 
function order_mat = randomize_frequency_orders(n_NDFs, n_measures, n_contrast_levels, n_frequencies)
    % Let's first initialize an array to hold the frequency order
    % for each measure and each contrast level at each measurement
    order_mat = zeros([n_NDFs, n_measures, n_contrast_levels, n_frequencies]); 

    % Then, let's generate the orders 
    for nn = 1:n_NDFs
        for mm = 1:n_measures
            for cc = 1:n_contrast_levels
                % Generate a random ordering of the frequencies for 
                % this measurement and contrast level 
                order_mat(nn, mm, cc, :) = randperm(n_frequencies);
            end 
        end 
    end 

    return ; 
end 

function value = py_module_attr(module, attr_name)
    value = py.getattr(module, attr_name);
end

function value = py_call_module_attr(module, attr_name, varargin)
    callable = py_module_attr(module, attr_name);
    value = callable(varargin{:});
end

function value = py_getitem(py_obj, key)
    getitem = py.getattr(py_obj, "__getitem__");
    value = getitem(key);
end

% Local function to upload the CalibrationData struct. Used on error 
% and at the end of the script 
function upload_calibration_data(CalibrationData, dropbox_savedir, upload_mode)
    caldata_dropbox_path = fullfile(dropbox_savedir, "calibration_metadata.mat.bytes"); 
    py_call_module_attr(upload_mode, "upload_file", py.bytes(getByteStreamFromArray(CalibrationData)), caldata_dropbox_path);

end 

function raise_calibration_failure(prefix, calibration_substruct)
    detailed_message = "";
    if(isstruct(calibration_substruct) && isfield(calibration_substruct, "last_error_message"))
        detailed_message = string(calibration_substruct.last_error_message);
    end

    if(strlength(strtrim(detailed_message)) > 0)
        error("ERROR: %s\n%s", prefix, detailed_message);
    end

    error("ERROR: %s", prefix);
end
