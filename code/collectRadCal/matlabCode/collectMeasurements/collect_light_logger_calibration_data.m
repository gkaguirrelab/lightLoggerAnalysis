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
%   linearity of the light sensing chips of the MS, and 
%   2. measure the temporal sensitivity of all sensors.
%   3. Measure the offset between sensors in phase/time of all sensors 
%   4. Measure the contrast gamma of the sensors of the world camera. 
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
%   NDF                   - Int. The NDF that is present for 
%                           this calibration 
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
   
    % Ensure the sensor ids is a vector with at most 3 elements 
    % [W_id, P_id, M_id]
    sensor_ids = sensor_ids(:); % Flatten the id vector 
    assert(numel(sensor_ids) > 0 && numel(sensor_ids) <= 3); % Ensure it has the 3 elements

    % Build the experiment name from the sensor ids by interveaving 
    % the sensor id number after the sensor initial
    sensor_ids = reshape( [ 'WPM' ; num2str(sensor_ids)' ], 1 , []); 

    disp("Calibration | Importing libraries...");
    tbUseProject('lightLogger'); 

    % Import the Python library used to communicate 
    % via Bluetooth with the peripheral 
    bluetooth_central = import_pyfile(getpref("lightLogger", "bluetooth_central_path"));

    % Import the Python library used for uploading to Dropbox 
    upload_mode = import_pyfile(getpref("lightLogger", "upload_mode_path"));

    % Import the utility files of the world and pupil cameras 
    world_util = import_pyfile(getpref("lightLogger", "world_util_path"));
    pupil_util = import_pyfile(getpref("lightLogger", "pupil_util_path"));
    Pi_util = import_pyfile(getpref("lightLogger", "Pi_util_path")); 

    % Verify that the devices can be found 
    disp("Calibration | Initializing devices...");
    
    % Set up and initialize the CombiLED 
    % Save the combiExperiments path and add CombiLED Toolbox to the path 
    combiExperiments_path = "~/Documents/MATLAB/combiExperiments/";

    % Return to the lightLogger project
    tbUseProject('lightLogger'); 

    % Then, make sure the bluetooth device is available to connected 
    if(~bluetooth_central.peripheral_is_available_matlab_wrapper(device_num)) 
        error("ERROR: Light logger bluetooth probe failed");
    end 

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
    end 

    % Step 2: Initialize output directories 
    dropbox_savedir = fullfile(getpref('lightLogger', 'dropbox_calibration_dir'), sensor_ids, experiment_name); 
    cloud_output_dir = dropbox_savedir; 
    local_savedir = "";
    if(upload_video_data == false)
        cloud_output_dir = "";
    end 
        
    % Gather the target local savedir 
    default_external_output_dir = fullfile("/media/rpiControl/T7 Shield/", experiment_name);
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
                                                                             "",...
                                                                             cloud_output_dir,...
                                                                             external_output_dir...
                                                                            ); 
    calibration_metadata.ms_linearity = ms_linearity_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured, error cleanly
    if(~success)
        error("ERROR: MS linearity calibration quit early due to an unknown error on the light logger");
    end 
    

    tbUseProject('lightLogger');

    % B. Calibrate the temporal sensitivity of the different sensors 
    %    Only collect this measurement at NDF [0, 5) 
    disp("Calibration | Collecting Temporal Sensitivity with struct: ");
    disp(calibration_metadata.temporal_sensitivity) ;
    
    % Attempt to collect this data
    [success, temporal_sensitivity_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                             calibration_metadata.temporal_sensitivity,...
                                                                                             bluetooth_central,...
                                                                                             "",...
                                                                                             cloud_output_dir,...
                                                                                             external_output_dir...
                                                                                            ); 
    calibration_metadata.temporal_sensitivity = temporal_sensitivity_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured error cleanly                                                                                         
    if(~success)       
        error("ERROR: temporal sensitivity calibration quit early due to an unknown error on the light logger");
    end

    tbUseProject('lightLogger');

    % C. Characterize the phase offset between the pairs of sensors 
    %    Only collect this measurement at the 1 NDF level 
    disp("Calibration | Collecting Phase Fitting with struct: ");
    disp(calibration_metadata.phase_fitting);

    % Attempt to collect this measurement 
    [succces, phase_fitting_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                      calibration_metadata.phase_fitting,...
                                                                                      bluetooth_central,...
                                                                                      "PhaseFitting_",...
                                                                                      cloud_output_dir,...
                                                                                      external_output_dir...
                                                                                     ); 
    calibration_metadata.phase_fitting = phase_fitting_calibration_metadata;    
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we got on this recording

    % If an error occured error cleanly                                                                                   
    if(~succces) 
        error("ERROR: phase fitting calibration quit early due to an unknown error on the light logger");
    end 

    % D. Calculate the contrast gamma function
    % Only collect this measurement at the 1 NDF level
    disp("Calibration | Collecting Contrast Gamma with struct: ");
    disp(calibration_metadata.contrast_gamma);
    
    % Attempt to collect this measurement
    [success, contrast_gamma_calibration_metadata] = collect_temporal_sensitivity_data(device_num,...
                                                                                       calibration_metadata.contrast_gamma,...
                                                                                       bluetooth_central,...
                                                                                       "ContrastGamma_",...
                                                                                       cloud_output_dir,...
                                                                                       external_output_dir...
                                                                                      ); 
    calibration_metadata.contrast_gamma = contrast_gamma_calibration_metadata;
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode); % Upload how far we've gotten in the recording 

    % If an error occured error cleanly                                                                                      
    if(~success)                                                                                    
        error("ERROR: phase fitting calibration quit early due to an unknown error on the light logger");
    end

    % Step 4: Save object with all of the parameters to the folder on DropBox
    upload_calibration_data(calibration_metadata, dropbox_savedir, upload_mode)

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
    sensor_mode = world_util.WORLD_CAMERA_CUSTOM_MODES{1}; 

    % A. Set up the substructs we will use for each of the Calibration measures 
    CalibrationData.ms_linearity = struct; 
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
        fixed_settings = double(world_util.WORLD_NDF_LEVEL_SETTINGS{NDF});
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2);
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  % TODO: Turn this back on 
        sensors.W.save_agc_metadata = true; 
        sensors.W.sensor_mode = sensor_mode;
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
        fixed_settings = double(world_util.WORLD_NDF_LEVEL_SETTINGS{NDF});
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2); 
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  
        sensors.W.save_agc_metadata = true;
        sensors.W.sensor_mode = sensor_mode;
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
    contrast_levels_orders = randomize_contrast_orders(numel(NDFs.contrast_gamma), numel(contrast_levels), n_measures); % The order we will expose contrasts in, randomized per measure 
    frequencies_orders = randomize_frequency_orders(numel(NDFs.contrast_gamma), n_measures, numel(contrast_levels), numel(frequencies)); % The orders in which we will expose the frequencies at a given contrast, randomized per contrast
    n_measures = 3; % Define the number of measurements to make at each combination of contrast + frequency 
    completed_measurements = false(numel(NDFs.contrast_gamma), numel(contrast_levels), numel(frequencies), n_measures); % Define a boolean matrix of completed measurements. We will use this to resume recording on failure

    % Initailize the sensors for this measurement. Use only world and turn agc OFF 
    clear sensors_and_settings; 

    % Build the sensors and settings per NDF
    sensors_and_settings = {}; 
    for ii = 1:numel(NDFs.contrast_gamma)
        % Retrieve the given NDF 
        NDF = NDFs.contrast_gamma(ii); 

        % Build the settings for this NDF level 
        fixed_settings = double(world_util.WORLD_NDF_LEVEL_SETTINGS{NDF});
        sensors.W.Again = fixed_settings(1);
        sensors.W.Dgain = fixed_settings(2);
        sensors.W.exposure = int32(fixed_settings(3)); 
        sensors.W.agc = false;  
        sensors.W.save_agc_metadata = true; 
        sensors.W.sensor_mode = sensor_mode;
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


    return ; 

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

% Local function to upload the CalibrationData struct. Used on error 
% and at the end of the script 
function upload_calibration_data(CalibrationData, dropbox_savedir, upload_mode)
    caldata_dropbox_path = fullfile(dropbox_savedir, "calibration_metadata.mat.bytes"); 
    upload_mode.upload_file(py.bytes(getByteStreamFromArray(CalibrationData)), caldata_dropbox_path);

end 

