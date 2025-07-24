function make_calibration_figures

% get path to directory with calibration data
drop_box_base_dir = getpref('combiExperiments','dropboxBaseDir');
cal_data_path = '/FLIC_data/LightLogger_RadCal/W1P1M1/';

%% 1. Plot AS & TS linearity
% load the AS data
light_logger_calibration_data = load([drop_box_base_dir, cal_data_path, 'MSOnly_TS_range_slightly_weird/calibration_data.mat']);
% Struct that has both metadata and parsed readings from the experiment

% First, extract the broad subfields of the calibration information and the parsed readings
% from the input struct
MS_metadata = light_logger_calibration_data.MSOnly_imperfect_TS_data.metadata; % will need to be updated if using a different file
MS_parsed_readings = light_logger_calibration_data.MSOnly_imperfect_TS_data.readings;

% Analyze the MS linearity readings if there are any to analyze.
if(numel(MS_metadata.ms_linearity.NDFs > 0))
    analyze_ms_linearity_data(MS_metadata.ms_linearity,...
        MS_parsed_readings.ms_linearity,...
        'plotSettingLevel', false...
        );
end

%% 2. Plot gamma data
% load gamma data
gamma_calibration_data = load([drop_box_base_dir, cal_data_path, 'contrast_gamma_data/contrast_gamma_data.mat']);
gamma_metadata = gamma_calibration_data.contrast_gamma_data.metadata.contrast_gamma; % will need to be updated if using a different file
gamma_parsed_readings = gamma_calibration_data.contrast_gamma_data.readings.contrast_gamma;

% plot gamma data
analyze_contrast_gamma_data(gamma_metadata, gamma_parsed_readings);

%% 3. Plot phase data
% load gamma data
% currently does not make plots :(
phase_calibration_data = load([drop_box_base_dir, cal_data_path, 'phase_alignment_data/phase_alignment_data.mat']);
phase_metadata = phase_calibration_data.phase_alignment_calibration_data.metadata.phase_fitting; % will need to be updated if using a different file
phase_parsed_readings = phase_calibration_data.phase_alignment_calibration_data.readings.phase_fitting;
analyze_phase_fit_data(phase_metadata, phase_parsed_readings)

%% 4. Plot temporal sensitivity
world_ttf_data = load([drop_box_base_dir, cal_data_path, 'world_ttf_data_noAGC/world_ttf_data.mat']);
world_ttf_metadata = world_ttf_data.world_ttf_data.metadata.temporal_sensitivity; % will need to be updated if using a different file
world_ttf_readings = world_ttf_data.world_ttf_data.readings.temporal_sensitivity;

if(numel(world_ttf_metadata.NDFs) > 0)
    analyze_temporal_sensitivity_data(world_ttf_metadata, ...
        world_ttf_readings,...
        'plotAGC', false , 'plotEachMeasure', false ...
        );
end

end