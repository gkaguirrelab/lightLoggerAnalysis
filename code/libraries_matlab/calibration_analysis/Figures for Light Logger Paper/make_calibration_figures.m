function make_calibration_figures

%% 1. Plot MS & TS linearity
% load the data
dropBoxBaseDir = getpref('combiExperiments','dropboxBaseDir');
light_logger_calibration_data = load([dropBoxBaseDir, '/FLIC_data/LightLogger_RadCal/W1P1M1/MSOnly_TS_range_slightly_weird/calibration_data.mat']);
% Struct that has both metadata and parsed readings from the experiment

% First, extract the broad subfields of the calibration information and the parsed readings
% from the input struct
calibration_metadata = light_logger_calibration_data.MSOnly_imperfect_TS_data.metadata; % will need to be updated if using a different file
parsed_readings = light_logger_calibration_data.MSOnly_imperfect_TS_data.readings;

% Analyze the MS linearity readings if there are any to analyze.
if(numel(calibration_metadata.ms_linearity.NDFs > 0))
    analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
        parsed_readings.ms_linearity,...
        'plotSettingLevel', false...
        );
end

end