function analyze_light_logger_calibration_data(light_logger_calibration_data)
% Analyze the results of a light logger calibration measurement (post-conversion)
%
% Syntax:
%  analyze_light_logger_calibration_data(light_logger_calibration_data)
%
% Description:
%   Given a converted light logger calibration data struct, 
%   analyze the various components of this measurement. That is, 
%   the ms linearity component, the temporal sensitivity component, 
%   the phase alignment component, and the contrast gamma component. 
%
%
% Inputs:
%   light_logger_calibration_data   - Struct. Converted metadata + parsed_readings
%                                     for the calibration measurement. 
%                              
%
% Examples:
%{
    path_to_experiment = "/example/path"; 
    converted_light_logger_data = convert_light_logger_calibration_data(path_to_experiment, true, true true, true); 
    analyze_light_logger_calibration(converted_light_logger_data);
%}
    arguments 
        light_logger_calibration_data; % Struct that has both converted metadata and readings from the experiment 
    end 

    % First, extract the broad subfields of the calibration information and the parsed readings
    % from the input struct
    calibration_metadata = light_logger_calibration_data.metadata;
    parsed_readings = light_logger_calibration_data.readings;

    % 1. Analyze the MS linearity readings if there are any to analyze.
    if(numel(calibration_metadata.ms_linearity.NDFs) > 0) 
        analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
                                  parsed_readings.ms_linearity...
                                 );
    end                             


    % 2. Analyze the Temporal Sensitivity readings
    if(numel(calibration_metadata.temporal_sensitivity.NDFs) > 0)
        analyze_temporal_sensitivity_data(calibration_metadata.temporal_sensitivity, ...
                                          parsed_readings.temporal_sensitivity...
                                         );
    end

    % 3. Analyze the Phase fitting readings
    if(numel(calibration_metadata.phase_fitting.NDFs) > 0)
        analyze_phase_fit_data(calibration_metadata.phase_fitting,...
                               parsed_readings.phase_fitting...
                              );
    end

    % 4. Analyze the contrast gamma readings
    if(numel(calibration_metadata.contrast_gamma.NDFs) > 0)
        analyze_contrast_gamma_data(calibration_metadata.contrast_gamma,...
                                    parsed_readings.contrast_gamma...
                                   );
    end
    

end

