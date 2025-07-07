function analyze_light_logger_calibration_data(light_logger_calibration_data)
    arguments 
        light_logger_calibration_data; % Struct that has both metadata and parsed readings from the experiment 
    end 



    % First, extract the broad subfields of the calibration information and the parsed readings
    % from the input struct
    calibration_metadata = light_logger_calibration_data.metadata;
    parsed_readings = light_logger_calibration_data.readings;

    % 1. Analyze the MS linearity readings if there are any to analyze.
    if(numel(calibration_metadata.ms_linearity.NDFs > 0)) 
        analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
                                  parsed_readings.ms_linearity...
                                 );
    end                             


    % 2. Analyze the Temporal Sensitivity readings
    if ~isempty(calibration_metadata.temporal_sensitivity.frequencies)
        analyze_temporal_sensitivity_data(calibration_metadata.temporal_sensitivity, ...
                                          parsed_readings.temporal_sensitivity...
                                         );
    end

    %{
    % 3. Analyze the Phase fitting readings
    if ~isempty(calibration_metadata.phase_fitting.frequencies)
        temporal_offsets_secs = analyze_phase_fit_data(calibration_metadata.phase_fitting, parsed_readings.phase_fitting);

        % Report the temporal offset values and standard deviations
        fprintf('Mean and std of world - pupil temporal offset [ms]: %2.5f, %2.5f\n',1000*mean(temporal_offsets_secs('W-P')),1000*std(temporal_offsets_secs('W-P')))
        fprintf('Mean and std of world - AS temporal offset [ms]: %2.5f, %2.5f\n',1000*mean(temporal_offsets_secs('W-AS')),1000*std(temporal_offsets_secs('W-AS')))
    end

    % 4. Analyze the contrast gamma readings
    if ~isempty(calibration_metadata.contrast_gamma.frequencies)
        analyze_contrast_gamma_data(calibration_metadata.contrast_gamma, parsed_readings.contrast_gamma, calibration_metadata.ndf);
    end
    %}

end

