function [msLinearityFigHandles, ttfFigHandle] = analyze_light_logger_calibration_data(LightLoggerCalibrationData,...
                                                                                       whichMSChannelsToPlot,...
                                                                                       msLinearityFigHandles,...
                                                                                       ttfFigHandle...
                                                                                      )
%{
    whichMSChannelsToPlot = {[],[]};
    msLinearityFigHandles = {};
    ttfFigHandle = [];
    load('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_data/LightLogger_RadCal/W1P1M1/NDF_0/nd0.mat','nd0');
    [msLinearityFigHandles, ttfFigHandle] = analyze_light_logger_calibration_data(nd0,whichMSChannelsToPlot,msLinearityFigHandles,ttfFigHandle);
%}

    arguments 
        LightLoggerCalibrationData; % Struct that has both metadata and parsed readings from the experiment 
        whichMSChannelsToPlot = {[], []}; % TODO: Explain this : 
        msLinearityFigHandles = {}; % TODO: Explain this: 
        ttfFigHandle = []; % TODO: Explain this: 
    end 

    % First, extract the broad subfields of the calibration information and the parsed readings
    % from the input struct
    calibration_metadata = LightLoggerCalibrationData.metadata;
    parsed_readings = LightLoggerCalibrationData.readings;

    % 1. Analyze the MS linearity readings.
    msLinearityFigHandles = analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
                                                    parsed_readings.ms_linearity,...
                                                    calibration_metadata.cal,...
                                                    calibration_metadata.ndf,...
                                                    whichMSChannelsToPlot,...
                                                    msLinearityFigHandles...
                                                    );

    % 2. Analyze the Temporal Sensitivity readings
    if ~isempty(calibration_metadata.temporal_sensitivity.frequencies)
        ttfFigHandle = analyze_temporal_sensitivity_data(...
            calibration_metadata.temporal_sensitivity, ...
            parsed_readings.temporal_sensitivity, calibration_metadata.ndf,...
            ttfFigHandle);
    end

    % 3. Analyze the Phase fitting readings
    if ~isempty(calibration_metadata.phase_fitting.frequencies)
        temporal_offsets_secs = analyze_phase_fit_data(calibration_metadata.phase_fitting, parsed_readings.phase_fitting);

        % Report the temporal offset values and standard deviations
        fprintf('Mean and std of world - pupil temporal offset [ms]: %2.5f, %2.5f\n',1000*mean(temporal_offsets_secs('W-P')),1000*std(temporal_offsets_secs('W-P')))
        fprintf('Mean and std of world - AS temporal offset [ms]: %2.5f, %2.5f\n',1000*mean(temporal_offsets_secs('W-AS')),1000*std(temporal_offsets_secs('W-AS')))
        fprintf('Mean and std of world - TS temporal offset [ms]: %2.5f, %2.5f\n',1000*mean(temporal_offsets_secs('W-TS')),1000*std(temporal_offsets_secs('W-TS')))
    end

    % 4. Analyze the contrast gamma readings
    if ~isempty(calibration_metadata.contrast_gamma.frequencies)
        analyze_contrast_gamma_data(calibration_metadata.contrast_gamma, parsed_readings.contrast_gamma, calibration_metadata.ndf);
    end

end

