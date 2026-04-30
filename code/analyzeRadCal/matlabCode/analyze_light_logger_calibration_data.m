function figureHandles = analyze_light_logger_calibration_data(light_logger_calibration_data, opts)
% Analyze the results of a light logger calibration measurement (post-conversion)
%
% Syntax:
%  figureHandles = analyze_light_logger_calibration_data(light_logger_calibration_data)
%  figureHandles = analyze_light_logger_calibration_data(light_logger_calibration_data, opts)
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
%   opts.figure_output_path         - false or path. If not logical, export
%                                     all generated figures to this folder.
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
        opts.figure_output_path = false;
    end 

    figureHandles = {};

    % First, extract the broad subfields of the calibration information and the parsed readings
    % from the input struct
    calibration_metadata = light_logger_calibration_data.metadata;
    parsed_readings = light_logger_calibration_data.readings;

    % 1. Analyze the MS linearity readings if there are any to analyze.
    if(numel(calibration_metadata.ms_linearity.NDFs) > 0) 
        figureHandles = [figureHandles; ...
            analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
                                      parsed_readings.ms_linearity,...
                                      "plotSettingLevel", true...
                                     )];
    end                             


    % 2. Analyze the Temporal Sensitivity readings
    if(numel(calibration_metadata.temporal_sensitivity.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_temporal_sensitivity_data(calibration_metadata.temporal_sensitivity, ...
                                              parsed_readings.temporal_sensitivity...
                                             )];
    end

    % 3. Analyze the Phase fitting readings
    if(numel(calibration_metadata.phase_fitting.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_phase_fit_data(calibration_metadata.phase_fitting,...
                                   parsed_readings.phase_fitting...
                                  )];
    end

    % 4. Analyze the contrast gamma readings
    if(numel(calibration_metadata.contrast_gamma.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_contrast_gamma_data(calibration_metadata.contrast_gamma,...
                                        parsed_readings.contrast_gamma...
                                       )];
    end

    if(~islogical(opts.figure_output_path))
        export_figure_handles(figureHandles, opts.figure_output_path);
    end

end

function export_figure_handles(figureHandles, figure_output_path)
    if(~exist(figure_output_path, "dir"))
        mkdir(figure_output_path);
    end

    pdf_output_path = fullfile(figure_output_path, "pdfs");
    eps_output_path = fullfile(figure_output_path, "eps");

    if(~exist(pdf_output_path, "dir"))
        mkdir(pdf_output_path);
    end

    if(~exist(eps_output_path, "dir"))
        mkdir(eps_output_path);
    end

    for ii = 1:numel(figureHandles)
        figHandle = figureHandles{ii};
        if(isempty(figHandle) || ~isgraphics(figHandle))
            continue;
        end

        is_ui_figure = isa(figHandle, 'matlab.ui.Figure');

        figName = string(figHandle.Name);
        if(strlength(strtrim(figName)) == 0)
            figName = sprintf("figure_%02d", ii);
        end

        figName = regexprep(figName, '[^\w\d-]+', '_');
        figName = regexprep(figName, '_+', '_');
        pdf_export_path = fullfile(pdf_output_path, sprintf("%02d_%s.pdf", ii, figName));
        eps_export_path = fullfile(eps_output_path, sprintf("%02d_%s.eps", ii, figName));

        try
            if(is_ui_figure)
                exportapp(figHandle, pdf_export_path);
            else
                exportgraphics(figHandle, pdf_export_path, 'ContentType', 'vector');
            end
        catch
            exportapp(figHandle, pdf_export_path);
        end

        if(~is_ui_figure)
            try
                exportgraphics(figHandle, eps_export_path, 'ContentType', 'vector');
            catch
                print(figHandle, eps_export_path, '-depsc', '-painters');
            end
        end

        if(isgraphics(figHandle))
            close(figHandle);
        end
    end
end
