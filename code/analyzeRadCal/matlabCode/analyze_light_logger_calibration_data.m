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
%   opts.export_eps                 - Logical. If true, also export EPS
%                                     versions of each figure.
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
        opts.export_eps logical = true;
    end 

    figureHandles = {};

    % First, extract the broad subfields of the calibration information and the parsed readings
    % from the input struct
    calibration_metadata = light_logger_calibration_data.metadata;
    parsed_readings = light_logger_calibration_data.readings;

    % 1. Analyze the MS linearity readings if there are any to analyze.
    if(isfield(calibration_metadata, "ms_linearity") && numel(calibration_metadata.ms_linearity.NDFs) > 0) 
        figureHandles = [figureHandles; ...
            analyze_ms_linearity_data(calibration_metadata.ms_linearity,...
                                      parsed_readings.ms_linearity,...
                                      "plotSettingLevel", true...
                                     )];
    end                             


    % 2. Analyze the Temporal Sensitivity readings
    if(isfield(calibration_metadata, "temporal_sensitivity") && numel(calibration_metadata.temporal_sensitivity.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_temporal_sensitivity_data(calibration_metadata.temporal_sensitivity, ...
                                              parsed_readings.temporal_sensitivity...
                                             )];
    end

    % 3. Analyze the Phase fitting readings
    if(isfield(calibration_metadata, "phase_fitting") && numel(calibration_metadata.phase_fitting.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_phase_fit_data(calibration_metadata.phase_fitting,...
                                   parsed_readings.phase_fitting...
                                  )];
    end

    % 4. Analyze the contrast gamma readings
    if(isfield(calibration_metadata, "contrast_gamma") && numel(calibration_metadata.contrast_gamma.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_contrast_gamma_data(calibration_metadata.contrast_gamma,...
                                        parsed_readings.contrast_gamma...
                                       )];
    end

    % 5. Analyze the world camera linearity readings
    if(isfield(calibration_metadata, "world_linearity") && numel(calibration_metadata.world_linearity.NDFs) > 0)
        figureHandles = [figureHandles; ...
            analyze_camera_linearity_data(calibration_metadata.world_linearity, ...
                                          parsed_readings.world_linearity...
                                         )];
    end

    if(~islogical(opts.figure_output_path))
        export_figure_handles(figureHandles, opts.figure_output_path, opts.export_eps);
    end

end

function export_figure_handles(figureHandles, figure_output_path, export_eps)
    if(~exist(figure_output_path, "dir"))
        mkdir(figure_output_path);
    end

    pdf_output_path = fullfile(figure_output_path, "pdfs");

    if(~exist(pdf_output_path, "dir"))
        mkdir(pdf_output_path);
    end

    if(export_eps)
        eps_output_path = fullfile(figure_output_path, "eps");
        if(~exist(eps_output_path, "dir"))
            mkdir(eps_output_path);
        end
    end

    for ii = 1:numel(figureHandles)
        figHandle = figureHandles{ii};
        if(isempty(figHandle) || ~isgraphics(figHandle))
            continue;
        end

        if(is_figure_blank(figHandle))
            if(isgraphics(figHandle))
                close(figHandle);
            end
            continue;
        end

        figName = string(figHandle.Name);
        if(strlength(strtrim(figName)) == 0)
            figName = sprintf("figure_%02d", ii);
        end

        figName = regexprep(figName, '[^\w\d-]+', '_');
        figName = regexprep(figName, '_+', '_');
        pdf_export_path = fullfile(pdf_output_path, sprintf("%02d_%s.pdf", ii, figName));

        pdf_export_succeeded = false;

        try
            exportgraphics(figHandle, pdf_export_path, 'ContentType', 'vector');
            pdf_export_succeeded = exist(pdf_export_path, "file") ~= 0;
        catch
            exportapp(figHandle, pdf_export_path);
            pdf_export_succeeded = exist(pdf_export_path, "file") ~= 0;
        end

        if(~pdf_export_succeeded)
            error("Failed to export PDF for figure %d (%s)", ii, figName);
        end

        if(export_eps)
            eps_export_path = fullfile(eps_output_path, sprintf("%02d_%s.eps", ii, figName));
            eps_export_succeeded = false;

            try
                exportgraphics(figHandle, eps_export_path, 'ContentType', 'vector');
                eps_export_succeeded = exist(eps_export_path, "file") ~= 0;
            catch
                print(figHandle, eps_export_path, '-depsc', '-painters');
                eps_export_succeeded = exist(eps_export_path, "file") ~= 0;
            end

            if(~eps_export_succeeded)
                error("Failed to export EPS for figure %d (%s)", ii, figName);
            end
        end

        if(isgraphics(figHandle))
            close(figHandle);
        end
    end
end

function tf = is_figure_blank(figHandle)
    tf = true;

    axesHandles = findall(figHandle, 'Type', 'axes');
    for ax = axesHandles'
        if(~isempty(allchild(ax)))
            tf = false;
            return;
        end
    end
end
