function avgSPDStruct = meanSPDs(spds, options)
% Compute mean and standard deviation of SPD results across subjects
%
% Syntax:
%   avgSPDStruct = meanSPDs(spds)
%   avgSPDStruct = meanSPDs(spds, options)
%
% Description:
%   Accepts a cell array of per-subject SPD result structures, each
%   containing one or more projection types (e.g., virtuallyFoveated,
%   justProjection). For each projection type, the function computes
%   element-wise means and standard deviations of the exponent maps,
%   variance maps, regional SPDs, median images, frequency vectors, and
%   best-fit line coefficients (slope, intercept, and fitted curves for
%   both center and periphery regions). The averaged results are returned
%   in a single struct and can optionally be saved to disk.
%
% Inputs:
%   spds                  - Cell array. Each element is a struct with
%                           fields indexed by projection type, where each
%                           projection type contains:
%                             .spd       - SPD result struct or filepath
%                             .best_fit  - Best-fit struct or filepath
%
% Optional key/value pairs:
%   output_path           - String. Path to save the averaged result as a
%                           .mat file. If empty (default), no file is
%                           written.
%
% Outputs:
%   avgSPDStruct          - Struct. Contains:
%                             .mean.(projection_type) - Mean values
%                             .std.(projection_type)  - Standard deviations
%                             .n                      - Number of input SPDs
%                           Each projection type subfield contains:
%                             .exponentMap, .varianceMap, .spdByRegion,
%                             .medianImage, .frq, .best_fit.center.*,
%                             .best_fit.periphery.*
%
% Examples:
%{
    spds = loadSPDS(src_dir);
    avgSPDStruct = meanSPDs(spds, "output_path", "mean_spds.mat");
%}

    arguments
        spds = {}
        options.output_path = "";
    end

    % The input argument of spds is a flat 
    % cell array of structs 
    % of the form {projection_type: SPDStruct (OR PATH TO SPD STRUCT): fields_to_average }
    spds = cell(spds);

    % The output avgSPDStruct should contain the following information 
    % should have the following form 
    % {mean: {projection_type: mean_of_fields_to_average}, std: {projection_type: std_of_fields_to_average}, n: numel(spds) }


    avgSPDStruct = struct();
    avgSPDStruct.mean = struct();
    avgSPDStruct.std = struct();
    avgSPDStruct.n = numel(spds);   

    % Retrieve the projection types from the first element in the spds 
    projection_types = fieldnames(spds{1}); 

    % Iterate over the projection types and compute means/stds for the
    % numeric SPD outputs independently within each projection family.
    for pp = 1:numel(projection_types)
        % Retrieve the name of the projection type we are loading 
        projection_type = projection_types{pp};

        % Initialize containers we will use to hold the data toaverage
        exponent_maps = [];
        variance_maps = [];
        spd_by_regions = [];
        median_images = [];
        frqs = [];
        best_fit_center_slope = [];
        best_fit_center_intercept = [];
        best_fit_center_log10_frequency = [];
        best_fit_center_log10_spd_fit = [];
        best_fit_center_frequency = [];
        best_fit_center_spd_fit = [];
        best_fit_periphery_slope = [];
        best_fit_periphery_intercept = [];
        best_fit_periphery_log10_frequency = [];
        best_fit_periphery_log10_spd_fit = [];
        best_fit_periphery_frequency = [];
        best_fit_periphery_spd_fit = [];

        % Iterate over the SPDs we have been provided
        for ii = 1:numel(spds)
            % Retrieve the current SPD
            current_spd = spds{ii};
            
            % Load in the SPD from the current spd
            loaded_projection_spd = iLoadProjectionEntry(current_spd.(projection_type).spd);
            loaded_projection_best_fit = iLoadBestFitEntry(current_spd.(projection_type).best_fit);
            

            % Save this spd's items 
            exponent_maps(:, :, ii) = loaded_projection_spd.exponentMap;
            variance_maps(:, :, ii) = loaded_projection_spd.varianceMap;
            spd_by_regions(:, :, :, ii) = loaded_projection_spd.spdByRegion;
            median_images(:, :, ii) = loaded_projection_spd.medianImage;
            frqs(:, ii) = loaded_projection_spd.frq(:);

            % Save this SPD's best-fit-line quantities for both the center
            % and periphery regions. Scalar coefficients are stacked across
            % columns, while vector-valued fit curves are aligned row-wise.
            best_fit_center_slope(1, ii) = loaded_projection_best_fit.center.slope;
            best_fit_center_intercept(1, ii) = loaded_projection_best_fit.center.intercept;
            best_fit_center_log10_frequency(:, ii) = loaded_projection_best_fit.center.log10_frequency(:);
            best_fit_center_log10_spd_fit(:, ii) = loaded_projection_best_fit.center.log10_spd_fit(:);
            best_fit_center_frequency(:, ii) = loaded_projection_best_fit.center.frequency(:);
            best_fit_center_spd_fit(:, ii) = loaded_projection_best_fit.center.spd_fit(:);

            best_fit_periphery_slope(1, ii) = loaded_projection_best_fit.periphery.slope;
            best_fit_periphery_intercept(1, ii) = loaded_projection_best_fit.periphery.intercept;
            best_fit_periphery_log10_frequency(:, ii) = loaded_projection_best_fit.periphery.log10_frequency(:);
            best_fit_periphery_log10_spd_fit(:, ii) = loaded_projection_best_fit.periphery.log10_spd_fit(:);
            best_fit_periphery_frequency(:, ii) = loaded_projection_best_fit.periphery.frequency(:);
            best_fit_periphery_spd_fit(:, ii) = loaded_projection_best_fit.periphery.spd_fit(:);
        end

        % Calculate the mean and STD of this projection type for all the SPDs
        avgSPDStruct.mean.(projection_type).exponentMap = mean(exponent_maps, 3, 'omitmissing');
        avgSPDStruct.mean.(projection_type).varianceMap = mean(variance_maps, 3, 'omitmissing');
        avgSPDStruct.mean.(projection_type).spdByRegion = mean(spd_by_regions, 4, 'omitmissing');
        avgSPDStruct.mean.(projection_type).medianImage = mean(median_images, 3, 'omitmissing');
        avgSPDStruct.mean.(projection_type).frq = mean(frqs, 2, 'omitmissing');

        avgSPDStruct.std.(projection_type).exponentMap = std(exponent_maps, 0, 3, 'omitmissing');
        avgSPDStruct.std.(projection_type).varianceMap = std(variance_maps, 0, 3, 'omitmissing');
        avgSPDStruct.std.(projection_type).spdByRegion = std(spd_by_regions, 0, 4, 'omitmissing');
        avgSPDStruct.std.(projection_type).medianImage = std(median_images, 0, 3, 'omitmissing');
        avgSPDStruct.std.(projection_type).frq = std(frqs, 0, 2, 'omitmissing');

        % Average the best-fit-line outputs for each spatial region and
        % store them alongside the main SPD summary fields.
        avgSPDStruct.mean.(projection_type).best_fit.center.slope = mean(best_fit_center_slope, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.center.intercept = mean(best_fit_center_intercept, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.center.log10_frequency = mean(best_fit_center_log10_frequency, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.center.log10_spd_fit = mean(best_fit_center_log10_spd_fit, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.center.frequency = mean(best_fit_center_frequency, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.center.spd_fit = mean(best_fit_center_spd_fit, 2, 'omitmissing');

        avgSPDStruct.mean.(projection_type).best_fit.periphery.slope = mean(best_fit_periphery_slope, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.periphery.intercept = mean(best_fit_periphery_intercept, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.periphery.log10_frequency = mean(best_fit_periphery_log10_frequency, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.periphery.log10_spd_fit = mean(best_fit_periphery_log10_spd_fit, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.periphery.frequency = mean(best_fit_periphery_frequency, 2, 'omitmissing');
        avgSPDStruct.mean.(projection_type).best_fit.periphery.spd_fit = mean(best_fit_periphery_spd_fit, 2, 'omitmissing');

        avgSPDStruct.std.(projection_type).best_fit.center.slope = std(best_fit_center_slope, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.center.intercept = std(best_fit_center_intercept, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.center.log10_frequency = std(best_fit_center_log10_frequency, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.center.log10_spd_fit = std(best_fit_center_log10_spd_fit, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.center.frequency = std(best_fit_center_frequency, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.center.spd_fit = std(best_fit_center_spd_fit, 0, 2, 'omitmissing');

        avgSPDStruct.std.(projection_type).best_fit.periphery.slope = std(best_fit_periphery_slope, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.periphery.intercept = std(best_fit_periphery_intercept, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.periphery.log10_frequency = std(best_fit_periphery_log10_frequency, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.periphery.log10_spd_fit = std(best_fit_periphery_log10_spd_fit, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.periphery.frequency = std(best_fit_periphery_frequency, 0, 2, 'omitmissing');
        avgSPDStruct.std.(projection_type).best_fit.periphery.spd_fit = std(best_fit_periphery_spd_fit, 0, 2, 'omitmissing');

        %%%% 
    end

    % If output path is not "", 
    % then we output to the target location
    if(~(options.output_path == ""))
        save(options.output_path, "avgSPDStruct")
    end 

end 


function loaded_spd = iLoadProjectionEntry(projection_entry)
% Internal helper to i load projection entry.
%
% Syntax:
%   loaded_spd = iLoadProjectionEntry(projection_entry)
%
% Description:
%   This local helper function internal helper to i load projection entry within its parent workflow.
% Inputs:
%   projection_entry         - Input used by the function.
%
% Outputs:
%   loaded_spd               - Output produced by the function.
%
% Examples:
%{
    % See meanSPDs.m for usage context.
%}

    if (isstring(projection_entry) || ischar(projection_entry))
        loaded_mat = load(projection_entry);

    % Otherwise, if we have a struct, then just go right ahead 
    elseif (isstruct(projection_entry))
        loaded_mat = projection_entry;

    % Otherwise, throw an error for an unsupported tpye 
    else
        error("Unsupported projection entry type: %s", class(projection_entry));
    end

    % Retrieve the activity name from the struct 
    activity_names = fieldnames(loaded_mat.activityData);
    
    % Splice just this field so it is the right format and return 
    loaded_spd = loaded_mat.activityData.(activity_names{1});
    
    return;
end


function loaded_best_fit = iLoadBestFitEntry(best_fit_entry)
% Internal helper to i load best fit entry.
%
% Syntax:
%   loaded_best_fit = iLoadBestFitEntry(best_fit_entry)
%
% Description:
%   This local helper function internal helper to i load best fit entry within its parent workflow.
% Inputs:
%   best_fit_entry           - Input used by the function.
%
% Outputs:
%   loaded_best_fit          - Output produced by the function.
%
% Examples:
%{
    % See meanSPDs.m for usage context.
%}

    if (isstring(best_fit_entry) || ischar(best_fit_entry))
        loaded_mat = load(best_fit_entry);

    % Otherwise, if we already have a struct, use it directly.
    elseif (isstruct(best_fit_entry))
        loaded_mat = best_fit_entry;

    % Otherwise, throw an error for an unsupported type.
    else
        error("Unsupported best-fit entry type: %s", class(best_fit_entry));
    end

    % Read in the loaded best fit field
    loaded_best_fit = loaded_mat.bestFit;
    return;

end
