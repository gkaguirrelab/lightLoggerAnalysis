function avgSPDStruct = meanSPDs(spds, options)
    arguments 
        spds = {}
        options.output_path = ""; 
    end 

    % The input argument of spds is a flat 
    % cell array of structs 
    % of the form {projection_type: SPDStruct (OR PATH TO SPD STRUCT): fields_to_average }
    
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

        % Iterate over the SPDs we have been provided
        for ii = 1:numel(spds)
            % Retrieve the current SPD
            current_spd = spds{ii};
            
            % Load in the SPD from the current spd
            loaded_projection_spd = iLoadProjectionEntry(current_spd.(projection_type).spd);
   

            % Save this spd's items 
            exponent_maps(:, :, ii) = loaded_projection_spd.exponentMap;
            variance_maps(:, :, ii) = loaded_projection_spd.varianceMap;
            spd_by_regions(:, :, :, ii) = loaded_projection_spd.spdByRegion;
            median_images(:, :, ii) = loaded_projection_spd.medianImage;
            frqs(:, ii) = loaded_projection_spd.frq(:);
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

        %%%% 
    end

    % If output path is not "", 
    % then we output to the target location
    if(~(options.output_path == ""))
        save(options.output_path, "avgSPDStruct")

    end 


end 


function loaded_spd = iLoadProjectionEntry(projection_entry)
% Normalize a projection leaf into the single-activity SPD struct that
% contains exponentMap / varianceMap / spdByRegion / frq / medianImage.

    % If we passed in a string path to the spd, load it in 
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
