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
    

    % If output path is not "", 
    % then we output to the target location
    if(~(options.output_path == ""))
        save(output_path, "avgSPDStruct")

    end 


end 