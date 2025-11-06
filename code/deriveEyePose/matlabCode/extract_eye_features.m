function eye_features = extract_eye_features(video, is_grayscale, visualize_results, pupil_feature_extraction_method)
% Extract eye features from a pupil video
%
% Syntax:
%   eye_features = extract_eye_features(video, is_grayscale, visualize_results, pupil_feature_extraction_method)
%
% Description:
%   Given a path to a video or an array of the frames of a video, 
%   extract eye features from each frame (pupil position, gaze angle, etc).  
%
% Inputs:
%   video                           - Text or Double. Path to the video to analyze 
%                                     or the video as an array of frames. 
%   is_grayscale                    - Logical. If the video is grayscale 
%
%   visualize_results               - Logical. Display plots of the tracked 
%                                     pupil per frame. 
%   
%   pupil_feature_extraction_method - String. Either pupil-labs or pylids
%                                     to use to extract pupil features
%
% Outputs:
%   eye_features          - Cell. Extarcted eye features per frame
% 
% Examples:
%{
    video = "path_to_video"; 
    is_grayscale = true; 
    visualize_results = false; 
    pupil_feature_extraction_method = "pylids"; 
    eye_features = extract_eye_features(video, is_grayscale, visualize_results, pupil_feature_extraction_method)
%}

    arguments
        video; % Video, either path to video or an array representing its frames 
        is_grayscale {mustBeNumericOrLogical} = true; % Whether or not the video is grayscale
        visualize_results {mustBeNumericOrLogical} = false; % Whether or not to generate plots of each frame overlaid with eyelid/pupil
        pupil_feature_extraction_method {mustBeText} = "pylids"; % Which pupil feature extraction method to use
    end

    % First, we will import the Python utility library
    extract_eye_features_lib = import_pyfile(fullfile(getpref("lightLoggerAnalysis", "light_logger_analysis_libraries_python_path"), "extract_eye_features.py"));
    
    % Next, call the Python helper function to extract the feature 
    eye_features_py = cell(extract_eye_features_lib.extract_eye_features(video, is_grayscale, visualize_results, pupil_feature_extraction_method)); 

    % Now, convert the python results into MATLAB
    eye_features = cellfun(@(x) convert_feature_dict(x, pupil_feature_extraction_method), eye_features_py, 'UniformOutput', false); 

    return ; 
end

% Local function to convert a Python feature dict to MATLAB
function converted = convert_feature_dict(feature_dict_py, pupil_feature_extraction_method)
    % First, convert the feature dict to a struct 
    converted = struct(feature_dict_py); 

    % Convert the pupil field of the struct 
    converted.pupil = convert_pupil_features(converted.pupil); 

    % Convert the eyelid field of the struct 
    converted.pupil = convert_eyelids_features(converted.eyelids); 

    return ; 
end 


% Local function to convert the pupil subdictionary 
function converted_pupil_features = convert_pupil_features(pupil_features, pupil_feature_extraction_method)
    % Initialize return variable 
    converted_pupil_features = struct(pupil_features); 

    % Parse the dict according to the method that was used to generate it
    if(pupil_feature_extraction_method == "pupil-labs")
        % Let's first convert the easy ones 
        converted_pupil_features.location = double(converted_pupil_features.location); 

        % Next, let's convert the others (which are all dicts themselves)
        fields_to_convert = {"sphere", "projected_sphere", "circle_3d", "ellipse"};
        for ff = 1:numel(fields_to_convert)
            % Retrieve the field name 
            field_name = fields_to_convert{ff}; 

            % Retrive the field (which is a subdict)
            field = struct(converted.(field_name));
            subfields = fieldnames(field);
            for sf = 1:numel(subfields)
                % Retrieve the subfield name 
                subfield_name = subfields{sf}; 
                
                % Convert the value to double 
                field.(subfield_name) = double(field.(subfield_name)); 
            end 

            % Save the converted field 
            converted_pupil_features.(field_name) = field; 
        end 
    
    else    
        % Gather the subfield names 
        subfields = fieldnames(converted_pupil_features); 

        % Iterate over the subfield names and convert all (np.ndarray) to double array 
        for sf = 1:numel(subfields)
            % Retrieve the subfield name 
            subfield_name = subfields{sf}; 

            % Convert the subfield 
            converted_pupil_features.(subfield) = double(converted_pupil_features.(subfield));
        end 
    end

    return  ; 
end 

% Local function to convert the eyelids subdictionary 
function converted_eyelid_features = convert_eyelid_features(eyelid_features)
    % Initialize return variable 
    converted_eyelid_features = struct(eyelid_features); 

    % Iterate over the the subfields of this new struct and convert all np.ndarrays 
    % to double 
    subfields = fieldnames(converted_eyelid_features);
    for sf = 1:numel(subfields)
        % Retrieve the subfield name 
        subfield_name = subfields{sf}; 

        % Convert the subfield np.ndarray to double array 
        convert_eyelid_features.(subfield_name) = double(converted_eyelid_features.(subfield_name));
    end 

    return ; 
end 