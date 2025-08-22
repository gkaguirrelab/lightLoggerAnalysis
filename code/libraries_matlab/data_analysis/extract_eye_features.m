function eye_features = extract_eye_features(video, is_grayscale, visualize_results)
% Extract eye features from a pupil video
%
% Syntax:
%   eye_features = extract_eye_features(video, is_grayscale, visualize_results)
%t
% Description:
%   TODO 
%
% Inputs:
%   video                 - tx480x640 array of 8 bit unsigned integers.
%                           This is a "chunk" of the world camera video
%   is_grayscale          - Logical. TODO
%   visualize_results     - Logical. TODO 
%
% Outputs:
%   eye_features          - Cell. The spectral power density in units of
% 
% Examples:
%{
    fps = 120;
    num_frames = 30 * fps;
    v = rand(num_frames, 480, 640)*255;
    [spd, frq] = calcTemporalSPD(v, fps);
    plotSPD(spd, frq);
%}


    arguments
        video; % Video, either path to video or an array representing its frames 
        is_grayscale {mustBeNumericOrLogical} = true;
        visualize_results {mustBeNumericOrLogical} = false; 
    end

    % First, we will import the Python utility library
    extract_eye_features_lib = import_pyfile(fullfile(getpref("lightLoggerAnalysis", "light_logger_analysis_libraries_python_path"), "extract_eye_features.py"));
    
    % Next, call the Python helper function to extract the feature 
    eye_features_py = cell(extract_eye_features_lib.extract_eye_features(video, is_grayscale, visualize_results)); 

    % Now, convert the python results into MATLAB
    eye_features = cellfun(@(x) convert_feature_dict(x), eye_features_py, 'UniformOutput', false); 

end

% Local function to convert a Python feature dict to MATLAB
function converted = convert_feature_dict(feature_dict_py)
    % First, convert the feature dict to a struct 
    converted = struct(feature_dict_py); 

    % Next, we will need to convert the fields of the struct. 
    % Let's first convert the easy ones 
    converted.location = double(converted.location); 

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
        converted.(field_name) = field; 
    end 

end 