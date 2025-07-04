function camera_intrinsics = convert_camera_calibration_data(input_path, output_path)


    arguments 
        input_path {mustBeText}; % Path to folder containing images (as raw .npy/.blosc) used for Camera Calibration in MATLAB
        output_path {mustBeText}; % Path to the folder where the converted images will be output as .tiff files
    end 

    % First, let's load the Python helper library 
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab")); 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path")); 

    % Retrieve the calibration images in the form 
    % { (title, image)  }
    calibration_data = cellfun(@(x) image_title_tuple_to_matlab(x),...
                               cell(world_util.load_calibration_images(input_path)),...
                               'UniformOutput', false...
                              );

    % Calculate the size needed for the preview grid 
    preview_grid_size = calculate_preview_grid_size(numel(calibration_data));

    % Display the images before we send them to be calibrated  
    t = tiledlayout(preview_grid_size(1), preview_grid_size(2)); 
    sgtitle("Image Preview");
    for ii = 1:numel(calibration_data)
        % Move to the next tile 
        nexttile;

        % Retrieve the image and title 
        title_and_image = calibration_data{ii}; 
        [title_str, image] = title_and_image{:}; 

        % Now, let's extract only the name part 
        % of the filename 
        [~, filename, ~] = fileparts(title_str);

        % Display the image 
        imshow(image);
        title(sprintf("%s", title_str));  
    
        % Save the image to the output path 
        % in an uncompressed file format 
        output_filepath = fullfile(output_path, sprintf("%s.tiff", filename));
        imwrite(uint8(image), output_filepath, 'Compression', 'none'); 
    end 

end 

% Local function to convert the py.tuple (title, image)
% to pure MATLAB type 
function converted = image_title_tuple_to_matlab(tuple)
    % Convert the tuple into a cell array 
    converted = cell(tuple); 

    % Next, the tuple is in the form (title, image)
    % so convert these to string and double 
    converted{1} = char(converted{1});
    converted{2} = uint8(converted{2}); 

    return ; 

end 

% Local function to calculate the size grid needed 
% to display the preview images 
% NOTE: You can do this more efficiently with binary search,
%       just doing it the lazy way for implementation speed 
function grid_size = calculate_preview_grid_size(num_images)
    % Iterate over ints up to the size of the images 
    for ii = 1:num_images 
        % If i^2 >= num_images, we found the size 
        if(ii ^ 2 >= num_images)
            grid_size = [ii, ii];
            return ;  
        end 
    end 
end 