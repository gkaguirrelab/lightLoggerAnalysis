function camera_intrinsics = analyze_camera_calibration_data(path_to_images)

    arguments 
        path_to_images {mustBeText}; % Path to folder containing images used for Camera Calibration in MATLAB
    end 

    % First, let's load the Python helper library 
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab")); 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path")); 

    % Retrieve the calibration images in the form 
    % { (title, image)  }
    calibration_data = cellfun(@(x) image_title_tuple_to_matlab(x),...
                               cell(world_util.load_calibration_images(path_to_images)),...
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

        % Display the image 
        imshow(image);
        title(sprintf("%s", title_str));  
    
    end 

    % Extract only the images from the array 

    % Calculate the camera intrinstics



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