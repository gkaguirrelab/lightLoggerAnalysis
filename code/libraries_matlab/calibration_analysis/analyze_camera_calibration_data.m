function camera_intrinsics = analyze_camera_calibration_data(path_to_images)

    arguments 
        path_to_images {mustBeText}; % Path to folder containing images used for Camera Calibration in MATLAB
    end 

    % First, let's load the Python helper library 
    addpath(getpref("lightLoggerAnalysis", "light_logger_libraries_matlab")); 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path")); 

    % Retrieve the calibration images in the form 
    % { (title, image)  }
    calibration_data_unconverted = world_util.load_calibration_images(path_to_images);





end 