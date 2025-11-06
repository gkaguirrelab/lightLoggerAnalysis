function calibration_metadata = load_calibration_metadata(path_to_metadata)
    % Open the file 
    calibration_metadata_struct_file = fopen(path_to_metadata, 'rb'); % Open the file whose bytes represent the CalibrationData struct 
    
    % Read the contents of the file as raw bytes 
    calibration_metadata_struct_bytes = fread(calibration_metadata_struct_file, Inf, '*uint8'); % Load in the bytes from that file
    
    % Parse back into a struct 
    calibration_metadata = getArrayFromByteStream(calibration_metadata_struct_bytes); % Parse the bytes back into a struct
    
    % Close the file we opened
    fclose(calibration_metadata_struct_file); % Close the file

    return ; 
end 