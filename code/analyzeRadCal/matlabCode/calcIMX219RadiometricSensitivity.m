function all_corrections_struct = calcIMX219RadiometricSensitivity(options)
    arguments 
        options.verbose logical = false; 
    end 
    verbose = options.verbose;


    % Estimate IMX219 radiometric sensitivity by channel

    % We need to obtain the spectrum of a source light and the IMX219 chip
    % spectral sensitivity functions.

    camera_spectral_sensitivity_filepath = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'),'data','IMX219_spectralSensitivity.mat');
    load(camera_spectral_sensitivity_filepath,'T');
    sensorS = WlsToS(T.wls);


    % Step 1: First, we will explore a 3NDF range of the world TTF measurements 
    % NDF 1, 2, and 3

    % To do this, we will first construct the path to the TTF data and the contrast gamma data
    world_ttf_data_dir = fullfile(getpref("lightLoggerAnalysis", "dropboxBaseDir"), "FLIC_data", "LightLogger_RadCal", "W1P1M1", "world_ttf_data_noAGC");
    world_ttf_data_path = fullfile(world_ttf_data_dir, "world_ttf_data.mat"); 
    world_ttf_data = load(world_ttf_data_path).world_ttf_data; 

    contrast_gamma_data_dir = fullfile(getpref("lightLoggerAnalysis", "dropboxBaseDir"), "FLIC_data", "LightLogger_RadCal", "W1P1M1", "contrast_gamma_data");
    contrast_gamma_data_path = fullfile(contrast_gamma_data_dir, "contrast_gamma_data.mat"); 
    contrast_gamma_data = load(contrast_gamma_data_path).contrast_gamma_data; 

    % Retrieve the metadata and the readings 
    world_ttf_metadata = world_ttf_data.metadata.temporal_sensitivity; 
    contrast_gamma_metadata = contrast_gamma_data.metadata.contrast_gamma;

    %ttf_all_corrections = calculate_corrections_by_experiment("TemporalSensitivity", world_ttf_data_dir, world_ttf_metadata, sensorS, T, verbose)
    contrast_gamma_all_corrections = calculate_corrections_by_experiment("ContrastGamma_TemporalSensitivity", contrast_gamma_data_dir, contrast_gamma_metadata, sensorS, T, verbose);
    
    all_corrections_struct.ttf_corrections = ttf_all_corrections; 
    all_corrections_struct.contrast_gamma_corrections = contrast_gamma_all_corrections;



end



function all_corrections_struct = calculate_corrections_by_experiment(experiment_class, data_dir, metadata, sensorS, T, verbose)
    % Import the Python utility libraries we wil luse 
    % Note that Pi util is from light logger, which should not 
    % be a dependency for this repo, but due to how integrated 
    % the function I need is, I include it for now and will fix this 
    % later 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"));

    % Retrieve the n_measures, frequencies and contrasts 
    NDFs = metadata.NDFs; 
    n_measures = metadata.n_measures; 
    contrast_levels = metadata.contrast_levels; 
    frequencies = metadata.frequencies;

    % Initialize the return list of BGR corrections 
    all_corrections_struct = struct; 
    all_corrections = nan(numel(NDFs), n_measures, numel(contrast_levels), numel(frequencies), 3);
    all_corrections_source_paths =  nan(numel(NDFs), n_measures, numel(contrast_levels), numel(frequencies));

    % Let's iterate through the desired NDF range 
    if(numel(NDFs) == 1)
        start = 1;
        endpoint = 1; 
    else 
        start = 4; 
        endpoint = 4; 
    end 

    for nn = start:endpoint % Index 1 0 is 0 NDF so skip it 
        % Retrieve the NDF 
        % and the cal file for this NDF 
        NDF = NDFs(nn); 
        if(verbose)
            fprintf("Processing NDF num: %d | NDF: %.2f\n", nn, NDF);
        end 

        cal = metadata.cal_files{nn};

        % Iterate over the 3 measurements we made
        for mm = 1:n_measures
            if(verbose)
                fprintf("\tProcessing Measure: %d\n", mm);
            end 
            
            % Retrieve the contrast level orders used for this 
            % NDF level and measurement 
            contrast_order = metadata.contrast_levels_orders(nn, mm, :);

            % Iterate over the contrast levels 
            for cc = 1:numel(contrast_order)
                % Retrieve the contrast level used for this measurement 
                contrast_idx = contrast_order(cc); 
                contrast = contrast_levels(contrast_idx);
                if(verbose)
                    fprintf("\t\tProcessing contrast idx: %d | contrast: %.2f\n", contrast_idx, contrast);
                end 

                % This cal file had already been adjusted for the presence of the NDF. We
                % need to account for the contrast that was used to present the stimulus.
                % We assume that the gamma correction was in place, so that the contrast
                % may be applied in a linear manner to the spd
                sourceSPD = sum(cal.processedData.P_device,2) * contrast;

                % Get the wavelength support for the source
                sourceS = cal.rawData.S;
                
                % Interpolate the source SPD to be in the same domain as the sensor
                % sensitivity functions. Need to correct for the change from 2 to 1 nm
                % sampling
                sourceSPDInterp = interp1(SToWls(sourceS),sourceSPD,SToWls(sensorS),"linear","extrap");
                sourceSPDInterp = sourceSPDInterp * (sensorS(2) / sourceS(2));

                % Obtain the weights for the three channels
                channelLabels = {'blue','green','red'};
                for ch = 1:length(channelLabels)
                    sensorSPD = T.(channelLabels{ch});
                    weightBGR(ch) = sourceSPDInterp' * sensorSPD;
                end

                % Find the frequencies order for this contrast 
                frequencies_order = metadata.frequencies_orders(nn, mm, contrast_idx, :);

                % Iterate over the frequencies 
                for ff = 1:numel(frequencies)
                    frequency_idx = frequencies_order(ff);
                    frequency = frequencies(frequency_idx); 
                    if(verbose)
                        fprintf("\t\t\tProcessing frequency idx: %d | frequency: %.2f\n", frequency_idx, frequency);
                    end 

                    % Now, let's find the path to this recording folder
                    measurement_folder_path = fullfile(data_dir, sprintf("NDF%f", NDF), sprintf("%s_%dcontrastIdx_%dfreqIdx_%dmeasurementIdx", experiment_class, cc, ff, mm));
                    if(~exist(measurement_folder_path, "dir"))
                        error(sprintf("Path: %s does not exist", measurement_folder_path)); 
                    end 
                    if(verbose)
                        fprintf("\t\t\t\tusing path: %s\n", measurement_folder_path);
                    end     
                    
                    % Now, first we need to parse the blosc file into chunks. This currently 
                    % uses something from the light logger repo, which technically should not be 
                    % a requirement for this repository, but it is so deeply engrained there 
                    % that I have not yet moved things yet
                    chunks = Pi_util.parse_chunks(measurement_folder_path, true, false, true, true);
                    observed_RGB = double(world_util.calculate_color_weights(chunks, true));
                    observed_BGR = permute(observed_RGB, [3 2 1]);
                    %py.gc.collect(); 

                    % We wish to calculate the factors by which we need to multiple the
                    % observed BGR weights to result in values that match the predicted BGR
                    % weights:
                    correctionBGR = weightBGR ./ observed_BGR;

                    % The absolute value of the correctionBGR vector is not important; instead
                    % we care about the relative correction weights. Set the mean of the
                    % correction weights to unity.
                    correctionBGR = correctionBGR / mean(correctionBGR);
                    
                    if(verbose)
                        fprintf("\t\t\t\t\tCorrection BGR");
                        disp(correctionBGR)
                    end 

                    % Save this BGR correction and the source path for it 
                    all_corrections(nn, mm, contrast_idx, frequency_idx, :) = correctionBGR; 
                    all_corrections_source_paths(nn, mm, contrast_idx, frequency_idx) = measurement_folder_path;
                end 
            end     
        end     
    end

    % Save the results into the struct 
    all_corrections_struct.corrections = all_corrections; 
    all_corrections_struct.source_paths = all_corrections_source_paths;

    return; 

end 