function world_ttf_data = calcIMX219RadiometricSensitivity(options)
    arguments 
        options.verbose logical = false; 
    end 
    verbose = options.verbose;

    % Import the Python utility libraries we wil luse 
    % Note that Pi util is from light logger, which should not 
    % be a dependency for this repo, but due to how integrated 
    % the function I need is, I include it for now and will fix this 
    % later 
    Pi_util = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 
    world_util = import_pyfile(getpref("lightLoggerAnalysis", "world_util_path"))


% Estimate IMX219 radiometric sensitivity by channel

    % We need to obtain the spectrum of a source light and the IMX219 chip
    % spectral sensitivity functions.


    % Step 1: First, we will explore a 3NDF range of the world TTF measurements 
    % NDF 1, 2, and 3

    % To do this, we will first construct the path to the TTF data 
    world_ttf_data_dir = fullfile(getpref("lightLoggerAnalysis", "dropboxBaseDir"), "FLIC_data", "LightLogger_RadCal", "W1P1M1", "world_ttf_data_noAGC");
    world_ttf_data_path = fullfile(world_ttf_data_dir, "world_ttf_data.mat"); 
    world_ttf_data = load(world_ttf_data_path).world_ttf_data; 

    % Retrieve the metadata and the readings 
    world_ttf_metadata = world_ttf_data.metadata.temporal_sensitivity; 
    world_ttf_readings = world_ttf_data.readings.temporal_sensitivity; 

    % Retrieve the n_measures, frequencies and contrasts 
    n_measures = world_ttf_metadata.n_measures; 
    contrast_levels = world_ttf_metadata.contrast_levels; 
    frequencies = world_ttf_metadata.frequencies;

    % Let's iterate through the desired NDF range 
    for nn = 2:4 % Index 1 0 is 0 NDF so skip it 
        % Retrieve the NDF 
        % and the cal file for this NDF 
        NDF = world_ttf_metadata.NDFs(nn); 
        cal = world_ttf_metadata.cal_files{nn};
        if(verbose)
            fprintf("Processing NDF: %.2f\n", NDF);
        end 

        % Iterate over the 3 measurements we made
        for mm = 1:n_measures
            if(verbose)
                fprintf("\tProcessing Measure: %d\n", mm);
            end 
            
            % Retrieve the contrast level orders used for this 
            % NDF level and measurement 
            contrast_order = world_ttf_metadata.contrast_levels_orders(nn, mm);

            % Iterate over the contrast levels 
            for cc = 1:numel(contrast_order)
                % Retrieve the contrast level used for this measurement 
                contrast = contrast_levels(contrast_order(cc));
                if(verbose)
                    fprintf("\tProcessing contrast: %.2f\n", contrast);
                end 

                % Find the frequencies order for this contrast 
                frequencies_order = world_ttf_metadata.frequencies_orders(nn, mm, cc, :);

                % Iterate over the frequencies 
                for ff = 1:numel(frequencies)
                    frequency = frequencies(frequencies_order(ff)); 
                    if(verbose)
                        fprintf("\tProcessing frequency: %.2f\n", frequency);
                    end 

                    % Now, let's find the path to this recording folder
                    measurement_folder_path = fullfile(world_ttf_data_dir, sprintf("NDF%f", NDF), sprintf("TemporalSensitivity_%dcontrastIdx_%dfreqIdx_%dmeasurementIdx", cc, ff, mm));
                    if(~exist(measurement_folder_path, "dir"))
                        error(sprintf("Path: %s does not exist", measurement_folder_path)); 
                    end 
                    
                    % Now, first we need to parse the blosc file into chunks. This currently 
                    % uses something from the light logger repo, which technically should not be 
                    % a requirement for this repository, but it is so deeply engrained there 
                    % that I have not yet moved things yet
                    chunks = Pi_util.parse_chunks(measurement_folder_path, true, false, true, true)
                    mean_RGB = world_util.calculate_color_weights(chunks);
                    

                end 


            end     
            


        end     
    





    end



    %{
    % Here we load a file that Zach made that contains a calibration file and
    % derive the SPD of the light source. This is all a bit idiosyncratic, but
    % that's OK for this one-time measurement
    load("/Users/aguirre/Desktop/world_ttf_data.mat");

    % We will use a particular measurement (NDF1)
    measIdx = 2; % This is the 2nd measurement which was NDF=1.
    cal = world_ttf_data.metadata.temporal_sensitivity.cal_files{measIdx};

    % This cal file had already been adjusted for the presence of the NDF. We
    % need to account for the contrast that was used to present the stimulus.
    % We assume that the gamma correction was in place, so that the contrast
    % may be applied in a linear manner to the spd
    contrast = world_ttf_data.metadata.temporal_sensitivity.contrast_levels;
    sourceSPD = sum(cal.processedData.P_device,2) * contrast;

    % Get the wavelength support for the source
    sourceS = cal.rawData.S;

    % Load the IMX219 spectral sensitivity functions
    filePath = fullfile(tbLocateProjectSilent('lightLoggerAnalysis'),'data','IMX219_spectralSensitivity.mat');
    load(filePath,'T');
    sensorS = WlsToS(T.wls);

    % Interpolate the source SPD to be in the same domain as the sensor
    % sensitivity functions. Need to correct for the change from 2 to 1 nm
    % sampling
    sourceSPDInterp = interp1(SToWls(sourceS),sourceSPD,SToWls(sensorS),"linear","extrap");
    sourceSPDInterp = sourceSPDInterp * (sensorS(2) / sourceS(2));

    % Obtain the weights for the three channels
    channelLabels = {'blue','green','red'};
    for cc = 1:length(channelLabels)
        sensorSPD = T.(channelLabels{cc});
        weightBGR(cc) = sourceSPDInterp' * sensorSPD;
    end

    % This is an example set of observed weights from a video recording of this
    % source spectrum
    observedBGR = [ 1.14006268, 0.7269785 , 1]; 

    % We wish to calculate the factors by which we need to multiple the
    % observed BGR weights to result in values that match the predicted BGR
    % weights:
    correctionBGR = weightBGR ./ observedBGR;

    % The absolute value of the correctionBGR vector is not important; instead
    % we care about the relative correction weights. Set the mean of the
    % correction weights to unity.
    correctionBGR = correctionBGR / mean(correctionBGR);
    %}
end

