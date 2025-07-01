function  analyze_ms_linearity_data(calibration_metadata, measurements)
% Analyze linearity calibration data collected from the MS
%
% Syntax:
%   analyze_ms_linearity_data(calibration_metadata, measurements)
%
% Description:
%   TODO 
%
% Inputs:
%   label                 - String. Optional notes to attach
%                           to the front of the saved filenames
%
%   dropbox_savedir       - String. The path to the directory
%                           on DropBox where the files will be
%                           output
%
% Outputs:
%
%   success               - Int. Returns 1 on success, 0 on failure.
%
% Examples:
%{

%}

    % Save the path to CombiExperiments. We will use this as a relative
    % path to find other files
    combiExperiments_path = getpref('lightLoggerAnalysis', 'combiExperiments_path');

    % Load the minispect SPDs
    spectral_sensitivity_map = containers.Map({'ASM7341', 'TSL2591'},...
                                              {fullfile(combiExperiments_path,'data','ASM7341_spectralSensitivity.mat');
                                              fullfile(combiExperiments_path,'data','TSL2591_spectralSensitivity.mat')...
                                              }...
                                            );

    % Create a map for the filters used to select good indices from the resulting curves for each chip
    as_chip_point_filter = @(x, y) and(and(~isinf(y), ~isinf(x)), y >= 0.25); % AS chip we want to exclude points in the mud
    ts_chip_point_filter= @(x, y) and(and(~isinf(y), ~isinf(x)), y < max(y)); % TS chip we want to exclude points that are saturated
    goodIdxFilterMap = containers.Map({'ASM7341', 'TSL2591'},...
                                      {as_chip_point_filter, ts_chip_point_filter}...
                                     );

                                    
    % Create a map for the limits for the chips' associated curves
    lim_map = containers.Map({'ASM7341', 'TSL2591'},...
                             {[-1, 5], [-1, 6]}...
                            );

    % Initialize a map between chips and the number of channels that they have
    n_channels_map = containers.Map({'ASM7341', 'TSL2591'},...
                                    {size(measurements{1,1}.M.v.AS, 2),...
                                     size(measurements{1,1}.M.v.TS, 2)...
                                    }...
                                   );


    % Get the background that was modified by the settings
    % scalar and shown to the MS
    background = calibration_metadata.background;

    % Retrieve the list of background scalars
    background_scalars = calibration_metadata.background_scalars;

    % Determine the number of settings levels that were exposed 
    n_settings_levels = numel(background_scalars);  

    % Save a map between the predicted and measured counts 
    % across NDF levels 
    predicted_measured_map = containers.Map( ...
                                            {'ASM7341', 'TSL2591'}, ...
                                            { ...
                                                {zeros(0, n_channels_map('ASM7341')), zeros(0, n_channels_map('ASM7341'))}, ...
                                                {zeros(0, n_channels_map('TSL2591')), zeros(0, n_channels_map('TSL2591'))} ...
                                            } ...
                                          );

    % First, let's iterate over the NDF levels 
    for nn = 1:numel(calibration_metadata.NDFs)
        % Retrieve the current NDF 
        NDF = calibration_metadata.NDFs(nn); 

        % Let's retrieve the cal file for this NDF
        cal = calibration_metadata.cal_files{nn};

        % Get the source from the cal file, as we need this to resample the
        % detector spectral sensitivity functions
        sourceS = cal.rawData.S; % This is the baseCal adjusted for transmitence

        % Extract information regarding the light source that was used to
        % calibrate the minispect
        sourceP_abs = cal.processedData.P_device;

        % Retrieve the wavelengths
        wls = SToWls(sourceS);

        % Reformat the minispect SPDs to be in the space of 
        % the source SPDs 
        minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS);

        % Next, let's iterate over the chips at this NDF level 
        chips = keys(spectral_sensitivity_map);
        for cc = 1:numel(chips)
            % Retrieve the name of the chip we are analyzing
            chip = chips{cc};

            % Initialize a summation variable for the detector counts
            % as we are going to average them over the reps
            sum_detector_counts = 0;

            % Initialize the predictedCounts variable to some value.
            % this is to have it in scope for use later.
            predictedCounts = 0;

            % Find the associated detectorP_rel for this chip
            detectorP_rel = minipspectP_rels_map(chip);

            % Grab the channels of the chip we are fitting;
            n_detector_channels = n_channels_map(chip);

            % Extract all of the counts from the sorted measurements
            % This gives you a settings x measurement x readings x channels
            % matrix
            detector_counts = extract_detector_counts(nn, measurements, chip);

            % Now, let's take the mean across measurements (the number of times)
            % we exposed this settings level as well as the mean across readings
            % the number of readings at each measurement. We may later be
            % interested in the variability across the set of measures.
            % This gives you a mat of settings level x count
            detector_counts = squeeze(mean(detector_counts, [2, 3]));

            % Plot the detector counts across settings levels
            figure ;
            title(sprintf("Averaged Counts by Settings Level | NDF %.3f | C: %s", NDF, chip));
            hold on ;
            for ch = 1:n_detector_channels
                plot(background_scalars, detector_counts(:, ch), "-x", 'DisplayName', sprintf("CH%d", ch));
            end

            % Label the plot
            xlabel("Settings Level");
            ylabel("Averged Count");

            % Set the limit for the mean detector counts as an unsigned 16 bit range 
            ylim([-5, 2^16+1000])

            % Show the legend for the plot
            legend show ;

            % Initialize some variables to hold loop results
            sphereSPDs = nan(n_settings_levels, sourceS(3));
            predictedCounts = nan(n_settings_levels, n_detector_channels);

            % Iterate over the primary steps, that is, the number of scalars
            for ss = 1:n_settings_levels
                % Get the source settings by multiplying background
                % by the scalar value at primaryStep kk
                source_settings = background * background_scalars(ss);

                % Derive the sphereSPD for this step in units of W/m2/sr/nm. We divide
                % by the nanometer sampling given in S to cast the units as nm, as
                % opposed to (e.g.) per 2 nm.
                sphereSPDs(ss,:) = ( (sourceP_abs*source_settings')/sourceS(2) );

                % Derive the prediction of the relative counts based upon the sphereSPD
                % and the minispectP_rel.
                predictedCounts(ss,:) = sphereSPDs(ss,:)*detectorP_rel;

            end % End n_primary steps

            % Find the limits for this chip
            limits = lim_map(chip);

            % Retrieve the filter function used to exclude points
            % from fitting LBF
            goodIdxFilter = goodIdxFilterMap(chip);

            % Retrieve the measured/predicted counts for this chip
            % at this NDF level
            measured_local = detector_counts;
            predicted_local = predictedCounts;

            % Append these to the across NDF measures 
            predicted_measured_global = predicted_measured_map(string(chip)); 
            predicted_measured_global{1} = [ predicted_measured_global{1} ; predicted_local]; 
            predicted_measured_global{2} = [ predicted_measured_global{2} ; measured_local]; 
            predicted_measured_map(chip) = predicted_measured_global; 

        end % Chip loop 

    end  % NDF loop 


    % Now, we can plot the results for each chip over all NDF levels 
    for cc = 1:numel(chips)
        % Retrieve the name of the chip 
        chip = chips{cc}; 

        % Retrieve the number of channels for this chip 
        n_detector_channels = n_channels_map(chip); 

        % Retrieve the predicted/measured values for this chip across NDF levels 
        predicted_measured = predicted_measured_map(chip); 
        predicted = predicted_measured{1}; 
        measured = predicted_measured{2};   

        % Generate a figure for this chip 
        [rows, cols] = find_min_figsize(n_detector_channels);
        figure; 
        t = tiledlayout(rows, cols);
        sgtitle(sprintf("%s Measured vs Predicted Counts", chip), 'FontWeight', 'bold', 'FontSize', 14); 

        % Loop across the channels of the chip and show predicted vs measured
        for ch = 1:n_detector_channels
            % Begin plotting on the tiles 
            nexttile; 

            % Draw a reference line
            plot([limits(1),limits(2)],[limits(1),limits(2)],':k', "DisplayName", "ReferenceLine");
            hold on;

            % Log transform the measured and predicted counts
            predicted_channel_readings = log10(predicted(:,ch));
            measured_channel_readings = log10(measured(:, ch));

            % Fit a linear model, but only to the "good" points (i.e., those
            % that are finite and not at the ceiling or floor. We also exclude
            % the points measured using the ND6 filter, as we do not have an
            % independent measure of the spectral transmittance of these.
            p = polyfit(predicted_channel_readings, measured_channel_readings,1);
            fitY = polyval(p, predicted_channel_readings);

            % Plot the original points
            plot(predicted_channel_readings, measured_channel_readings, '-x', 'DisplayName', 'Measured');
            plot(predicted_channel_readings, fitY, '-k', 'LineWidth', 1.5, 'DisplayName', "Fit");

            % Pretty up the plot
            xlim(limits);
            xlabel(sprintf('%s predicted counts [log]', chip));
            ylim(limits);
            ylabel(sprintf('%s measured counts [log]', chip));

            legend('Location','southeast');

            title(sprintf('channel %d, [slope intercept] = %2.2f, %2.2f',ch,p));

        end % Channel loop 
        hold off; 

    end 



end

% Local function to reformat the minispect SPDs to be in the space of 
% the source SPDs 
function minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS)
    % Next, let's iterate over the chips at this NDF level 
    % and reformat the minispect SPDs to be in the space of 
    % the source SPDs 
    minipspectP_rels_map = containers.Map();
    chips = keys(spectral_sensitivity_map); 
    for cc = 1:numel(chips)  
        % Retrieve the current chip 
        chip = chips{cc}; 

        miniSpectSPDPath = spectral_sensitivity_map(chip);
        load(miniSpectSPDPath,'T');
        minispectS = WlsToS(T.wl);
        minispectP_rel = T{:,2:end};

        detectorP_rel = [];
        for jj = 1:size(minispectP_rel,2)
            detectorP_rel(:,jj) = interp1(SToWls(minispectS),minispectP_rel(:,jj),SToWls(sourceS));
        end
        
        % Save the new SPDs 
        minipspectP_rels_map(chip) = detectorP_rel;
    end 

    return ; 
end 

% Local function to find the min square figsize required to plot data 
function [rows, cols] = find_min_figsize(num_plots)
    % Iterate over the ints between 1 and num_plots 
    for ii = 1:num_plots
        rows = ii; 
        cols = ii; 

        % Determine if we have reached the target 
        if(rows * cols >= num_plots)
            return ; 
        end 

    end 

end 

% Local function to extract solely the values of a given MS
% sensor from the measurements
function counts_mat = extract_detector_counts(NDF_num, measurements, chip)
    % Extract some information about the experiment
    [num_NDF_levels, num_settings_levels, n_measures] = size(measurements);

    % First, let's find the max number of readings we have and the number 
    % of channels that we have so that we can allocate a matrix 
    max_num_readings = 0; 
    n_channels = 0; 
    for ss = 1:num_settings_levels
        for nn = 1:n_measures
            % Extract this measurement struct
            measurement = measurements{NDF_num, ss, nn};

            % Add the number of readings 
            if(chip == "ASM7341")
                max_num_readings = max(max_num_readings, size(measurement.M.v.AS, 1));
                
                % Save the total number of channels if we have not already
                if(n_channels == 0)
                    n_channels = size(measurement.M.v.AS, 2);
                end 

            else 
                max_num_readings = max(max_num_readings, size(measurement.M.v.TS, 1));

                 % Save the total number of channels if we have not already
                if(n_channels == 0)
                    n_channels = size(measurement.M.v.TS, 2);
                end 

            end 
        end 
    end 


    % Initialize a matrix to store the values for this NDF 
    % TODO: Change this, can't have nan. It messes things up. 
    counts_mat = nan(num_settings_levels, n_measures, max_num_readings, n_channels);

    % Next, we will go over each measurement and extract the channels
    % for the desired chip 
    for ss = 1:num_settings_levels
        for nn = 1:n_measures
            % Extract this measurement struct
            measurement = measurements{NDF_num, ss, nn};

            % Retrieve the readings from the MS 
            readings = 0; 
            if(chip == 'ASM7341')
                readings = measurement.M.v.AS; 
            else    
                readings = measurement.M.v.TS; 
            end 

            % Insert these readings into the matrix 
            counts_mat(ss, nn, 1:size(readings, 1), :) = readings;

        
        end
    end

    return ;

end