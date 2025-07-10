function  analyze_ms_linearity_data(calibration_metadata, measurements)
% Analyze linearity calibration data collected from the MS
%
% Syntax:
%   analyze_ms_linearity_data(calibration_metadata, measurements)
%
% Description:
%   Analyze linearity calibration data collected from the MS. 
%   Illustrate several plots showing the raw sensor counts 
%   from all channels of each chip at various settings/NDF 
%   levels, as well as the linearity of these in comparison 
%   to the predicted counts at a given NDF level.  
%
% Inputs:
%   calibration_metadata  - Struct. ms_linearity substruct 
%                           of the light logger metadata 
%                           struct containing only 
%                           ms_linearity related metadata 
%
%   measurements         - Cell. Parsed recordings made from 
%                          the light logger and converted to MATLAB 
%                          type.  
%
% Outputs:
%
%   NONE             
%
% Examples:
%{

%}
    arguments 
        calibration_metadata; % Struct representing the metadata for the ms_linearity calibration measurement 
        measurements; % Parsed and converted recordings from the light logger 
    end 


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
                             {[-1, 5], [-6, 6]}...
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

    % Initialize a matrix of starting/ending indices 
    % for each NDF for each chip 
    NDF_start_end_map = containers.Map( ...
                                        {'ASM7341', 'TSL2591'}, ...
                                        { zeros(numel(calibration_metadata.NDFs), 2), ...
                                          zeros(numel(calibration_metadata.NDFs), 2) ...
                                        } ...
                                      );

    % Make a list of colors for each ND level for the conjoined plot
    colorList = [
        0, 0.4470, 0.7410;   % Blue
        0.8500, 0.3250, 0.0980;  % Orange
        0.9290, 0.6940, 0.1250;  % Yellow
        0.4940, 0.1840, 0.5560;  % Purple
        0.4660, 0.6740, 0.1880;  % Green
        0.3010, 0.7450, 0.9330;  % Light Blue
        0.6350, 0.0780, 0.1840   % Red
        ];                                      

    % add path to isetBIO CIE luminous efficiency function
    addpath('~/Documents/MATLAB/toolboxes/Psychtoolbox-3/Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles');
    load("T_CIE_Y2.mat");
    wls_CIE_Y2 = SToWls(S_CIE_Y2); % convert to wavelength
    
    % First, let's iterate over the chips 
    chips = keys(spectral_sensitivity_map);
    for cc = 1:numel(chips)
        % Retrieve the name of the chip we are analyzing
        chip = chips{cc};

        % Grab the channels of the chip we are fitting;
        n_detector_channels = n_channels_map(chip);

        % Retrieve the limits for this chip 
        limits = lim_map(chip); 

        % Create a tiledlayout figure we will use to show the counts by settings 
        % level for this chip across NDF levels 
        [rows, cols] = find_min_figsize(numel(calibration_metadata.NDFs)); 
        figure; 
        counts_by_NDF_tiled = tiledlayout(rows, cols); 
        title(counts_by_NDF_tiled, sprintf("Counts by NDF Level | C: %s", chip), 'FontWeight', 'Bold'); 

        % Create a tabbed figure with one tab per NDF level. 
        % Each tab will have all of a given chip's channels on it 
        % for that NDF 
        uif = uifigure('Name', sprintf("%s Channel Linearity Per NDF", chip)); 
        tab_group = uitabgroup(uif);  

        % Initialize a cell array that will hold the measured/predicted value by NDF 
        measured_predicted_by_NDF = {}; 

        % Next, iterate over the NDF levels
        for nn = 1:numel(calibration_metadata.NDFs)
            % Retrieve the current NDF 
            NDF = calibration_metadata.NDFs(nn); 

            % Creat a new tab for this NDF for measured/predicted plot  
            tab = uitab(tab_group, 'Title', sprintf("NDF %.2f", NDF)); 
            
            % Define a tiled layout that will go in this tab 
            [rows, cols] = find_min_figsize(n_channels_map(chip)); 
            measured_predicted_tiled = tiledlayout(tab, rows, cols); 

            % Let's retrieve the cal file for this NDF
            cal = calibration_metadata.cal_files{nn};

            % Get the source from the cal file, as we need this to resample the
            % detector spectral sensitivity functions
            % TO DO update code so rawData is available
            %sourceS = cal.rawData.S; % This is the baseCal adjusted for transmitence
            %in the meantime, hardcoded 
            sourceS = [380 2 201];

            % Extract information regarding the light source that was used to
            % calibrate the minispect
            sourceP_abs = cal.processedData.P_device;

            % Retrieve the wavelengths
            wls = SToWls(sourceS);

            % Reformat the minispect SPDs to be in the space of 
            % the source SPDs 
            minipspectP_rels_map = reformat_SPDs(spectral_sensitivity_map, sourceS);

            % Initialize a summation variable for the detector counts
            % as we are going to average them over the reps
            sum_detector_counts = 0;

            % Initialize the predictedCounts variable to some value.
            % this is to have it in scope for use later.
            predictedCounts = 0;

            % Find the associated detectorP_rel for this chip
            detectorP_rel = minipspectP_rels_map(chip);

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
            counts_ax = nexttile(counts_by_NDF_tiled) ;

            title(counts_ax, sprintf("Averaged Counts by Settings Level | NDF %.3f", NDF));
            hold(counts_ax, 'on'); 
            for ch = 1:n_detector_channels
                plot(counts_ax, background_scalars, detector_counts(:, ch), "-x", 'DisplayName', sprintf("Ch%d", ch));
            end
            
            % Label the plot
            xlabel(counts_ax, "Settings Level");
            ylabel(counts_ax, "Averged Count");

            % Show the legend for the plot
            legend(counts_ax, 'Location', 'best'); 
            hold(counts_ax, 'off');

            % Initialize variables to calculate the predicted counts
            sphereSPDs = nan(n_settings_levels, sourceS(3));
            predictedCounts = nan(n_settings_levels, n_detector_channels);

            %resample luminous efficacy function in same space as S (wavelengths)    
            T_CIE_Y2_resamp = interp1(wls_CIE_Y2, T_CIE_Y2, wls, 'linear', 0); 

            % Iterate over the primary steps, that is, the number of scalars
            % in radiance
            for ss = 1:n_settings_levels
                % Get the source settings by multiplying background
                % by the scalar value at primaryStep kk
                source_settings = background * background_scalars(ss);

                % Derive the sphereSPD for this step in units of W/m2/sr/nm. We divide
                % by the nanometer sampling given in S to cast the units as nm, as
                % opposed to (e.g.) per 2 nm.
                sphereSPDs(ss,:) = ( (sourceP_abs*source_settings')/sourceS(2) );

                %Calc illuminance
                %first, convert to irradiance
                sphereIrrad(ss,:) = sphereSPDs(ss,:).*pi;
                irradXCIE(ss,:) = sphereIrrad(ss,:) .* T_CIE_Y2_resamp';
                integratedIrradXCIE(ss) = sum(irradXCIE(ss,:));
                % Multiply this integrated value by the “maximum luminous efficacy” value,
                % which converts from Watts to lumens. This value is 683 lumen/Watt.
                %result is the illuminance in units of lux (lumen / m2).
                illum(ss) = integratedIrradXCIE(ss)*683;

                % Derive the prediction of the relative counts based upon the sphereSPD
                % and the minispectP_rel.
                predictedCounts(ss,:) = sphereSPDs(ss,:)*detectorP_rel;

            end % End n_primary steps

            % Retrieve the measured/predicted counts for this chip
            % at this NDF level
            measured = detector_counts;
            predicted = predictedCounts;

            % Next, we will plot the measured and predicted counts per channel
            for ch = 1:n_detector_channels
                % Retrieve the ax to plot on 
                measured_predicted_ax = nexttile(measured_predicted_tiled); 
                
                % Plot the measured vs predicted data 
                plot(measured_predicted_ax, log10((predicted(:, ch)).*pi), log10(measured(:, ch)), '-x', 'DisplayName', 'Data'); 
                hold(measured_predicted_ax, 'on'); 

                % Fit a linear model, but only to the "good" points (i.e., those
                % that are finite and not at the ceiling or floor. We also exclude
                % the points measured using the ND6 filter, as we do not have an
                % independent measure of the spectral transmittance of these.
                p = polyfit(log10(predicted(:, ch)), log10(measured(:, ch)), 1);
                fitY = polyval(p, log10(predicted(:, ch)));
                plot(measured_predicted_ax, (log10(predicted(:, ch)).*pi), fitY, '-r', 'DisplayName', 'Fit'); 

                % Add a reference line
                plot(measured_predicted_ax, [limits(1),limits(2)],[limits(1),limits(2)],':k', "DisplayName", "ReferenceLine");

                % Pretty up the plot
                xlim(measured_predicted_ax, limits);
                xlabel(measured_predicted_ax, sprintf('%s predicted irradiance [log]', chip));
                ylim(measured_predicted_ax, limits);
                ylabel(measured_predicted_ax, sprintf('%s measured counts [log]', chip));

                legend(measured_predicted_ax, 'Location','best');

                title(measured_predicted_ax, sprintf('channel %d, [slope intercept] = %2.2f, %2.2f',ch, p));
                
                hold(measured_predicted_ax, 'off');
            end

            % Save the measured and predicted at this NDF
            measured_predicted_by_NDF{nn} = {measured, predicted};
            predicted_Illum_by_NDF(cc,nn,:) = illum;

        end  % NDF loop

        % Now, plot linearity across all NDF levels for a given chip
        [rows, cols] = find_min_figsize(n_detector_channels);
        % figure;
        % across_NDF_figure = tiledlayout(rows, cols);
        % title(sprintf("Measured vs Predicted across NDF | C: %s", chip), 'FontWeight', 'Bold');

        
        % Iterate over the channels 
        for ch = 1:n_detector_channels
            figure;
            % Retrieve the axes to plot on
            %across_NDF_channel_ax = nexttile(across_NDF_figure);
            across_NDF_channel_ax = axes;
            hold(across_NDF_channel_ax, 'on');

            % Retrieve the values per NDF level
            NDF_data = zeros(1, numel(calibration_metadata.NDFs));

            % Plot the data 
            for nn = 1:numel(calibration_metadata.NDFs)
                NDF_measured_predicted = measured_predicted_by_NDF{nn}; 
                measured = NDF_measured_predicted{1};
                predicted = NDF_measured_predicted{2} * pi; % converts radiance to irradiance

                h = scatter(across_NDF_channel_ax,...
                            log10(predicted(:, ch)), log10(measured(:, ch)),...
                            'o','MarkerFaceColor', colorList(nn,:), 'DisplayName', sprintf("NDF%.1f", calibration_metadata.NDFs(nn))...
                            );
                            
                h.MarkerFaceAlpha = 0.2;
                h.MarkerEdgeAlpha = 0.2;
                hold(across_NDF_channel_ax, 'on'); 
            end 

            % Add a reference line
            plot(across_NDF_channel_ax, [limits(1),limits(2)],[limits(1),limits(2)],':k', "DisplayName", "ReferenceLine");

            % Label the plot 
            title(across_NDF_channel_ax, sprintf("Ch: %d", ch)); 
            
            xlim(across_NDF_channel_ax, limits);
            xlabel(across_NDF_channel_ax, sprintf('%s predicted irradiance [log]', chip));
            ylim(across_NDF_channel_ax, limits);
            ylabel(across_NDF_channel_ax, sprintf('%s measured counts [log]', chip));

            legend(across_NDF_channel_ax, 'Location', 'best');
            %Add a second x-axis for: Illuminance vs Measured Counts
            % Create overlay axes (top x-axis, same Y, transparent background)
            % illum_ax = axes('Position', get(across_NDF_channel_ax, 'Position'), ...
            %     'Color', 'none', ...
            %     'XAxisLocation', 'top', ...
            %     'YAxisLocation', 'right', ...
            %     'XColor', [0.2 0.2 1], ...
            %     'YColor', 'none',...
            %     'Box', 'off',...
            %     'Parent', gcf);
            % hold(illum_ax, 'on');

            % % Plot the same data using illuminance
            % for nn = 1:numel(calibration_metadata.NDFs)
            %     NDF_measured_predicted = measured_predicted_by_NDF{nn};
            %     measured = NDF_measured_predicted{1};
            %     p = scatter(illum_ax, log10(predicted_Illum_by_NDF(ch,nn,:)), log10(measured(:, ch)), ...
            %         'o','MarkerFaceColor', colorList(nn,:));
            %     p.MarkerFaceAlpha = 0.2;
            %     p.MarkerEdgeAlpha = 0.2;
            % end
            % 
            % %Set axis limits to match original Y
            % set(illum_ax, 'YLim', get(across_NDF_channel_ax, 'YLim'));
            % set(illum_ax, 'XLimMode', 'auto');  % independent x-limits OK
            % xlabel(illum_ax, 'log Illuminance [lux]');

                hold(across_NDF_channel_ax, 'off');
                % hold(illum_ax, 'off');


        end
    
    end % Chip loop
predicted_Illum_by_NDF(1,:,5)

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
    min_num_readings = inf; 
    n_channels = 0; 
    for ss = 1:num_settings_levels
        for nn = 1:n_measures
            % Extract this measurement struct
            measurement = measurements{NDF_num, ss, nn};

            % Add the number of readings 
            if(chip == "ASM7341")
                min_num_readings = min(min_num_readings, size(measurement.M.v.AS, 1));
                
                % Save the total number of channels if we have not already
                if(n_channels == 0)
                    n_channels = size(measurement.M.v.AS, 2);
                end 

            else 
                min_num_readings = min(min_num_readings, size(measurement.M.v.TS, 1));

                 % Save the total number of channels if we have not already
                if(n_channels == 0)
                    n_channels = size(measurement.M.v.TS, 2);
                end 

            end 
        end 
    end 


    % Initialize a matrix to store the values for this NDF 
    counts_mat = nan(num_settings_levels, n_measures, min_num_readings, n_channels);

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
            counts_mat(ss, nn, 1:min_num_readings, :) = readings(1:min_num_readings, :);
        end

    end
    return ;

end
