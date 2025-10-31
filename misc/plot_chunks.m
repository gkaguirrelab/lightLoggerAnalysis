function plot_chunks(chunks)
% Plot all of the sensors readings for all of the chunks in a parsed_chunks
% cell. Treat this as one large timeseries and video and label sensors 
% accordingly. Convert the counts to contrast units to put the sensors 
% on the same scale
%
% Syntax:
%   plot_chunks(chunks)
%
% Description:
%   Takes a cell of parsed chunks from a video recorded recording 
%   on the Raspberry Pi and plots the sensors together 
%   for the entire duration of the video in contrast units.
%
%
% Inputs:
%   chunks                - Cell. Cell array of parsed chunk 
%                           structs. 
%
% Outputs:
%
%   NONE                
%
% Examples:
%{
    path_to_experiment = '/Volumes/EXTERNAL1/test_folder_0';
    chunks = parse_chunks(path_to_experiment, true, true, true);
    plot_chunks(chunks, matlab_analysis_libraries_path, path_to_ms_util);
%}

    arguments 
        chunks; % Cell array of parsed chunks 
    end 
    
    % First, let's add the MATLAB helper libraries to the path 
    addpath(getpref('lightLoggerAnalysis', 'light_logger_libraries_matlab')); 

    % First, we will retrieve the Python MS util library to help 
    % us more easily plot it 
    ms_util = import_pyfile(getpref("lightLoggerAnalysis", "ms_util_path"));

    % Initialize containers for the flattened readings
    % of world and pupil
    world_chunks_t = [];
    world_chunks_v = [];
    world_chunks_v_contrast = [];

    pupil_chunks_t = [];
    pupil_chunks_v = [];
    pupil_chunks_v_contrast = [];

    % Initialize containers for the readings of 
    % MS sensors across time. 
    MS_light_sensing_chunks_t = [];
    MS_AS_chunks_v = [];
    MS_TS_chunks_v = [];
    MS_accelerometer_chunks_t = [];
    MS_chunkstarts = []; % Denote the start of each chunk for the MS
    MS_LS_chunks_v = [];
    MS_LS_chunks_v_std = [];
    MS_sunglasses_chunks_v = [];
    MS_LS_window_length = double(ms_util.MS_LS_BUFFER_SIZE); % This defines the sliding window used to
                                                             % plot the std of accelerometer readings.
                                                             % This is default set to 1 packet (1 second)

    % Initialize containers for the approximate FPS values 
    world_chunks_fps = [];
    pupil_chunks_fps = [];
    MS_chunks_fps = []; 

    % Initialize containers to track the world and pupil sensors' performance 
    world_frame_means = [];
    pupil_frame_means = [];
    
    world_chunks_agc_settings = [];
    pupil_chunks_agc_settings = []; 

    % Initialize an array to track start and 
    % end times for all of the chunks of the sensors. 
    % We will use this to find the points at which 
    % we have all sensor readings
    world_start_ends = {};
    pupil_start_ends = {}; 

    % Iterate over the chunks
    for cc = 1:numel(chunks)
        fprintf('Plotting chunk: %d\n', cc);

        % Retrieve the sensor readings for this chunk
        world_chunk_t = chunks{cc}.W.t;
        world_chunk_v = chunks{cc}.W.v;

        % Convert raw frames to array of mean pixels if we have been passed mean frames 
        if(numel(size(world_chunk_v)) > 2)
            world_chunk_v = mean(world_chunk_v, [2, 3])'; 
        end 

        world_chunk_settings = chunks{cc}.W.settings; 

        pupil_chunk_t = chunks{cc}.P.t;
        pupil_chunk_v = chunks{cc}.P.v;
        % Convert raw frames to array of mean pixels if we have been passed mean frames 
        if(numel(size(pupil_chunk_v)) > 2)
            pupil_chunk_v = mean(pupil_chunk_v, [2, 3])'; 
        end 

        pupil_chunk_settings = chunks{cc}.P.settings; 

        MS_light_sensing_chunk_t = chunks{cc}.M.t.AS;
        MS_accelerometer_chunk_t = chunks{cc}.M.t.LS;
        MS_AS_chunk_v = chunks{cc}.M.v.AS; 
        MS_TS_chunk_v = chunks{cc}.M.v.TS; 
        MS_LS_chunk_v = chunks{cc}.M.v.LS; 

        MS_sunglasses_chunk_v = chunks{cc}.S.v; 

        % Only interact with those sensors for which we have readings 
        % and follow the following proccedure: 
        %   1. Retrieve the start and end times for the current chunk 
        %      for that sensor
        %   2. Convert all the sensors to contrast units 
        %   3. Calculate the average FPS 
        %   4. Append them to the flattened arrays
        %   Note: this procedure differs for the MS as there are two separate t vecotrs, 
        %         one for the instaneous readings of the light sensors as well as one for 
        %         the buffer of accelerometer values. In addition, we do not convert these 
        %         to contrast
        if(numel(world_chunk_t) > 0) 
            % Retrieve the start and end times for the current chunk 
            world_start_ends{cc} = [world_chunk_t(1); world_chunk_t(end)];

            % Convert all the readings to contrast units 
            world_chunk_v_contrast = (world_chunk_v - mean(world_chunk_v)) / mean(world_chunk_v);
            
            % Collect the AGC settigns
            world_chunk_agc_settings = [world_chunk_settings.Again', world_chunk_settings.Dgain', world_chunk_settings.exposure']; 

            % Calculate the average FPS 
            world_chunk_fps = numel(world_chunk_t) / (world_chunk_t(end) - world_chunk_t(1)) ;

            % Append them to the flattened arrays
            world_chunks_t = [world_chunks_t, world_chunk_t];
            world_chunks_v = [world_chunks_v, world_chunk_v]; 
            world_chunks_v_contrast = [world_chunks_v_contrast, world_chunk_v_contrast]; 
            world_chunks_fps = [world_chunks_fps, world_chunk_fps]; 
            world_chunks_agc_settings = [world_chunks_agc_settings ; world_chunk_agc_settings]; 
            world_frame_means = [world_frame_means, world_chunk_v];

        end 

        if(numel(pupil_chunk_t) > 0) 
            % Retrieve the start and end times for the current chunk 
            pupil_start_ends{cc} = [pupil_chunk_t(1); pupil_chunk_t(end)];
            
            % Convert all the readings to contrast units 
            pupil_chunk_v_contrast = (pupil_chunk_v - mean(pupil_chunk_v)) / mean(pupil_chunk_v);

            % Collect the AGC settigns
            pupil_chunk_agc_settings = [pupil_chunk_settings.Again', pupil_chunk_settings.Dgain', pupil_chunk_settings.exposure']; 

            % Calculate the average FPS 
            pupil_chunk_fps = numel(pupil_chunk_t) / (pupil_chunk_t(end) - pupil_chunk_t(1));

            % Append them to the flattened arrays
            pupil_chunks_t = [pupil_chunks_t, pupil_chunk_t];
            pupil_chunks_v = [pupil_chunks_v, pupil_chunk_v]; 
            pupil_chunks_v_contrast = [pupil_chunks_v_contrast, pupil_chunk_v_contrast]; 
            pupil_chunks_fps = [pupil_chunks_fps, pupil_chunk_fps]; 
            pupil_chunks_agc_settings = [pupil_chunks_agc_settings ; pupil_chunk_agc_settings]; 
            pupil_frame_means = [pupil_frame_means, pupil_chunk_v];
        end 

        if(numel(MS_light_sensing_chunk_t) > 0)
            % Calculate the average FPS 
            MS_chunk_fps = numel(MS_light_sensing_chunk_t) / (MS_light_sensing_chunk_t(end) - MS_light_sensing_chunk_t(1));

            % Append the chunk readings to the flattened arrays
            MS_chunks_fps = [MS_chunks_fps, MS_chunk_fps]; 
            MS_light_sensing_chunks_t = [MS_light_sensing_chunks_t, MS_light_sensing_chunk_t];
            MS_AS_chunks_v = [MS_AS_chunks_v; MS_AS_chunk_v];
            MS_TS_chunks_v = [MS_TS_chunks_v; MS_TS_chunk_v];
        
            MS_accelerometer_chunks_t = [MS_accelerometer_chunks_t, MS_accelerometer_chunk_t];
            MS_LS_chunks_v = [MS_LS_chunks_v; MS_LS_chunk_v];
    
            MS_sunglasses_chunks_v = [MS_sunglasses_chunks_v, MS_sunglasses_chunk_v];

            MS_chunkstarts = [MS_chunkstarts, MS_light_sensing_chunk_t(1)];


        end 

    end

    % Create copies of the arrays so we can display the waves only where 
    % all sensor readings are present 
    present_world_chunks_v_contrast = ones(size(world_chunks_v_contrast));
    present_pupil_chunks_v_contrast = ones(size(pupil_chunks_v_contrast)); 
    present_MS_AS_chunks_v = ones(size(MS_AS_chunks_v));

    % Iterate over the starts and ends of each of 
    % the chunks and set the gap periods to nan in the plot
    for gg = 2:numel(world_start_ends)
        previous_start_end = world_start_ends{gg-1};
        current_start_end = world_start_ends{gg};

        previous_end = previous_start_end(2);
        current_start = current_start_end(1);

        % Find the indices for each sensor that are present during the world gap
        world_missing_indices = find(previous_end < world_chunks_t & world_chunks_t < current_start);
        pupil_missing_indices = find(previous_end < pupil_chunks_t & pupil_chunks_t < current_start);
        MS_AS_missing_indices = find(previous_end < MS_light_sensing_chunks_t & MS_light_sensing_chunks_t < current_start);

        % Set the values during a sensor gap to NAN
        present_world_chunks_v_contrast(world_missing_indices) = 0; 
        present_pupil_chunks_v_contrast(pupil_missing_indices) = 0; 
        present_MS_AS_chunks_v(MS_AS_missing_indices) = 0; 

    end

    for gg = 2:numel(pupil_start_ends)
        previous_start_end = pupil_start_ends{gg-1};
        current_start_end = pupil_start_ends{gg};

        previous_end = previous_start_end(2);
        current_start = current_start_end(1);

        % Find the indices for each sensor that are present during the world gap
        world_missing_indices = find(previous_end < world_chunks_t & world_chunks_t < current_start);
        pupil_missing_indices = find(previous_end < pupil_chunks_t & pupil_chunks_t < current_start);
        MS_AS_missing_indices = find(previous_end < MS_light_sensing_chunks_t & MS_light_sensing_chunks_t < current_start);

        % Set the values during a sensor gap to NAN
        present_world_chunks_v_contrast(world_missing_indices) = 0; 
        present_pupil_chunks_v_contrast(pupil_missing_indices) = 0; 
        present_MS_AS_chunks_v(MS_AS_missing_indices) = 0; 

    end 

    % Replace the RAW (hard to interpret) accelerometer readings with 
    % the std of the last N seconds for each point 
    MS_LS_chunks_v_std = zeros(size(MS_LS_chunks_v));
    for rr = 1:size(MS_LS_chunks_v, 1) % Iterate over the rows of the accelerometer 
        % For the first N readings, ignore since we do not have enough readings
        % for the context window 
        if(rr <= MS_LS_window_length)
            continue;
        end 

        % Otherwise, we calculate the standard deviation of the last N readings 
        % for each channel via the context window 
        context_window = MS_LS_chunks_v(rr-MS_LS_window_length:rr, :); 
        context_window_std = std(context_window);    
        MS_LS_chunks_v_std(rr, :) = context_window_std; 

    end
    
    % Find the min/max time (x-axis) of the sensors we will plot together
    min_common_time = min([min(world_chunks_t), min(pupil_chunks_t), min(MS_light_sensing_chunks_t), min(MS_accelerometer_chunks_t)]); 
    max_common_time = max([max(world_chunks_t), max(pupil_chunks_t), max(MS_light_sensing_chunks_t), max(MS_accelerometer_chunks_t)]); 
    
    %%%%%%%%%%%{ Plot Summary Card with important information at a glance }%%%%%%%%%%%
    plot_summary_card(world_chunks_t, world_chunks_v, world_chunks_v_contrast, present_world_chunks_v_contrast, world_chunks_fps,...
                      pupil_chunks_t, pupil_chunks_v, pupil_chunks_v_contrast, present_pupil_chunks_v_contrast, pupil_chunks_fps,... 
                      MS_light_sensing_chunks_t, MS_AS_chunks_v, present_MS_AS_chunks_v, MS_TS_chunks_v, MS_chunks_fps,...
                      MS_accelerometer_chunks_t, MS_LS_chunks_v, MS_LS_chunks_v_std,... 
                      MS_chunkstarts,...  
                      MS_sunglasses_chunks_v,...
                      min_common_time, max_common_time...
                     )

    % Generate per sensor figures with expanded information

    %%%%%%%%%%%{ Plot World Card }%%%%%%%%%%%
    plot_world_card(world_chunks_t, world_frame_means, world_chunks_agc_settings,...
                    min_common_time, max_common_time...
                   );

    %%%%%%%%%%%{ Plot Pupil Card }%%%%%%%%%%%
    plot_pupil_card(pupil_chunks_t, pupil_frame_means, pupil_chunks_agc_settings,...
                    min_common_time, max_common_time...
                   );

    %%%%%%%%%%%{ Plot MS Card }%%%%%%%%%%%
    plot_ms_and_sunglasses_card(MS_light_sensing_chunks_t, MS_AS_chunks_v, MS_TS_chunks_v,...
                                MS_accelerometer_chunks_t, MS_LS_chunks_v, MS_LS_chunks_v_std,...
                                MS_chunkstarts,...
                                MS_sunglasses_chunks_v,...
                                min_common_time, max_common_time... 
                               ); 


end

% Local function to plot the summarized figure 
function plot_summary_card(world_chunks_t, world_chunks_v, world_chunks_v_contrast, present_world_chunks_v_contrast, world_chunks_fps,...
                           pupil_chunks_t, pupil_chunks_v, pupil_chunks_v_contrast, present_pupil_chunks_v_contrast, pupil_chunks_fps,... 
                           MS_light_sensing_chunks_t, MS_AS_chunks_v, present_MS_AS_chunks_v, MS_TS_chunks_v, MS_chunks_fps,...
                           MS_accelerometer_chunks_t, MS_LS_chunks_v, MS_LS_chunks_v_std,... 
                           MS_chunkstarts,...  
                           MS_sunglasses_chunks_v,...
                           min_common_time, max_common_time...
                          )
     % Open a figure for plotting 
    % information about the performance of sensors combined 
    figure ; 
    t = tiledlayout(3,3); 
    title(t, "Sensors' Performance Summarized", "FontWeight", "bold");

    % First, plot the world sensor counts 
    nexttile; 
    title('World Sensor Counts');
    hold on ; 
    xlabel('Time [s]');
    xlim([min_common_time, max_common_time]); 

    % Plot the contrast on the left side 
    yyaxis left; 
    ylabel('Contrast');
    hold on ; 

    plot(world_chunks_t, world_chunks_v_contrast, 'x', 'DisplayName', 'WorldContrast');
    ylim([-0.5, 1]); 

    % Plot the raw mean sensor counts on the right side 
    yyaxis right; 
    ylabel('Raw Mean');
    hold on ; 

    plot(world_chunks_t, world_chunks_v, 'x', 'DisplayName', 'WorldMeanCount');
    ylim([-5, 260]); 

    
    % Next, plot the pupil sensor counts 
    nexttile; 
    title("Pupil Sensor Counts"); 
    hold on ; 
    xlabel('Time [s]');
    xlim([min_common_time, max_common_time]); 

    % Plot the contrast on the left side  
    yyaxis left; 
    ylabel('Contrast'); 
    hold on; 
    
    plot(pupil_chunks_t, pupil_chunks_v_contrast, 'x', 'DisplayName', 'Pupil');
    ylim([-0.5, 1]);

    yyaxis right; 
    ylabel('Raw Mean');
    hold on; 

    plot(pupil_chunks_t, pupil_chunks_v, 'x', 'DisplayName', 'PupilMeanCount');
    ylim([-5, 260]); 

    
    % Then plot the MS sensor counts. 
    % Note: Need to do dual x axis plot here 
    %       to put both light sensing chips 
    %       on the same plot
    nexttile; 
    title("MS Light Sensors Counts"); 
    hold on; 
    xlabel('Time [s]');
    
    yyaxis left; 

    % First, we plot the AS channels
    for cc = 1:size(MS_AS_chunks_v, 2)
        plot(MS_light_sensing_chunks_t, MS_AS_chunks_v(:, cc), 'x', "DisplayName", sprintf("ch%d", cc));
        hold on; 
    end 
    ylabel("AS Counts"); 
    ylim([-5, 67000]); % Define an unsigned 16 bit range 
    legend show; 
    
    yyaxis right; 
    for cc = 1:size(MS_TS_chunks_v, 2)
        plot(MS_light_sensing_chunks_t, MS_TS_chunks_v(:, cc), 'x', "DisplayName", sprintf("ch%d", cc-1));
    end 
    ylabel("TS Counts");
    ylim([-5, 67000]) % Define an unsigned 16 bit range 
    xlim([min_common_time, max_common_time]); 
    legend show; 

    % Plot the RAW accelerometer readings
    nexttile; 
    title("MS Accelerometer Readings [Raw]"); 
    hold on ; 
    % LS accelerometer has meaningful channel label names, 
    % not just numbers 
    channel_labels = {'X', 'Y', 'Z', 'ΩP', 'ΩR', 'ΩY'}; 
    for cc = 1:size(MS_LS_chunks_v, 2)
        plot(MS_accelerometer_chunks_t, MS_LS_chunks_v(:, cc), 'x', "DisplayName", channel_labels{cc});
    end 

    % Now, let's add vertical lines to denote chunk start points 
    if(numel(MS_chunkstarts) > 0)
        xline(MS_chunkstarts, '--r', 'Chunk Start', 'HandleVisibility', 'Off')
    end 

    legend show ; 
    xlabel("Time [s]");
    ylabel("Count");
    ylim([-33000, 33000]); % Set the ylim to show between 16 bit signed int range smoothly 
    xlim([min_common_time, max_common_time]); 

    % Plot the window calcualted STD accelerometer readings
    nexttile; 
    title("MS Accelerometer Readings [Window std]"); 
    hold on ; 
    % LS accelerometer has meaningful channel label names, 
    % not just numbers 
    channel_labels = {'X', 'Y', 'Z', 'ΩP', 'ΩR', 'ΩY'}; 
    for cc = 1:size(MS_LS_chunks_v, 2)
        plot(MS_accelerometer_chunks_t, MS_LS_chunks_v_std(:, cc), 'x', "DisplayName", channel_labels{cc});
    end 
    
    % Now, let's add vertical lines to denote chunk start points 
    if(numel(MS_chunkstarts) > 0)
        xline(MS_chunkstarts, '--r', 'Chunk Start', 'HandleVisibility', 'Off');
    end 

    legend show ; 
    xlabel("Time [s]");
    ylabel("In-Window std");
    ylim([-33000, 33000]); % Set the ylim to show between 16 bit signed int range smoothly 
    xlim([min_common_time, max_common_time]); 

    % Plot the sunglasses readings
    nexttile; 
    title("Sunglasses (Hall Effect Sensor) Readings"); 
    hold on; 
    plot(MS_light_sensing_chunks_t, MS_sunglasses_chunks_v, 'x', "DisplayName", "Hall Effect Sensor"); 
    legend show; 
    xlabel("Time [s]");
    ylabel("Count");
    ylim([-5, (2^12)-1 + 100]); % Set the ylim to show between 12 bit unsigned int range smoothly 
    xlim([min_common_time, max_common_time]); 

    % Create another figure to just show the waves where all sensor values are present overlaid 
    nexttile;
    title('Sensor Reading Overlap');
    hold on ; 
    xlabel('Time [s]');
    ylabel('All Sensors Present');
    xlim([min_common_time, max_common_time]); 

    % Add the values to the plot
    plot(world_chunks_t, present_world_chunks_v_contrast, 'kx', 'DisplayName', 'World');
    % Put a small constant offset here so we can see the gaps between readings better. 
    plot(pupil_chunks_t, present_pupil_chunks_v_contrast, 'rx', 'DisplayName', 'Pupil');

    if(size(present_MS_AS_chunks_v, 2) > 1) % If we have a proper reading, that is, we have all the channels 
        plot(MS_light_sensing_chunks_t, present_MS_AS_chunks_v(:, 1), 'bx', 'DisplayName', "AS1");
    end
    ylim([-0.25, 1.25]);
    legend show; 


    % Generate a figure to display the average FPS per sensor 
    nexttile;   
    title("Avg FPS by Sensor");
    hold on;
    xlabel("Chunk #");
    ylabel("FPS"); 
    ylim([-5, 205]); 

    plot(world_chunks_fps, "DisplayName", "World");
    plot(pupil_chunks_fps, "DisplayName", "Pupil");
    plot(MS_chunks_fps, "DisplayName", "MS");

    legend show; 

end 


% Local function to plot the world camera's figure 
function plot_world_card(world_chunks_t, world_frame_means, world_chunks_agc_settings,...
                         min_common_time, max_common_time...
                        )
    % Generate a figure to display the World Camera mean and settings over time 
    figure ; 
    t = tiledlayout(1, 2); 
    title(t, "World Camera", "FontWeight", "bold");

    % Access the first subplot, the leftmost one
    % to display mean frame count over time
    nexttile; 
    title("Mean Pixel Intensity by Frame"); 
    hold on ; 
    plot(world_chunks_t, world_frame_means, "-x", "DisplayName", "World"); 
    xlim([min_common_time, max_common_time]);
    ylim([-5, 260]);
    ylabel("Mean Pixel Intensity"); 
    xlabel("Time [s]");
    legend show; 
    
    
    % Access the second subplot, the rightmost one, 
    % to display how the settings control this behavior 
    nexttile;
    title("AGC Settings");
    hold on;
    xlabel("Time [s]"); 
    xlim([min_common_time, max_common_time]);

    yyaxis left; 
    
    if(size(world_chunks_agc_settings, 1) > 0)
        plot(world_chunks_t, world_chunks_agc_settings(:, 1), "bx", "DisplayName", "AnalogueGain");
    end
    hold on; 
    if(size(world_chunks_agc_settings, 1) > 0)
        plot(world_chunks_t, world_chunks_agc_settings(:, 2), "gx", "DisplayName", "DigitalGain");
    end
    ylabel("Gain"); 
    ylim([-1, 11.5]); 

    yyaxis right; 
    if(size(world_chunks_agc_settings, 1) > 0)
        plot(world_chunks_t, world_chunks_agc_settings(:, 3), "o", "DisplayName", "Exposure [s]");
    end
    hold on; 
    ylabel("Exposure Time"); 
    ylim([0, 10000]); 

    legend show;  

end 

% Local function to plot the pupil camera's figure 
function plot_pupil_card(pupil_chunks_t, pupil_frame_means, pupil_chunks_agc_settings,...
                         min_common_time, max_common_time...
                        )
    % Generate a figure to display the Pupil Camera mean and settings over time 
    figure ; 
    t = tiledlayout(1, 2); 
    title(t, "Pupil Camera", "FontWeight", "bold");

    % Access the first subplot, the leftmost one
    % to display mean frame count over time
    nexttile; 
    title("Mean Pixel Intensity"); 
    hold on ; 
    plot(pupil_chunks_t, pupil_frame_means, "-x", "DisplayName", "Pupil"); 
    xlim([min_common_time, max_common_time]); 
    ylim([-5, 260]);
    ylabel("Mean Pixel Intensity"); 
    xlabel("Time [s]");
    legend show; 

        % Access the second subplot, the rightmost one, 
    % to display how the settings control this behavior 
    nexttile;
    title("AGC Settings");
    hold on;
    xlim([min_common_time, max_common_time]);

    yyaxis left; 
    
    if(size(pupil_chunks_agc_settings, 1) > 0)
        plot(pupil_chunks_t, pupil_agc_settings(:, 1), "bx", "DisplayName", "AnalogueGain");
    end
    hold on; 
    if(size(pupil_chunks_agc_settings, 1) > 0)
        plot(pupil_chunks_t, pupil_agc_settings(:, 2), "gs", "DisplayName", "DigitalGain");
    end
    ylabel("Gain"); 
    ylim([-1, 11]); 

    yyaxis right; 
    if(size(pupil_chunks_agc_settings, 1) > 0)
        plot(pupil_chunks_t, pupil_agc_settings(:, 3), "o", "DisplayName", "Exposure [units]");
    end
    hold on; 
    ylabel("Exposure [units]"); 
    ylim([-2, 515]); 

    legend show; 

end 

% Local function to plot the MS's figure 
function plot_ms_and_sunglasses_card(MS_light_sensing_chunks_t, MS_AS_chunks_v, MS_TS_chunks_v,...
                                     MS_accelerometer_chunks_t, MS_LS_chunks_v, MS_LS_chunks_v_std,...
                                     MS_chunkstarts,...
                                     MS_sunglasses_chunks_v,...
                                     min_common_time, max_common_time...
                                    )
    % Open a figure for plotting 
    % information about the MS 
    figure ; 
    t = tiledlayout(3,3); 
    title(t, "MS + Sunglasses", "FontWeight", "bold"); 

    % Plot all of the channels of the AS chip 

    %%%%%%%%%%%{ Plot AS Chip }%%%%%%%%%%%
    nexttile; 
    title("AS Channel Readings"); 
    hold on; 
    for cc = 1:size(MS_AS_chunks_v, 2)
        plot(MS_light_sensing_chunks_t, log10(MS_AS_chunks_v(:, cc)), 'x', "DisplayName", sprintf("ch%d", cc));
    end 
    xlabel("Time [s]");
    ylabel("Count [log]");
    legend show ; 

    %%%%%%%%%%%{ Plot TS Chip }%%%%%%%%%%%
    nexttile;
    title("TS Channel Readings"); 
    hold on ; 
    for cc = 1:size(MS_TS_chunks_v, 2)
        plot(MS_light_sensing_chunks_t, log10(MS_TS_chunks_v(:, cc)), 'x', "DisplayName", sprintf("ch%d", cc-1));
    end 
    xlabel("Time [s]");
    ylabel("Count [log]");
    legend show ; 


    %%%%%%%%%%%{ Plot LS Chip [RAW] }%%%%%%%%%%%
    nexttile; 
    title("LS Readings [Raw]"); 
    hold on ; 
    % LS accelerometer has meaningful channel label names, 
    % not just numbers 
    channel_labels = {'X', 'Y', 'Z', 'ΩP', 'ΩR', 'ΩY'}; 
    for cc = 1:size(MS_LS_chunks_v, 2)
        plot(MS_accelerometer_chunks_t, MS_LS_chunks_v(:, cc), 'x', "DisplayName", channel_labels{cc});
    end 

    % Now, let's add vertical lines to denote chunk start points 
    if(numel(MS_chunkstarts) > 0)
        xline(MS_chunkstarts, '--r', 'Chunk Start', 'HandleVisibility', 'Off')
    end 

    legend show ; 
    xlabel("Time [s]");
    ylabel("Count");
    ylim([-33000, 33000]); % Set the ylim to show between 16 bit signed int range smoothly 
    xlim([min_common_time, max_common_time]); 

    %%%%%%%%%%%{ Plot LS Chip [WINDOW STD]}%%%%%%%%%%%    
    nexttile; 
    title("LS Readings [Window std]"); 
    hold on ; 
    % LS accelerometer has meaningful channel label names, 
    % not just numbers 
    channel_labels = {'X', 'Y', 'Z', 'ΩP', 'ΩR', 'ΩY'}; 
    for cc = 1:size(MS_LS_chunks_v, 2)
        plot(MS_accelerometer_chunks_t, MS_LS_chunks_v_std(:, cc), 'x', "DisplayName", channel_labels{cc});
    end 
    
    % Now, let's add vertical lines to denote chunk start points 
    if(numel(MS_chunkstarts) > 0)
        xline(MS_chunkstarts, '--r', 'Chunk Start', 'HandleVisibility', 'Off');
    end

    legend show ; 
    xlabel("Time [s]");
    ylabel("In-Window std");
    ylim([-33000, 33000]); % Set the ylim to show between 16 bit signed int range smoothly 
    xlim([min_common_time, max_common_time]); 

    %%%%%%%%%%%{ Plot Sunglasses }%%%%%%%%%%%
    nexttile; 
    title("Sunglasses Readings"); 
    hold on; 
    plot(MS_light_sensing_chunks_t, MS_sunglasses_chunks_v, "DisplayName", "Hall Effect Sensor"); 
    legend show; 
    xlabel("Time [s]");
    ylabel("Count");
    % Set the ylim to show between 12 bit unsigned int range smoothly 
    ylim([-5, (2^12)-1 + 100]); 



end 
