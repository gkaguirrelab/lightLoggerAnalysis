function barPlotSpdSlopes(activity_name, options)
% PLOT_SPD_SLOPES Creates a vertically stacked 4-tile bar plot with SEM error bars and vibrant, semi-transparent data points.
%
% Inputs:
%   - activity_name: Name of the activity to plot (e.g., 'walkIndoor' or 'walkOutdoor')
%   - options: 
%       - src_dir: path to data
%       - color_mode: "a", "c_lm", "c_s" or an array of them
%       - show_beeswarm: true/false (default true)
% Example
%{
    barPlotSpdSlopes('chat', "color_mode", ["a"])
%}  
    arguments 
        activity_name 
        options.src_dir = fullfile(getpref("lightLogger", "dropboxBaseDir"), "FLIC_analysis", "lightLogger", "NEWscriptedIndoorOutdoorVideos2026") 
        options.color_mode = ["a", "c_lm", "c_s"]
        options.show_beeswarm = true
    end 
    
    % Load data
    spds_and_best_fits = loadSPDs(options.src_dir, "activities_to_process", {activity_name}, "color_modes_to_process", cellstr(options.color_mode)); 
    
    % Adapt structure to legacy plotting layout
    best_fit_lines = struct();
    color_modes = fieldnames(spds_and_best_fits);
    for cc = 1:numel(color_modes)
        color_mode = color_modes{cc};
        color_mode_struct = spds_and_best_fits.(color_mode);
        subjects_for_color = fieldnames(color_mode_struct);
        for ss = 1:numel(subjects_for_color)
            subject_id = subjects_for_color{ss};
            subject_struct = color_mode_struct.(subject_id);
            if (~isfield(subject_struct, activity_name)), continue; end
            activity_struct = subject_struct.(activity_name);
            projection_types = fieldnames(activity_struct);
            for pp = 1:numel(projection_types)
                projection_type = projection_types{pp};
                projection_struct = activity_struct.(projection_type);
                if (~isfield(projection_struct, "best_fit")), continue; end
                best_fit_struct = projection_struct.best_fit;
                if (isfield(best_fit_struct, "bestFit"))
                    best_fit_struct = best_fit_struct.bestFit;
                end
                best_fit_lines.(subject_id).(activity_name).(color_mode).(projection_type) = best_fit_struct;
            end
        end
    end
    
    subjects = fieldnames(best_fit_lines);
    num_subjects = length(subjects);
    conditions = color_modes;
    condition_names = cell(size(conditions));
    for cc = 1:numel(conditions)
        condition_names{cc} = iFormatConditionName(conditions{cc});
    end
    num_conditions = length(conditions);
    
    % Define subject colors
    subj_colors = turbo(num_subjects); 
    
    % Define Bar colors: {Darker, Lighter}
    colors_map = {
        [0.35 0.35 0.35; 0.75 0.75 0.75], % Achromatic
        [0.75 0.20 0.20; 0.90 0.55 0.55], % L-M
        [0.05 0.35 0.75; 0.55 0.75 1.00]  % S
    };
    
    % --- Data Extraction ---
    data_cell = cell(num_conditions, 1);
    for c = 1:num_conditions
        cond = conditions{c};
        raw = struct('fov_c', nan(1, num_subjects), 'fov_p', nan(1, num_subjects), ...
                     'non_c', nan(1, num_subjects), 'non_p', nan(1, num_subjects));
        
        for s = 1:num_subjects
            subj = subjects{s};
            if isfield(best_fit_lines, subj) && isfield(best_fit_lines.(subj).(activity_name), cond)
                data = best_fit_lines.(subj).(activity_name).(cond);
                raw.fov_c(s) = -data.virtuallyFoveated.center.slope;
                raw.fov_p(s) = -data.virtuallyFoveated.periphery.slope;
                raw.non_c(s) = -data.justProjection.center.slope;
                raw.non_p(s) = -data.justProjection.periphery.slope;
            end
        end
        data_cell{c}.raw = raw;
        data_cell{c}.means = [nanmean(raw.non_c), nanmean(raw.non_p), nanmean(raw.fov_c), nanmean(raw.fov_p)];
        data_cell{c}.sems = [nanstd(raw.non_c), nanstd(raw.non_p), nanstd(raw.fov_c), nanstd(raw.fov_p)] ./ sqrt(sum(~isnan(raw.non_c)));
    end
    
    % --- Plotting ---
    figure('Color', 'w', 'Position', [100, 100, 750, 950]);
    t = tiledlayout(num_conditions + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Spd Slopes: %s', activity_name), 'FontSize', 14, 'FontWeight', 'bold');
    x_pos = [1, 2, 3.5, 4.5];
    
    for c = 1:num_conditions
        cond_data = data_cell{c};
        nexttile; hold on;
        
        % Add dotted reference line at y = 2
        yline(2, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
        
        % 1. Plot Solid Bars
        for i = 1:4
            face_col = colors_map{c}(1 + mod(i-1,2), :);
            l_style = '-'; if i < 3, l_style = '--'; end
            
            bar(x_pos(i), cond_data.means(i), 0.8, ...
                'FaceColor', face_col, 'EdgeColor', face_col * 0.6, ...
                'LineWidth', 1.5, 'LineStyle', l_style, 'FaceAlpha', 1.0);
        end
        
        % 2. Semi-transparent jittered points
        if options.show_beeswarm
            rng(1); 
            jitter_width = 0.3;
            raw_fields = {'non_c', 'non_p', 'fov_c', 'fov_p'};
            
            for i = 1:4
                y_pts = cond_data.raw.(raw_fields{i});
                x_pts = x_pos(i) + (rand(size(y_pts)) - 0.5) * jitter_width;
                
                for s = 1:num_subjects
                    if isnan(y_pts(s)), continue; end
                    scatter(x_pts(s), y_pts(s), 55, 'filled', ...
                        'MarkerFaceColor', subj_colors(s,:), ...
                        'MarkerEdgeColor', [0.1 0.1 0.1], ...
                        'LineWidth', 0.5, ...
                        'MarkerFaceAlpha', 0.5);
                end
            end
        end
        
        % 3. Error Bars
        errorbar(x_pos, cond_data.means, cond_data.sems, 'k', 'LineWidth', 2.2, 'LineStyle', 'none');
        
        % Formatting
        set(gca, 'Layer', 'top', 'Box', 'off', 'XColor', 'none', 'TickDir', 'out');
        xlim([0.5, 5.0]);
        ylabel(['1/f exponent, ' condition_names{c}], 'FontSize', 13);
        ylim([0, 2.5]);
        yticks(0:0.5:2.5);
        grid on;
        ax = gca; ax.GridAlpha = 0.1;
    end
    
    % --- Labels Tile ---
    nexttile; hold on;
    xlim([0.5, 5.0]); ylim([-1, 1]);
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none');
    
    itemNames = {'Center', 'Periphery', 'Center', 'Periphery'};
    for i = 1:4
        text(x_pos(i), 0.9, itemNames{i}, 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    text(1.5, 0.4, 'Not foveated', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
    text(4.0, 0.4, 'Foveated', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
end

function condition_name = iFormatConditionName(condition)
    switch string(condition)
        case "a",   condition_name = 'Achromatic';
        case "c_lm", condition_name = 'L-M';
        case "c_s",  condition_name = 'S';
        otherwise,  condition_name = char(condition);
    end
end