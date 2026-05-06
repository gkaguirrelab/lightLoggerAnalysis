function barPlotSpdSlopes(best_fit_lines, activity_name)
% PLOT_SPD_SLOPES Creates a vertically stacked 4-tile bar plot with SEM error bars.
%
% Inputs:
%   - best_fit_lines: The data struct containing subject fields
%   - activity_name: Name of the activity to plot (e.g., 'walkIndoor' or 'walkOutdoor')
% Example
%{
    barPlotSpdSlopes(best_fit_lines, 'walkOutdoor')
%}
    subjects = fieldnames(best_fit_lines);
    num_subjects = length(subjects);
    
    conditions = {'a', 'c_lm', 'c_s'};
    condition_names = {'Achromatic', 'L-M', 'S'};
    num_conditions = length(conditions);

    % Define colors: {Dark, Light}
    % Row 1: Achromatic (Gray), Row 2: L-M (Red/Pink), Row 3: S (Blue)
    colors = {
        [0.3 0.3 0.3; 0.7 0.7 0.7], % Achromatic
        [0.7 0.2 0.2; 0.9 0.6 0.6], % L-M
        [0.0 0.3 0.7; 0.6 0.8 1.0]  % S
    };

    % --- Precomputation to find global maximum for consistent Y-axis limits ---
    all_max_vals = zeros(1, num_conditions);
    data_cell = cell(num_conditions, 1);
    
    for c = 1:num_conditions
        cond = conditions{c};
        fov_c_vals = zeros(1, num_subjects);
        fov_p_vals = zeros(1, num_subjects);
        non_c_vals = zeros(1, num_subjects);
        non_p_vals = zeros(1, num_subjects);

        for s = 1:num_subjects
            subj = subjects{s};
            if isfield(best_fit_lines.(subj), activity_name)
                data = best_fit_lines.(subj).(activity_name).(cond);
                
                % Invert the sign to plot alpha (-slope)
                fov_c_vals(s) = -data.virtuallyFoveated.center.slope;
                fov_p_vals(s) = -data.virtuallyFoveated.periphery.slope;
                non_c_vals(s) = -data.justProjection.center.slope;
                non_p_vals(s) = -data.justProjection.periphery.slope;
            end
        end

        m_fov_c = mean(fov_c_vals);
        m_fov_p = mean(fov_p_vals);
        m_non_c = mean(non_c_vals);
        m_non_p = mean(non_p_vals);

        se_fov_c = std(fov_c_vals) / sqrt(num_subjects);
        se_fov_p = std(fov_p_vals) / sqrt(num_subjects);
        se_non_c = std(non_c_vals) / sqrt(num_subjects);
        se_non_p = std(non_p_vals) / sqrt(num_subjects);

        fov_means = [m_fov_c, m_fov_p];
        fov_sem = [se_fov_c, se_fov_p];
        non_means = [m_non_c, m_non_p];
        non_sem = [se_non_c, se_non_p];

        all_means = [fov_means, non_means];
        all_max_vals(c) = max(all_means) + max([fov_sem, non_sem]) * 1.2;
        
        data_cell{c} = struct(...
            'fov_means', fov_means, ...
            'fov_sem', fov_sem, ...
            'non_means', non_means, ...
            'non_sem', non_sem ...
        );
    end
    
    global_max = max(all_max_vals);
    y_max = max(ceil(global_max), 1); % Round up to next whole number, at least 1

    figure('Color', 'w', 'Position', [100, 100, 700, 950]);
    t = tiledlayout(num_conditions + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    title(t, sprintf('Spd Slopes: %s', activity_name), 'FontSize', 14, 'FontWeight', 'bold');

    for c = 1:num_conditions
        cond_data = data_cell{c};

        % Create tile
        nexttile;
        hold on;

        % 1. Plot Non-Foveated Bars (x = 1, 2) with darker edges
        bar(1, cond_data.non_means(1), ...
            'FaceColor', colors{c}(1, :), ...
            'EdgeColor', colors{c}(1, :) * 0.7, ...
            'LineWidth', 1.5, ...
            'LineStyle', '--');
        
        bar(2, cond_data.non_means(2), ...
            'FaceColor', colors{c}(2, :), ...
            'EdgeColor', colors{c}(2, :) * 0.7, ...
            'LineWidth', 1.5, ...
            'LineStyle', '--');

        % 2. Plot Foveated Bars (x = 3.5, 4.5) with darker edges
        bar(3.5, cond_data.fov_means(1), ...
            'FaceColor', colors{c}(1, :), ...
            'EdgeColor', colors{c}(1, :) * 0.7, ...
            'LineWidth', 1.5, ...
            'LineStyle', '-');
        
        bar(4.5, cond_data.fov_means(2), ...
            'FaceColor', colors{c}(2, :), ...
            'EdgeColor', colors{c}(2, :) * 0.7, ...
            'LineWidth', 1.5, ...
            'LineStyle', '-');

        % 3. Add SEM Error Bars
        errorbar([1, 2], cond_data.non_means, cond_data.non_sem, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');
        errorbar([3.5, 4.5], cond_data.fov_means, cond_data.fov_sem, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

        % Format plot and axis styling
        set(gca, 'Layer', 'top', 'Box', 'off');
        xlim([0.5, 5.0]);
        xticks([]); % Turn off default ticks to make room for custom labels
        ylabel(condition_names{c});
        
        % Set uniform Y limits and whole number ticks (starting at 0)
        ylim([0, y_max]);
        yticks(0:y_max);
        
        % Remove the horizontal axis lines to avoid cluttered lower bounds
        set(gca, 'XColor', 'none');
        
        grid on;
        box off; % Turn off the bounding outline
    end

    % --- 4th Tile for Labels ---
    nexttile;
    hold on;
    xlim([0.5, 5.0]);
    ylim([-1, 1]);
    
    % Hide axes properties for the bottom tile so only text shows
    set(gca, 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
    set(gca, 'Color', 'none');
    
    % Item Labels 
    x_coords = [1, 2, 3.5, 4.5];
    itemNames = {'Center', 'Periphery', 'Center', 'Periphery'};
    
    for i = 1:4
        text(x_coords(i), 0.9, itemNames{i}, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 9);
    end

    % Group (Condition) Labels
    text(1.5, 0.5, 'Not foveated', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11);
    text(4.0, 0.5, 'Foveated', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 11);
        
    box off;
end