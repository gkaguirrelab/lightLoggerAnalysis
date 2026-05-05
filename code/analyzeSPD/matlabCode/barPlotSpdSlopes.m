function barPlotSpdSlopes(best_fit_lines, activity_name)
% PLOT_SPD_SLOPES Creates a vertically stacked 3-tile bar plot with SEM error bars.
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

    figure('Color', 'w', 'Position', [100, 100, 700, 900]);
    t = tiledlayout(num_conditions, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Spd Slopes: %s', activity_name), 'FontSize', 14, 'FontWeight', 'bold');

    for c = 1:num_conditions
        cond = conditions{c};

        % Initialize arrays to hold subject values
        fov_c_vals = zeros(1, num_subjects);
        fov_p_vals = zeros(1, num_subjects);
        non_c_vals = zeros(1, num_subjects);
        non_p_vals = zeros(1, num_subjects);

        % Loop over subjects
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

        % Compute Mean and SEM
        m_fov_c = mean(fov_c_vals);
        m_fov_p = mean(fov_p_vals);
        m_non_c = mean(non_c_vals);
        m_non_p = mean(non_p_vals);

        se_fov_c = std(fov_c_vals) / sqrt(num_subjects);
        se_fov_p = std(fov_p_vals) / sqrt(num_subjects);
        se_non_c = std(non_c_vals) / sqrt(num_subjects);
        se_non_p = std(non_p_vals) / sqrt(num_subjects);

        % Create tile
        nexttile;
        hold on;

        fov_means = [m_fov_c, m_fov_p];
        fov_sem = [se_fov_c, se_fov_p];
        non_means = [m_non_c, m_non_p];
        non_sem = [se_non_c, se_non_p];

        % 1. Plot Foveated Bars (x = 1, 2) with solid edges
        b1 = bar([1, 2], fov_means, 'FaceColor', 'flat');
        b1.CData(1, :) = colors{c}(1, :); % Dark version (Center)
        b1.CData(2, :) = colors{c}(2, :); % Light version (Periphery)
        b1.LineStyle = '-';
        b1.LineWidth = 1.5;

        % 2. Plot Non-Foveated Bars (x = 4, 5) with dashed edges
        b2 = bar([4, 5], non_means, 'FaceColor', 'flat');
        b2.CData(1, :) = colors{c}(1, :); % Dark version (Center)
        b2.CData(2, :) = colors{c}(2, :); % Light version (Periphery)
        b2.LineStyle = '--';
        b2.LineWidth = 1.5;

        % 3. Add SEM Error Bars
        errorbar([1, 2], fov_means, fov_sem, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');
        errorbar([4, 5], non_means, non_sem, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

        % Format plot and axis styling
        set(gca, 'Layer', 'top', 'Box', 'off');
        xlim([0.5, 5.5]);
        xticks([]); % Turn off default ticks to make room for custom labels
        ylabel(condition_names{c});
        
        % Legend
        legend([b1, b2], {'Foveated', 'Not foveated'}, 'Location', 'Northwest');

        % Dynamically compute plot range limits
        all_means = [fov_means, non_means];
        max_val = max(all_means) + max([fov_sem, non_sem]) * 1.2;
        
        % Scale Y limits to give enough space at the bottom for labels
        ylim([-max_val * 0.45, max_val]);
        
        % Calculate placement positions for the text (within the new Y range)
        yContrast = -0.1 * max_val;
        yLight = -0.32 * max_val;

        % Item Labels (Centered under each bar)
        x_coords = [1, 2, 4, 5];
        itemNames = {'Center', 'Periphery', 'Center', 'Periphery'};
        
        for i = 1:4
            text(x_coords(i), yContrast, itemNames{i}, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 9, ...
                'Clipping', 'off');
        end

        % Group (Condition) Labels
        text(1.5, yLight, 'Foveated', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 11, 'Clipping', 'off');
        text(4.5, yLight, 'Not foveated', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 11, 'Clipping', 'off');
        
        grid on;
        box on;
    end
end