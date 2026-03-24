function groupedSpdHandle = groupSPDs(groupedActivityData, options)
% Plot grouped SPD-by-region summaries across multiple activities
%
% Syntax:
%   groupedSpdHandle = groupSPDs(groupedActivityData)
%   groupedSpdHandle = groupSPDs(groupedActivityData, options)
%
% Description:
%   This function groups multiple activity-level SPD summaries into a
%   single figure arranged as a tiled layout.
%
%   The input groupedActivityData is expected to be organized as:
%
%       groupedActivityData.(groupName).(activityName).virtuallyFoveated
%       groupedActivityData.(groupName).(activityName).justProjection
%
%   where each projection field contains either:
%       - a filepath to an activityData .mat file, or
%       - an already-loaded activityData struct
%
%   Each row of the output figure corresponds to one group, and each column
%   corresponds to one activity within that group.
%
%   Internally, this function calls plotSPDs for each activity and directs
%   the SPD-by-region plot to draw directly into the target tile axes. This
%   is more robust than generating temporary SPD figures and copying their
%   graphics objects into the grouped figure.
%
%   Optionally, activities within each row can be reordered so that sitting
%   activities are shown first and walking activities are shown second.
%
% Inputs:
%   groupedActivityData - Struct or filepath. If a filepath is provided,
%       it should point to a .mat file containing a variable named
%       groupedActivityData.
%
% Optional key/value pairs:
%   exponent_clim               - false or 2-element numeric vector for exponent CLim.
%   variance_clim               - false or 2-element numeric vector for variance CLim.
%   spd_xlim                    - false or 2-element numeric vector for SPD x-limits.
%   spd_ylim                    - false or 2-element numeric vector for SPD y-limits.
%   title                       - Title stem used in the figure super-title and export.
%   output_dir                  - Output directory for saving the grouped PDF.
%   overwrite_existing          - Logical. Whether to overwrite an existing export.
%   sort_sitting_before_walking - Logical. If true, each row is reordered so
%                                 activities with "sit" in the name appear first,
%                                 then activities with "walk" in the name, then all
%                                 remaining activities. Within each category, names
%                                 are alphabetically sorted.
%
% Outputs:
%   groupedSpdHandle   - Figure handle for the grouped SPD figure.
%
% Example:
%   groupedFig = groupSPDs(groupedActivityData, ...
%       "title", "2001", ...
%       "spd_xlim", [0.5 60], ...
%       "spd_ylim", [1e-6 1e2], ...
%       "sort_sitting_before_walking", true, ...
%       "output_dir", "/path/to/output");

    arguments 
        groupedActivityData
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.title = ""
        options.output_dir = ""
        options.overwrite_existing = false
        options.sort_sitting_before_walking = false
    end 

    % If groupedActivityData was passed in as a filepath, load it now
    if (isstring(groupedActivityData) || ischar(groupedActivityData))
        groupedActivityData = load(groupedActivityData).groupedActivityData; 
    end     

    % The top-level field names correspond to the row groups
    group_names = fieldnames(groupedActivityData); 
    n_graph_rows = numel(group_names);

    % Determine the number of columns required by taking the largest group
    n_graph_cols = -inf; 
    for gg = 1:n_graph_rows
        group_name = group_names{gg};
        group_activities_struct = groupedActivityData.(group_name); 
        group_activities = fieldnames(group_activities_struct);

        % Optionally reorder activities within this row
        if (options.sort_sitting_before_walking)
            group_activities = iSortActivitiesSitThenWalk(group_activities);
        else
            group_activities = sort(group_activities);
        end

        n_graph_cols = max([n_graph_cols, numel(group_activities)]);
    end 

    % Create a larger grouped figure so the tiles have more room and the
    % resulting PDF is less visually cramped
    groupedSpdHandle = figure('Units', 'pixels', 'Position', [100 100 2200 1300]);

    % Use looser spacing than before to give labels, legends, and titles
    % more breathing room
    tiledlayout(n_graph_rows, n_graph_cols, "TileSpacing", "loose", "Padding", "loose");

    % Populate each tile directly by asking plotSPDs to draw the SPD plot
    % into the destination axes for that tile
    for gg = 1:n_graph_rows
        group_name = group_names{gg};
        group_activities_struct = groupedActivityData.(group_name); 
        group_activities = fieldnames(group_activities_struct);

        % Optionally reorder activities within this row
        if (options.sort_sitting_before_walking)
            group_activities = iSortActivitiesSitThenWalk(group_activities);
        else
            group_activities = sort(group_activities);
        end

        for aa = 1:numel(group_activities)
            activity_name = group_activities{aa};
            activity_struct = group_activities_struct.(activity_name); 
            projection_types = fieldnames(activity_struct); 

            % Virtually foveated data is required
            virtuallyFoveatedActivityData = activity_struct.virtuallyFoveated; 
            assert(~strcmp(virtuallyFoveatedActivityData, ""));

            % justProjection data is optional, but if present, use it for
            % direct comparison in the SPD plot
            justProjectionActivityData = false; 
            if (ismember("justProjection", projection_types))
                justProjectionActivityData = activity_struct.justProjection; 
                assert(~strcmp(justProjectionActivityData, "")); 
            end 

            % Compute the tile index for this row / column
            tile_idx = (gg - 1) * n_graph_cols + aa;
            targetAxes = nexttile(tile_idx);
            hold(targetAxes, 'on');

            % Draw the SPD plot directly into the target axes. The exponent
            % and variance figures are still created internally by plotSPDs,
            % but we immediately close them since groupSPDs only needs the
            % SPD panel.
            [exponentMapHandle, varianceMapHandle, ~] = plotSPDs( ...
                virtuallyFoveatedActivityData, ...
                "exponent_clim", options.exponent_clim, ...
                "variance_clim", options.variance_clim, ...
                "spd_xlim", options.spd_xlim, ...
                "spd_ylim", options.spd_ylim, ...
                "justProjectionActivityData", justProjectionActivityData, ...
                "spd_target_axes", targetAxes ...
            );

            % Replace the default plotSPDs title with a compact per-tile
            % activity title better suited for grouped viewing
            title(targetAxes, activity_name, 'Interpreter', 'none', 'FontWeight', 'bold');

            % Keep axis labels consistent across grouped plots
            xlabel(targetAxes, 'Frequency [Hz]');
            ylabel(targetAxes, 'Power [contrast^2/Hz]');

            % On the first column of each row, prepend the group name to the
            % y-axis label so each row is clearly labeled by group
            if (aa == 1)
                ylabel(targetAxes, sprintf('\\bf%s\\rm\n%s', group_name, 'Power [contrast^2/Hz]'), 'Interpreter', 'tex');
            end

            % Close the temporary exponent / variance map figures created by
            % plotSPDs so they do not accumulate during grouped plotting
            if (~isempty(exponentMapHandle) && isgraphics(exponentMapHandle))
                close(exponentMapHandle);
            end
            if (~isempty(varianceMapHandle) && isgraphics(varianceMapHandle))
                close(varianceMapHandle);
            end
        end

        % If this row has fewer activities than the maximum number of
        % columns, fill the remaining tiles with blank axes
        for aa = (numel(group_activities) + 1):n_graph_cols
            tile_idx = (gg - 1) * n_graph_cols + aa;
            blankAxes = nexttile(tile_idx);
            axis(blankAxes, 'off');
        end
    end 

    % Add a bold super-title to the grouped figure
    sgtitle(sprintf('%s | SPDs Grouped By Subject', options.title), ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    % Export the grouped figure if requested
    if (~strcmp(options.output_dir, ""))
        assert(~strcmp(options.title, ""), ...
            'groupSPDs:MissingTitle', ...
            'options.title must be provided when saving output.');

        output_filepath = fullfile(options.output_dir, sprintf('%s_groupedSPDs.pdf', options.title));

        if (~exist(output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(groupedSpdHandle, output_filepath, 'ContentType', 'vector');
        end
    end
end


function sorted_activity_names = iSortActivitiesSitThenWalk(activity_names)
% Reorder activity names so sitting activities come first, walking
% activities come second, and all remaining activities come last.
% Within each category, names are sorted alphabetically.

    % Normalize to a row cell array of char/string-like names
    activity_names = sort(activity_names);

    sitting_names = {};
    walking_names = {};
    other_names = {};

    for ii = 1:numel(activity_names)
        current_name = activity_names{ii};
        lower_name = lower(current_name);

        if contains(lower_name, 'sit')
            sitting_names{end+1} = current_name; %#ok<AGROW>
        elseif contains(lower_name, 'walk')
            walking_names{end+1} = current_name; %#ok<AGROW>
        else
            other_names{end+1} = current_name; %#ok<AGROW>
        end
    end

    sorted_activity_names = [sort(sitting_names), sort(walking_names), sort(other_names)];
end