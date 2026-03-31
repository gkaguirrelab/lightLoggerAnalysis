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
%   Optionally, activities within each row can be reordered according to a
%   fixed preferred order:
%
%       indoor  -> work, chat, walkIndoor
%       outdoor -> sitBiopond, walkBiopond, walkOutdoor
%
%   Any remaining activities not listed above are appended to the end of
%   the row in alphabetical order.
%
% Inputs:
%   groupedActivityData - Struct or filepath. If a filepath is provided,
%       it should point to a .mat file containing a variable named
%       groupedActivityData.
%
% Optional key/value pairs:
%   exponent_clim           - false or 2-element numeric vector for exponent CLim.
%   variance_clim           - false or 2-element numeric vector for variance CLim.
%   spd_xlim                - false or 2-element numeric vector for SPD x-limits.
%   spd_ylim                - false or 2-element numeric vector for SPD y-limits.
%   title                   - Title stem used in the figure super-title and export.
%   output_dir              - Output directory for saving the grouped PDF.
%   overwrite_existing      - Logical. Whether to overwrite an existing export.
%   sort_by_preferred_order - Logical. If true, each row is reordered using
%                             the hard-coded preferred activity order for
%                             that group, with all other activities appended
%                             alphabetically afterward.
%
% Outputs:
%   groupedSpdHandle   - Figure handle for the grouped SPD figure.
%
% Example:
%   groupedFig = groupSPDs(groupedActivityData, ...
%       "title", "2001", ...
%       "spd_xlim", [0.5 60], ...
%       "spd_ylim", [1e-6 1e2], ...
%       "sort_by_preferred_order", true, ...
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
        options.sort_by_preferred_order = false
    end 

    % If groupedActivityData was passed in as a filepath, load it now
    if (isstring(groupedActivityData) || ischar(groupedActivityData))
        groupedActivityData = load(groupedActivityData).groupedActivityData; 
    end     

    disp(groupedActivityData)

    % The top-level field names correspond to the row groups
    group_names = fieldnames(groupedActivityData); 
    n_graph_rows = numel(group_names);

    % Determine the number of columns required by taking the largest group
    n_graph_cols = -inf; 
    for gg = 1:n_graph_rows
        group_name = group_names{gg};
        group_activities_struct = groupedActivityData.(group_name); 
        group_activities = fieldnames(group_activities_struct);

        % Optionally reorder activities within this row using the
        % group-specific preferred ordering
        if (options.sort_by_preferred_order)
            group_activities = iSortActivitiesByPreferredOrder(group_activities, group_name);
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

        % Optionally reorder activities within this row using the
        % group-specific preferred ordering
        if (options.sort_by_preferred_order)
            group_activities = iSortActivitiesByPreferredOrder(group_activities, group_name);
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


function sorted_activity_names = iSortActivitiesByPreferredOrder(activity_names, group_name)
% Reorder activity names according to a fixed preferred order for each
% group. Any activity names not explicitly listed are appended afterward in
% alphabetical order.
%
% Preferred order:
%   indoor  -> work, chat, walkIndoor
%   outdoor -> sitBiopond, walkBiopond, walkOutdoor

    % Start from alphabetically sorted names so "other" activities are
    % appended in a stable, predictable order
    activity_names = sort(activity_names);

    switch lower(group_name)
        case 'indoor'
            preferred_order = {'work', 'chat', 'walkIndoor'};
        case 'outdoor'
            preferred_order = {'sitBiopond', 'walkBiopond', 'walkOutdoor'};
        otherwise
            preferred_order = {};
    end

    sorted_activity_names = {};

    % First, add activities that appear in the preferred order list and
    % actually exist in this row
    for ii = 1:numel(preferred_order)
        preferred_name = preferred_order{ii};

        for jj = 1:numel(activity_names)
            if strcmp(activity_names{jj}, preferred_name)
                sorted_activity_names{end+1} = activity_names{jj}; %#ok<AGROW>
                break
            end
        end
    end

    % Next, append any remaining activities that were not explicitly listed
    for ii = 1:numel(activity_names)
        current_name = activity_names{ii};

        if ~ismember(current_name, sorted_activity_names)
            sorted_activity_names{end+1} = current_name; %#ok<AGROW>
        end
    end
end