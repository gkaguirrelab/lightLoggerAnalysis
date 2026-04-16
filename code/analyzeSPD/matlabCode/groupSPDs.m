function [groupedExponentHandle, groupedVarianceHandle, groupedSpdHandle] = groupSPDs(groupedActivityData, options)
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
        options.n_participants = 1; 
        options.groups_across_participant_std = false; 
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

        % Optionally reorder activities within this row using the
        % group-specific preferred ordering
        if (options.sort_by_preferred_order)
            group_activities = iSortActivitiesByPreferredOrder(group_activities, group_name);
        else
            group_activities = sort(group_activities);
        end

        n_graph_cols = max([n_graph_cols, numel(group_activities)]);
    end 

    % Find the STD by group 
    if(options.n_participants > 1)
        groups_across_participant_std = find_across_participant_std(groupedActivityData, options.n_participants);
    end 

    % Create grouped exponent-map figure and KEEP the tiledlayout handle
    groupedExponentHandle = figure('Units', 'pixels', 'Position', [100 100 2200 1300]);
    groupedExponentLayout = tiledlayout(groupedExponentHandle, n_graph_rows, n_graph_cols, ...
        "TileSpacing", "loose", "Padding", "loose");

    % Create grouped variance-map figure and KEEP the tiledlayout handle
    groupedVarianceHandle = figure('Units', 'pixels', 'Position', [120 120 2200 1300]);
    groupedVarianceLayout = tiledlayout(groupedVarianceHandle, n_graph_rows, n_graph_cols, ...
        "TileSpacing", "loose", "Padding", "loose");

    % Create grouped SPD figure and KEEP the tiledlayout handle
    groupedSpdHandle = figure('Units', 'pixels', 'Position', [140 140 2200 1300]);
    groupedSpdLayout = tiledlayout(groupedSpdHandle, n_graph_rows, n_graph_cols, ...
        "TileSpacing", "loose", "Padding", "loose");

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
            % direct comparison in the SPD plot and map plots
            justProjectionActivityData = false; 
            if (ismember("justProjection", projection_types))
                justProjectionActivityData = activity_struct.justProjection; 
                assert(~strcmp(justProjectionActivityData, "")); 
            end 

            % Compute the tile index for this row / column
            tile_idx = (gg - 1) * n_graph_cols + aa;

            % -----------------------------------------------------------------
            % SPD tile
            % -----------------------------------------------------------------
            spdTargetAxes = nexttile(groupedSpdLayout, tile_idx);
            hold(spdTargetAxes, 'on');

            % Retrieve the STD for this group 
            if(options.n_participants > 1)
                group_std = groups_across_participant_std(char(group_name)); 
            else
                group_std = false; 
            end 

            [exponentMapHandle, varianceMapHandle, ~] = plotSPDs( ...
                virtuallyFoveatedActivityData, ...
                "exponent_clim", options.exponent_clim, ...
                "variance_clim", options.variance_clim, ...
                "spd_xlim", options.spd_xlim, ...
                "spd_ylim", options.spd_ylim, ...
                "justProjectionActivityData", justProjectionActivityData, ...
                "spd_target_axes", spdTargetAxes, ...
                "acrossParticipantRegionSTD", group_std,... 
                "num_participants", options.n_participants... 
            );

            title(spdTargetAxes, activity_name, 'Interpreter', 'none', 'FontWeight', 'bold');
            xlabel(spdTargetAxes, 'Frequency [Hz]');

            if (aa == 1)
                ylabel(spdTargetAxes, sprintf('\\bf%s\\rm\n%s', group_name, 'Power [contrast^2/Hz]'), 'Interpreter', 'tex');
            else
                ylabel(spdTargetAxes, 'Power [contrast^2/Hz]');
            end

             % -----------------------------------------------------------------
            % Exponent tile
            % Recreate the full plotSPDs map layout directly in the grouped
            % figure using ordinary axes/colorbars so exportgraphics works.
            % -----------------------------------------------------------------
            exponentTileAxes = nexttile(groupedExponentLayout, tile_idx);
            axis(exponentTileAxes, 'off');
            drawnow();

            exponentTilePosition = exponentTileAxes.Position;
            delete(exponentTileAxes);

            iCopyEntireMapFigureFromPlotSPDs( ...
                exponentMapHandle, ...
                groupedExponentHandle, ...
                exponentTilePosition, ...
                activity_name, ...
                group_name, ...
                aa == 1);

            % -----------------------------------------------------------------
            % Variance tile
            % Same idea as exponent tile.
            % -----------------------------------------------------------------
            varianceTileAxes = nexttile(groupedVarianceLayout, tile_idx);
            axis(varianceTileAxes, 'off');
            drawnow();

            varianceTilePosition = varianceTileAxes.Position;
            delete(varianceTileAxes);

            iCopyEntireMapFigureFromPlotSPDs( ...
                varianceMapHandle, ...
                groupedVarianceHandle, ...
                varianceTilePosition, ...
                activity_name, ...
                group_name, ...
                aa == 1);
        end 

        % If this row has fewer activities than the maximum number of
        % columns, fill the remaining tiles with blank axes in all 3 figures
        for aa = (numel(group_activities) + 1):n_graph_cols
            tile_idx = (gg - 1) * n_graph_cols + aa;

            blankAxes = nexttile(groupedSpdLayout, tile_idx);
            axis(blankAxes, 'off');

            blankAxes = nexttile(groupedExponentLayout, tile_idx);
            axis(blankAxes, 'off');

            blankAxes = nexttile(groupedVarianceLayout, tile_idx);
            axis(blankAxes, 'off');
        end
    end

    % Add bold super-titles to all grouped figures
    figure(groupedExponentHandle);
    sgtitle(sprintf('%s | Exponent Maps Grouped By Subject', options.title), ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    figure(groupedVarianceHandle);
    sgtitle(sprintf('%s | Variance Maps Grouped By Subject', options.title), ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    figure(groupedSpdHandle);
    sgtitle(sprintf('%s | SPDs Grouped By Subject', options.title), ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    % Export grouped figures if requested
    if (~strcmp(options.output_dir, ""))
        assert(~strcmp(options.title, ""), ...
            'groupSPDs:MissingTitle', ...
            'options.title must be provided when saving output.');

        exponent_output_filepath = fullfile(options.output_dir, sprintf('%s_groupedExponentMaps.pdf', options.title));
        variance_output_filepath = fullfile(options.output_dir, sprintf('%s_groupedVarianceMaps.pdf', options.title));
        spd_output_filepath = fullfile(options.output_dir, sprintf('%s_groupedSPDs.pdf', options.title));

        if (~exist(exponent_output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(groupedExponentHandle, exponent_output_filepath, 'ContentType', 'vector');
        end

        if (~exist(variance_output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(groupedVarianceHandle, variance_output_filepath, 'ContentType', 'vector');
        end

        if (~exist(spd_output_filepath, 'file') || options.overwrite_existing)
            exportgraphics(groupedSpdHandle, spd_output_filepath, 'ContentType', 'vector');
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

function iCopyEntireMapFigureFromPlotSPDs(sourceFigureHandle, targetFigureHandle, tilePosition, activity_name, group_name, is_first_column)
% Copy the FULL map layout from a plotSPDs figure directly into a target
% figure region corresponding to one outer tile.
%
% Unlike the earlier uipanel-based approach, this recreates everything as
% ordinary axes/colorbar graphics objects parented directly to the figure,
% which is much more reliable for exportgraphics.

    % Find all real plotting axes in the source figure
    sourceAxes = findall(sourceFigureHandle, 'Type', 'axes');

    keepMask = true(size(sourceAxes));
    for ii = 1:numel(sourceAxes)
        thisTag = get(sourceAxes(ii), 'Tag');
        if strcmpi(thisTag, 'Colorbar') || strcmpi(thisTag, 'legend')
            keepMask(ii) = false;
        end
    end
    sourceAxes = sourceAxes(keepMask);

    % Sort source axes left-to-right
    if (numel(sourceAxes) > 1)
        xPositions = zeros(numel(sourceAxes), 1);
        for ii = 1:numel(sourceAxes)
            thisPos = get(sourceAxes(ii), 'Position');
            xPositions(ii) = thisPos(1);
        end
        [~, sortIdx] = sort(xPositions, 'ascend');
        sourceAxes = sourceAxes(sortIdx);
    end

    % ---------------------------------------------------------------------
    % Layout inside the tile
    % ---------------------------------------------------------------------
    tileLeft   = tilePosition(1);
    tileBottom = tilePosition(2);
    tileWidth  = tilePosition(3);
    tileHeight = tilePosition(4);

    % Reserve a small strip at top for row/activity labels
    titleFrac = 0.10;
    gapFrac = 0.02;

    contentBottom = tileBottom;
    contentHeight = tileHeight * (1 - titleFrac);
    titleBottom = tileBottom + contentHeight;
    titleHeight = tileHeight * titleFrac;

    % Add invisible annotation axes for labels
    annotationAxes = axes( ...
        'Parent', targetFigureHandle, ...
        'Units', 'normalized', ...
        'Position', [tileLeft, titleBottom, tileWidth, titleHeight], ...
        'Visible', 'off', ...
        'HitTest', 'off');

    if (is_first_column)
        text(annotationAxes, 0.00, 0.5, sprintf('%s', group_name), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    text(annotationAxes, 0.5, 0.5, activity_name, ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');

    % ---------------------------------------------------------------------
    % Build either a 1x2 inner layout or a single full-size map
    % ---------------------------------------------------------------------
    if (numel(sourceAxes) > 1)
        innerGap = tileWidth * gapFrac;
        innerWidth = (tileWidth - innerGap) / 2;

        leftPos = [tileLeft, contentBottom, innerWidth, contentHeight];
        rightPos = [tileLeft + innerWidth + innerGap, contentBottom, innerWidth, contentHeight];

        iCopySingleAxesWithColorbar(sourceAxes(1), targetFigureHandle, leftPos);
        iCopySingleAxesWithColorbar(sourceAxes(2), targetFigureHandle, rightPos);
    else
        fullPos = [tileLeft, contentBottom, tileWidth, contentHeight];
        iCopySingleAxesWithColorbar(sourceAxes(1), targetFigureHandle, fullPos);
    end
end

function iCopySingleAxesWithColorbar(sourceAxes, targetFigureHandle, targetPosition)
% Recreate one map axes and its colorbar directly in the destination figure.

    % Create destination axes
    targetAxes = axes( ...
        'Parent', targetFigureHandle, ...
        'Units', 'normalized', ...
        'Position', targetPosition);

    % Copy children
    childrenToCopy = allchild(sourceAxes);
    copyobj(childrenToCopy, targetAxes);

    % Copy important visual properties
    set(targetAxes, ...
        'XLim', get(sourceAxes, 'XLim'), ...
        'YLim', get(sourceAxes, 'YLim'), ...
        'CLim', get(sourceAxes, 'CLim'), ...
        'DataAspectRatio', get(sourceAxes, 'DataAspectRatio'), ...
        'PlotBoxAspectRatio', get(sourceAxes, 'PlotBoxAspectRatio'), ...
        'XDir', get(sourceAxes, 'XDir'), ...
        'YDir', get(sourceAxes, 'YDir'), ...
        'FontSize', get(sourceAxes, 'FontSize'), ...
        'Box', get(sourceAxes, 'Box'), ...
        'Visible', get(sourceAxes, 'Visible'));

    colormap(targetAxes, colormap(sourceAxes));

    % Copy the source map title over
    sourceTitle = get(sourceAxes, 'Title');
    title(targetAxes, sourceTitle.String, ...
        'Interpreter', sourceTitle.Interpreter, ...
        'FontWeight', sourceTitle.FontWeight, ...
        'FontSize', sourceTitle.FontSize, ...
        'Color', sourceTitle.Color);

    % Keep axis labels suppressed in the grouped mini-maps
    xlabel(targetAxes, '');
    ylabel(targetAxes, '');

    % Recreate colorbar as a normal figure graphics object
    cb = colorbar(targetAxes);

    % Keep the axes square-ish for map display
    axis(targetAxes, 'square');

    % Slightly shrink axes width so colorbar has room
    axPos = targetAxes.Position;
    axPos(3) = axPos(3) * 0.82;
    targetAxes.Position = axPos;

    % Reposition colorbar neatly inside this tile region
    cb.Units = 'normalized';
    cb.Position = [axPos(1) + axPos(3) + 0.01 * targetPosition(3), ...
                   axPos(2), ...
                   0.04 * targetPosition(3), ...
                   axPos(4)];
end

function groups_across_participant_std = find_across_participant_std(groupedActivityData, n_participants)
    % First, find the group names 
    group_names = fieldnames(groupedActivityData); 

    % Initialize map to save averages by group struct 
    groups_across_participant_std = containers.Map(); 

    % Go over the group names 
    for nn = 1:numel(group_names)
        group_name = group_names{nn};

        % Initialize containers for the data we will append over for this group 
        virtuallyFoveated_groupAvgs_center = {}; 
        virtuallyFoveated_groupAvgs_periphery = {}; 

        justProjection_groupAvgs_center = {}; 
        justProjection_groupAvgs_periphery = {};

        % Go over each of the activity in this group 
        group_activities = fieldnames(groupedActivityData.(group_name)); 
        for aa = 1:numel(group_activities)
            activity_name = group_activities{aa};
            
            % For each activity, there may be 1 or more projection types 
            projection_types = fieldnames(groupedActivityData.(group_name).(activity_name)); 
            
            % Start accumulating 
            if(ismember("virtuallyFoveated", projection_types))
                % Load in the virtually foveated data 
                virtually_foveated_data = load(groupedActivityData.(group_name).(activity_name).virtuallyFoveated).activityData.(activity_name); 

                % Iterate over the subjects we have for this activity 

                num_subjects = n_participants; 
                for ss = 1:num_subjects
                    center_region_avg = virtually_foveated_data.regionAveragesAcrossSubjects{ss}.center; 
                    periphery_region_avg = virtually_foveated_data.regionAveragesAcrossSubjects{ss}.periphery;
                
                    % Save these 
                    virtuallyFoveated_groupAvgs_center{end+1} = center_region_avg;
                    virtuallyFoveated_groupAvgs_periphery{end+1} = periphery_region_avg; 

                end 

       
            end 

            if(ismember("justProjection", projection_types))
                just_projection_data = load(groupedActivityData.(group_name).(activity_name).justProjection).activityData.(activity_name); 

                % Iterate over the subjects we have for this activity 
                num_subjects = n_participants; 
                for ss = 1:num_subjects
                    center_region_avg = just_projection_data.regionAveragesAcrossSubjects{ss}.center; 
                    periphery_region_avg = just_projection_data.regionAveragesAcrossSubjects{ss}.periphery;

                    justProjection_groupAvgs_center{end+1} = center_region_avg;
                    justProjection_groupAvgs_periphery{end+1} = periphery_region_avg; 
                end 
            end 
        end 

        % Now let's compute and save the STDs 
        group_std_struct = struct(); 
        if(numel(virtuallyFoveated_groupAvgs_center) > 0)
            center_mat = horzcat(virtuallyFoveated_groupAvgs_center{:}).';
            periphery_mat = horzcat(virtuallyFoveated_groupAvgs_periphery{:}).';

            center_std = std(center_mat, 0, 1, 'omitnan');
            periphery_std = std(periphery_mat, 0, 1, 'omitnan');

            group_std_struct.virtuallyFoveated.center = center_std;
            group_std_struct.virtuallyFoveated.periphery = periphery_std;
        end

        if(numel(justProjection_groupAvgs_center) > 0) 
            center_mat = horzcat(justProjection_groupAvgs_center{:}).';
            periphery_mat = horzcat(justProjection_groupAvgs_periphery{:}).';

            center_std = std(center_mat, 0, 1, 'omitnan');
            periphery_std = std(periphery_mat, 0, 1, 'omitnan');

            group_std_struct.justProjection.center = center_std; 
            group_std_struct.justProjection.periphery = periphery_std; 
        end 

        groups_across_participant_std(char(group_name)) = group_std_struct; 

    end 

    return ; 
end 