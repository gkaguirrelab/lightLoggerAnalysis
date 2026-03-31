function virtuallyFoveatedActivityDataAcrossSubjects = processSPDAcrossActivities(input_dir, output_dir, options)
% Aggregate per-activity temporal SPD analysis results across all activities
%
% Syntax:
%   virtuallyFoveatedActivityDataAcrossSubjects = processSPDAcrossActivities(input_dir, output_dir)
%   virtuallyFoveatedActivityDataAcrossSubjects = processSPDAcrossActivities(input_dir, output_dir, options)
%
% Description:
%   This function searches a directory of activity folders and aggregates
%   previously computed SPD analysis results across activities. Activity
%   folders are expected to be direct subdirectories of input_dir, and each
%   activity folder should contain a saved "*SPDResultsAcrossSubjects.mat"
%   file representing the SPD results already averaged across subjects.
%
%   The function loads these per-activity SPD outputs and stacks them across
%   activities. It then computes the mean exponent map, variance map,
%   regional SPDs, and median images across all activities. The aggregated
%   results are stored under a single field:
%
%       virtuallyFoveatedActivityDataAcrossSubjects.acrossAll
%
%   The function can optionally save the aggregated results to disk and
%   export summary figures showing the across-activity exponent map,
%   variance map, and regional SPD summaries.
%
% Inputs:
%   input_dir             - Char/string. Path to the directory containing
%                           activity folders.
%   output_dir            - Char/string. Path to the directory where
%                           aggregated results and optional figures should
%                           be saved. If empty, results are returned but
%                           not written to disk.
%
% Optional key/value pairs:
%   activities            - Cell array. List of activity folder names to
%                           include. If empty, all activity folders in
%                           input_dir are used.
%   verbose               - Logical. If true, print progress information
%                           while processing activities.
%   save_figures          - Logical. If true, generate and export summary
%                           figures to output_dir.
%   overwrite_existing    - Logical. If false, skip processing when a saved
%                           across-all results file already exists.
%   projection_type       - String. Either "justProjection" or
%                           "virtuallyFoveated".
%   fovDegrees            - Scalar. Field of view in degrees used for
%                           plotting aggregated SPD maps.
%   exponent_clim         - Passed through to plotSPDs.
%   variance_clim         - Passed through to plotSPDs.
%   spd_xlim              - Passed through to plotSPDs.
%   spd_ylim              - Passed through to plotSPDs.
%   combine_figures       - Logical. If true, also load justProjection data
%                           and pass it into plotSPDs for comparison plots.
%
% Outputs:
%   virtuallyFoveatedActivityDataAcrossSubjects
%                         - Struct. Structure containing across-activity SPD
%                           analysis results stored under:
%
%                               .acrossAll
%
%                           This field contains:
%                               .exponentMaps
%                               .varianceMaps
%                               .spdByRegions
%                               .medianImages
%                               .frameDropVector
%                               .frq
%                               .exponentMap
%                               .varianceMap
%                               .spdByRegion
%                               .medianImage
%                               .activityNames
%
% Example:
%{
    input_dir = "/path/to/SPD_by_activity";
    output_dir = "/path/to/SPD_across_all";

    virtuallyFoveatedActivityDataAcrossSubjects = processSPDAcrossActivities( ...
        input_dir, ...
        output_dir, ...
        "activities", {'walkIndoor', 'chat', 'work'}, ...
        "verbose", true, ...
        "save_figures", true ...
    );
%}

    arguments
        input_dir
        output_dir
        options.activities = {}
        options.verbose = false
        options.save_figures = false
        options.overwrite_existing = false
        options.projection_type {mustBeMemberString(options.projection_type, ["justProjection", "virtuallyFoveated"])} = "virtuallyFoveated"
        options.fovDegrees = 120
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.combine_figures = false
    end

    % Store frequently reused options in local variables for readability.
    fovDegrees = options.fovDegrees;
    combine_figures = options.combine_figures;

    % ---------------------------------------------------------------------
    % STEP 0: Initialize output structures
    % ---------------------------------------------------------------------
    % The final output is stored under:
    %
    %   virtuallyFoveatedActivityDataAcrossSubjects.acrossAll
    %
    % This mirrors the same indexing style you use elsewhere in the SPD
    % pipeline, except here the "activity name" is a single synthetic field
    % called "acrossAll".
    virtuallyFoveatedActivityDataAcrossSubjects = struct();
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll = struct();

    % Build results into a temporary struct, then assign it back into
    % .acrossAll once all fields are populated.
    dataStruct = virtuallyFoveatedActivityDataAcrossSubjects.acrossAll;

    % If combined figures are requested, maintain a second comparison struct
    % containing across-activity justProjection results.
    justProjectionDataAcrossAll = false;
    if (combine_figures)
        justProjectionDataAcrossAll = struct();
    end

    % ---------------------------------------------------------------------
    % STEP 1: Find all activity folders directly under input_dir
    % ---------------------------------------------------------------------
    % We assume input_dir contains one folder per activity.
    activityFiles = dir(input_dir);

    % Keep only directories.
    activityFiles = activityFiles([activityFiles.isdir]);

    % Remove "." and "..", since these are not real activity folders.
    activityNames = {activityFiles.name};
    isRealFolder = ~ismember(activityNames, {'.', '..'});
    activityFiles = activityFiles(isRealFolder);

    % If the caller specified a subset of activities, restrict to those.
    if (~isempty(options.activities))
        keepActivities = ismember({activityFiles.name}, options.activities);
        activityFiles = activityFiles(keepActivities);
    end

    % Count how many activity folders we will process.
    nActivities = length(activityFiles);

    % Save activity names into the output so the final struct records which
    % folders were included in the aggregation.
    for aa = 1:nActivities
        dataStruct.activityNames{aa} = activityFiles(aa).name;
    end

    % ---------------------------------------------------------------------
    % STEP 2: Check whether final output already exists
    % ---------------------------------------------------------------------
    % If overwrite_existing is false and the saved across-all result already
    % exists, skip the computation entirely.
    if (output_dir ~= "")
        output_filepath = fullfile(output_dir, ...
            sprintf("allActivities_%s_SPDResultsAcrossAll.mat", options.projection_type));

        if (~options.overwrite_existing && isfile(output_filepath))
            if (options.verbose)
                fprintf("Skipping because output already exists: %s\n", output_filepath);
            end
            return;
        end
    end

    % ---------------------------------------------------------------------
    % STEP 3: Loop over all activity folders and load each activity's
    %         SPDResultsAcrossSubjects file
    % ---------------------------------------------------------------------
    for aa = 1:nActivities

        if (options.verbose)
            fprintf("Activity: %d/%d\n", aa, nActivities);
        end

        % Retrieve the current activity's name and full path.
        activityName = activityFiles(aa).name;
        activityPath = fullfile(activityFiles(aa).folder, activityName);

        % Search within this activity folder for the requested projection
        % type's already-across-subjects SPD result.
        files = dir(fullfile(activityPath, ...
            sprintf('*%s_SPDResultsAcrossSubjects.mat', options.projection_type)));

        % If no file exists for this activity, warn and skip it.
        if isempty(files)
            fprintf('No SPDResultsAcrossSubjects file for %s\n', activityName);
            continue;
        end

        % Build the full path to the first matching file.
        thisFile = fullfile(files(1).folder, files(1).name);

        % Load the saved .mat file.
        loadedStruct = load(thisFile);

        % Each file is expected to contain a top-level struct named
        % activityData, and inside that a single field storing the actual
        % SPD result struct. To remain robust to exact field naming, we
        % dynamically retrieve the first field.
        fieldNames = fieldnames(loadedStruct.activityData);
        dataField = fieldNames{1};

        % Extract the SPD struct for this activity.
        spd = loadedStruct.activityData.(dataField);

        % -------------------------------------------------------------
        % Stack the per-activity results
        % -------------------------------------------------------------
        % We stack each activity along the final dimension so that later we
        % can compute a mean across activities.
        dataStruct.exponentMaps(:,:,aa) = spd.exponentMap;
        dataStruct.varianceMaps(:,:,aa) = spd.varianceMap;
        dataStruct.spdByRegions(:,:,:,aa) = spd.spdByRegion;
        dataStruct.medianImages(:,:,aa) = spd.medianImage;
        dataStruct.frameDropVector{aa} = spd.frameDropVector;

        % Frequency support is expected to be the same across all activity
        % files, so we simply keep the last loaded one.
        dataStruct.frq = spd.frq;

        % -------------------------------------------------------------
        % Optionally load justProjection data for comparison figures
        % -------------------------------------------------------------
        if (combine_figures)
            jp_files = dir(fullfile(activityPath, ...
                sprintf('*%s_SPDResultsAcrossSubjects.mat', "justProjection")));

            if isempty(jp_files)
                fprintf('No justProjection SPDResultsAcrossSubjects file for %s\n', activityName);
            else
                jpFile = fullfile(jp_files(1).folder, jp_files(1).name);
                jpStruct = load(jpFile);

                jpFieldNames = fieldnames(jpStruct.activityData);
                jpField = jpFieldNames{1};

                jp = jpStruct.activityData.(jpField);

                justProjectionDataAcrossAll.exponentMaps(:,:,aa) = jp.exponentMap;
                justProjectionDataAcrossAll.varianceMaps(:,:,aa) = jp.varianceMap;
                justProjectionDataAcrossAll.spdByRegions(:,:,:,aa) = jp.spdByRegion;
                justProjectionDataAcrossAll.medianImages(:,:,aa) = jp.medianImage;
                justProjectionDataAcrossAll.frameDropVector{aa} = jp.frameDropVector;
                justProjectionDataAcrossAll.frq = jp.frq;
            end
        end
    end

    % ---------------------------------------------------------------------
    % STEP 4: Compute averages across all activities
    % ---------------------------------------------------------------------
    % The stacked exponent/variance/median arrays use the 3rd dimension for
    % activity, while spdByRegions uses the 4th dimension for activity.
    %
    % We use 'omitmissing' to avoid missing/NaN entries contaminating the
    % average when some activities were absent or incomplete.
    dataStruct.exponentMap = squeeze(mean(dataStruct.exponentMaps, 3, 'omitmissing'));
    dataStruct.varianceMap = squeeze(mean(dataStruct.varianceMaps, 3, 'omitmissing'));
    dataStruct.spdByRegion = squeeze(mean(dataStruct.spdByRegions, 4, 'omitmissing'));
    dataStruct.medianImage = squeeze(mean(dataStruct.medianImages, 3, 'omitmissing'));

    % If comparison data was loaded, average that across activities too.
    if (combine_figures && isstruct(justProjectionDataAcrossAll))
        justProjectionDataAcrossAll.exponentMap = squeeze(mean(justProjectionDataAcrossAll.exponentMaps, 3, 'omitmissing'));
        justProjectionDataAcrossAll.varianceMap = squeeze(mean(justProjectionDataAcrossAll.varianceMaps, 3, 'omitmissing'));
        justProjectionDataAcrossAll.spdByRegion = squeeze(mean(justProjectionDataAcrossAll.spdByRegions, 4, 'omitmissing'));
        justProjectionDataAcrossAll.medianImage = squeeze(mean(justProjectionDataAcrossAll.medianImages, 3, 'omitmissing'));
    end

    % Assign the completed result struct back into the final output.
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll = dataStruct;

    % ---------------------------------------------------------------------
    % STEP 5: Save results and optionally export figures
    % ---------------------------------------------------------------------
    if (output_dir ~= "")
        % Create the output directory if needed.
        if (~isfolder(output_dir))
            mkdir(output_dir);
        end

        % Save the full wrapper struct so the saved file preserves the
        % ".acrossAll" indexing convention.
        activityData = virtuallyFoveatedActivityDataAcrossSubjects;
        save(output_filepath, 'activityData');

        % If requested, generate summary plots and export them.
        if (options.save_figures)

            % plotSPDs expects a struct containing SPD fields such as
            % exponentMap, varianceMap, spdByRegion, etc. Remove metadata
            % fields like activityNames before plotting.
            plotData = dataStruct;
            if (isfield(plotData, 'activityNames'))
                plotData = rmfield(plotData, 'activityNames');
            end

            % Prepare justProjection comparison data, if requested.
            justProjectionPlotData = justProjectionDataAcrossAll;
            if (combine_figures && isstruct(justProjectionPlotData) && ...
                    isfield(justProjectionPlotData, 'activityNames'))
                justProjectionPlotData = rmfield(justProjectionPlotData, 'activityNames');
            end

            % Generate the SPD summary figures. Use descriptive figure-handle
            % names so they match the rest of your codebase.
            if (combine_figures && isstruct(justProjectionPlotData))
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs( ...
                    plotData, ...
                    "fovDegrees", fovDegrees, ...
                    "exponent_clim", options.exponent_clim, ...
                    "variance_clim", options.variance_clim, ...
                    "spd_xlim", options.spd_xlim, ...
                    "spd_ylim", options.spd_ylim, ...
                    "justProjectionActivityData", justProjectionPlotData ...
                );
            else
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs( ...
                    plotData, ...
                    "fovDegrees", fovDegrees, ...
                    "exponent_clim", options.exponent_clim, ...
                    "variance_clim", options.variance_clim, ...
                    "spd_xlim", options.spd_xlim, ...
                    "spd_ylim", options.spd_ylim ...
                );
            end

            % Build the output paths for each exported figure.
            exponentMapPath = fullfile(output_dir, ...
                sprintf('allActivities_%s_exponentMapAcrossAll.pdf', options.projection_type));
            varianceMapPath = fullfile(output_dir, ...
                sprintf('allActivities_%s_varianceMapAcrossAll.pdf', options.projection_type));
            spdByRegionPath = fullfile(output_dir, ...
                sprintf('allActivities_%s_spdByRegionAcrossAll.pdf', options.projection_type));

            % Only overwrite exported figures if overwrite_existing is true,
            % matching the behavior of the across-subjects version.
            if (~exist(exponentMapPath, 'file') || options.overwrite_existing)
                exportgraphics(exponentMapHandle, exponentMapPath, 'ContentType', 'vector');
            end

            if (~exist(varianceMapPath, 'file') || options.overwrite_existing)
                exportgraphics(varianceMapHandle, varianceMapPath, 'ContentType', 'vector');
            end

            if (~exist(spdByRegionPath, 'file') || options.overwrite_existing)
                exportgraphics(spdByRegionHandle, spdByRegionPath, 'ContentType', 'vector');
            end

            % Close figure windows so repeated calls do not accumulate open
            % figures.
            close([exponentMapHandle varianceMapHandle spdByRegionHandle]);
        end
    end
end


function mustBeMemberString(value, allowed)
% mustBeMemberString
%
% Local custom validator for MATLAB argument blocks.
%
% This helper converts both the input and the list of allowed values to
% strings, then checks whether the input matches one of the permitted
% values.

    valueStr = string(value);
    allowedStr = string(allowed);

    if ~any(valueStr == allowedStr)
        error("Value must be one of: %s. Got: %s", ...
            strjoin(allowedStr, ", "), valueStr);
    end
end