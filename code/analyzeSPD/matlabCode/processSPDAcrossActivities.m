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
        options.n_participants = 1; 
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
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .activityNames{aa} = activityFiles(aa).name;
    end

    % ---------------------------------------------------------------------
    % STEP 2: Check whether final output already exists
    % ---------------------------------------------------------------------
    % If overwrite_existing is false and the saved across-all result already
    % exists, skip the computation entirely.
    if (output_dir ~= "")
        output_filepath = fullfile(output_dir, ...
            sprintf("acrossActivities_%s_SPDResultsAcrossAll.mat", options.projection_type));

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
    validActivityIdx = 0;

    for aa = 1:nActivities

        if (options.verbose)
            fprintf("Activity candidate: %d/%d\n", aa, nActivities);
        end

        activityName = activityFiles(aa).name;
        activityPath = fullfile(activityFiles(aa).folder, activityName);

        files = dir(fullfile(activityPath, ...
            sprintf('*%s_SPDResultsAcrossSubjects.mat', options.projection_type)));

        if (options.verbose)
            fprintf('Checking activity folder: %s\n', activityPath);
            fprintf('Looking for pattern: *%s_SPDResultsAcrossSubjects.mat\n', string(options.projection_type));
            fprintf('Found %d matching files\n', length(files));
        end

        if isempty(files)
            fprintf('No SPDResultsAcrossSubjects file for %s\n', activityName);
            continue;
        end

        thisFile = fullfile(files(1).folder, files(1).name);
        loadedStruct = load(thisFile);

        if (~isfield(loadedStruct, 'activityData'))
            fprintf('Skipping %s because variable "activityData" was not found\n', thisFile);
            continue;
        end

        fieldNames = fieldnames(loadedStruct.activityData);
        if isempty(fieldNames)
            fprintf('Skipping %s because activityData had no fields\n', thisFile);
            continue;
        end

        dataField = fieldNames{1};
        spd = loadedStruct.activityData.(dataField);

        if (~isfield(spd, 'exponentMap') || ~isfield(spd, 'varianceMap') || ...
                ~isfield(spd, 'spdByRegion') || ~isfield(spd, 'medianImage') || ...
                ~isfield(spd, 'frameDropVector') || ~isfield(spd, 'frq'))
            fprintf('Skipping %s because required SPD fields were missing\n', thisFile);
            continue;
        end

        validActivityIdx = validActivityIdx + 1;

        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .activityNames{validActivityIdx} = activityName;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .exponentMaps(:,:,validActivityIdx) = spd.exponentMap;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .varianceMaps(:,:,validActivityIdx) = spd.varianceMap;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .spdByRegions(:,:,:,validActivityIdx) = spd.spdByRegion;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .medianImages(:,:,validActivityIdx) = spd.medianImage;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .frameDropVector{validActivityIdx} = spd.frameDropVector;
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .frq = spd.frq;

        if (combine_figures)
            jp_files = dir(fullfile(activityPath, ...
                sprintf('*%s_SPDResultsAcrossSubjects.mat', "justProjection")));

            if isempty(jp_files)
                fprintf('No justProjection SPDResultsAcrossSubjects file for %s\n', activityName);
            else
                jpFile = fullfile(jp_files(1).folder, jp_files(1).name);
                jpStruct = load(jpFile);

                if (isfield(jpStruct, 'activityData'))
                    jpFieldNames = fieldnames(jpStruct.activityData);
                    if ~isempty(jpFieldNames)
                        jpField = jpFieldNames{1};
                        jp = jpStruct.activityData.(jpField);

                        if (isfield(jp, 'exponentMap') && isfield(jp, 'varianceMap') && ...
                                isfield(jp, 'spdByRegion') && isfield(jp, 'medianImage') && ...
                                isfield(jp, 'frameDropVector') && isfield(jp, 'frq'))

                            justProjectionDataAcrossAll.exponentMaps(:,:,validActivityIdx) = jp.exponentMap;
                            justProjectionDataAcrossAll.varianceMaps(:,:,validActivityIdx) = jp.varianceMap;
                            justProjectionDataAcrossAll.spdByRegions(:,:,:,validActivityIdx) = jp.spdByRegion;
                            justProjectionDataAcrossAll.medianImages(:,:,validActivityIdx) = jp.medianImage;
                            justProjectionDataAcrossAll.frameDropVector{validActivityIdx} = jp.frameDropVector;
                            justProjectionDataAcrossAll.frq = jp.frq;
                        end
                    end
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % STEP 3.5: Handle case where nothing valid was loaded
    % ---------------------------------------------------------------------
    if (validActivityIdx == 0)
        warning('No valid SPDResultsAcrossSubjects files were loaded from input_dir: %s', input_dir);
        virtuallyFoveatedActivityDataAcrossSubjects.acrossAll = virtuallyFoveatedActivityDataAcrossSubjects.acrossAll ;
        return;
    end

    % ---------------------------------------------------------------------
    % STEP 4: Compute averages across all successfully loaded activities
    % ---------------------------------------------------------------------
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .exponentMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .exponentMaps(:,:,1:validActivityIdx), 3, 'omitmissing'));
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .varianceMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .varianceMaps(:,:,1:validActivityIdx), 3, 'omitmissing'));
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .spdByRegion = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .spdByRegions(:,:,:,1:validActivityIdx), 4, 'omitmissing'));
    virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .medianImage = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.acrossAll .medianImages(:,:,1:validActivityIdx), 3, 'omitmissing'));

    if (combine_figures && isstruct(justProjectionDataAcrossAll) && ...
            isfield(justProjectionDataAcrossAll, 'exponentMaps'))
        justProjectionDataAcrossAll.exponentMap = squeeze(mean(justProjectionDataAcrossAll.exponentMaps, 3, 'omitmissing'));
        justProjectionDataAcrossAll.varianceMap = squeeze(mean(justProjectionDataAcrossAll.varianceMaps, 3, 'omitmissing'));
        justProjectionDataAcrossAll.spdByRegion = squeeze(mean(justProjectionDataAcrossAll.spdByRegions, 4, 'omitmissing'));
        justProjectionDataAcrossAll.medianImage = squeeze(mean(justProjectionDataAcrossAll.medianImages, 3, 'omitmissing'));
    end

    % ---------------------------------------------------------------------
    % STEP 5: Save results and optionally export figures
    % ---------------------------------------------------------------------
    if (output_dir ~= "")
        if (~isfolder(output_dir))
            mkdir(output_dir);
        end

        activityData = virtuallyFoveatedActivityDataAcrossSubjects;
        save(output_filepath, 'activityData');

        if (options.save_figures)

            % plotSPDs expects a wrapper struct whose first field name is the
            % activity name. For this across-all case, that synthetic activity
            % name is "acrossAll". Therefore, we must pass wrapper structs of
            % the form:
            %
            %   plotData.acrossAll = ...
            %
            % rather than passing the bare SPD struct directly.

            % If requested, wrap the justProjection comparison data in the same
            % way so plotSPDs can index both datasets using the same field name.
            justProjectionPlotData = false;
            if (combine_figures && isstruct(justProjectionDataAcrossAll) && ...
                    isfield(justProjectionDataAcrossAll, 'exponentMap'))
                justProjectionPlotData = struct();
                justProjectionPlotData.acrossAll = justProjectionDataAcrossAll;
            end

            % Generate the SPD summary figures. Use descriptive figure-handle
            % names to match the rest of the pipeline.
            if (~islogical(justProjectionPlotData))
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs( ...
                    virtuallyFoveatedActivityDataAcrossSubjects, ...
                    "fovDegrees", fovDegrees, ...
                    "exponent_clim", options.exponent_clim, ...
                    "variance_clim", options.variance_clim, ...
                    "spd_xlim", options.spd_xlim, ...
                    "spd_ylim", options.spd_ylim, ...
                    "justProjectionActivityData", justProjectionPlotData, ...
                    "num_participants", options.n_participants...
                );
            else
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs( ...
                    virtuallyFoveatedActivityDataAcrossSubjects, ...
                    "fovDegrees", fovDegrees, ...
                    "exponent_clim", options.exponent_clim, ...
                    "variance_clim", options.variance_clim, ...
                    "spd_xlim", options.spd_xlim, ...
                    "spd_ylim", options.spd_ylim, ...
                    "num_participants", options.n_participants...
                );
            end

            exponentMapPath = fullfile(output_dir, ...
                sprintf('acrossActivities_%s_exponentMapAcrossAll.pdf', options.projection_type));
            varianceMapPath = fullfile(output_dir, ...
                sprintf('acrossActivities_%s_varianceMapAcrossAll.pdf', options.projection_type));
            spdByRegionPath = fullfile(output_dir, ...
                sprintf('acrossActivities_%s_spdByRegionAcrossAll.pdf', options.projection_type));

            if (~exist(exponentMapPath, 'file') || options.overwrite_existing)
                exportgraphics(exponentMapHandle, exponentMapPath, 'ContentType', 'vector');
            end

            if (~exist(varianceMapPath, 'file') || options.overwrite_existing)
                exportgraphics(varianceMapHandle, varianceMapPath, 'ContentType', 'vector');
            end

            if (~exist(spdByRegionPath, 'file') || options.overwrite_existing)
                exportgraphics(spdByRegionHandle, spdByRegionPath, 'ContentType', 'vector');
            end

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