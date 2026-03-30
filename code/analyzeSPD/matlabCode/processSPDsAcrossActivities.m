function virtuallyFoveatedActivityDataAcrossActivities = processSPDsAcrossActivities(input_dir, output_dir, options)
% processSPDsAcrossActivities
%
% Aggregate temporal spatial power spectrum (SPD) results across activities.
%
% This function is designed to mirror the behavior of the corresponding
% across-subjects function, but instead of starting from subject folders,
% it starts from activity folders.
%
% Expected directory structure:
%   input_dir/
%       activityA/
%           ... many files ...
%           *_virtuallyFoveated_SPDResultsAcrossSubjects.mat
%           *_justProjection_SPDResultsAcrossSubjects.mat
%       activityB/
%           ...
%
% For each activity folder, this function loads the saved SPD results that
% were already averaged across subjects, and then averages those results
% across activities.
%
% In other words:
%   each activity folder contains: "across-subjects" SPD outputs
%   this function produces:        "across-activities" SPD outputs
%
% The output naming convention intentionally uses "AcrossActivities" in:
%   - the saved .mat file
%   - the exported figures
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% input_dir
%   Path to a folder whose direct subdirectories are activity folders.
%
% output_dir
%   Path where the final aggregated .mat file and figures should be saved.
%
% options
%   Name-value options controlling what gets loaded, plotted, and saved.
%
% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
% options.activities
%   Optional list of activity folder names to include. If empty, all
%   activity folders under input_dir are used.
%
% options.verbose
%   If true, print progress messages while processing.
%
% options.save_figures
%   If true, generate and export figures summarizing the averaged SPD data.
%
% options.overwrite_existing
%   If false, skip computation if the final saved result already exists.
%
% options.projection_type
%   Which projection type to load as the main data source:
%       "virtuallyFoveated" or "justProjection"
%
% options.fovDegrees
%   Passed to plotSPDs for display/plotting.
%
% options.exponent_clim
% options.variance_clim
% options.spd_xlim
% options.spd_ylim
%   Optional plotting limits passed directly to plotSPDs.
%
% options.combine_figures
%   If true, also load justProjection data and pass it to plotSPDs for
%   side-by-side comparison plots.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% virtuallyFoveatedActivityDataAcrossActivities
%   Struct containing the averaged-across-activities SPD outputs.
%
%   This struct includes:
%       exponentMaps     - all activity exponent maps stacked together
%       varianceMaps     - all activity variance maps stacked together
%       spdByRegions     - all activity regional SPD results stacked together
%       medianImages     - all activity median images stacked together
%       frameDropVector  - frame drop information for each activity
%       frq              - frequency vector
%       exponentMap      - average exponent map across activities
%       varianceMap      - average variance map across activities
%       spdByRegion      - average regional SPD summary across activities
%       medianImage      - average median image across activities
%       activityNames    - list of included activity names

    arguments
        input_dir
        output_dir
        options.activities = {}   % Optional list of activity names to keep
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

    % Store frequently used option values in local variables for readability.
    fovDegrees = options.fovDegrees;
    combine_figures = options.combine_figures;

    % Initialize the main output struct that will hold:
    %   1) the per-activity loaded results
    %   2) the final averaged-across-activities results
    virtuallyFoveatedActivityDataAcrossActivities = struct();

    % If we are making combined plots, initialize the corresponding
    % justProjection struct too. Otherwise keep it as false.
    justProjectionActivityDataAcrossActivities = false;
    if (combine_figures)
        justProjectionActivityDataAcrossActivities = struct();
    end

    % ---------------------------------------------------------------------
    % STEP 1: Find all activity folders directly under input_dir
    % ---------------------------------------------------------------------
    % We assume input_dir contains one folder per activity.
    activityFiles = dir(input_dir);

    % Keep only directory entries.
    activityFiles = activityFiles([activityFiles.isdir]);

    % Remove "." and "..", since those are not real activity folders.
    activityNames = {activityFiles.name};
    isRealFolder = ~ismember(activityNames, {'.', '..'});
    activityFiles = activityFiles(isRealFolder);

    % If the caller specified a subset of activities, keep only those.
    if (~isempty(options.activities))
        keepActivities = ismember({activityFiles.name}, options.activities);
        activityFiles = activityFiles(keepActivities);
    end

    % Count how many activity folders we will process.
    nActivities = length(activityFiles);

    % Save the activity names into the output struct so the final saved
    % results remember which activities were included.
    for aa = 1:nActivities
        virtuallyFoveatedActivityDataAcrossActivities.activityNames{aa} = activityFiles(aa).name;
    end

    % ---------------------------------------------------------------------
    % STEP 2: Check whether final output already exists
    % ---------------------------------------------------------------------
    % If overwrite_existing is false and the final .mat file is already
    % present, then skip the computation entirely.
    if (output_dir ~= "")
        output_filepath = fullfile(output_dir, sprintf("allActivities_%s_SPDResultsAcrossActivities.mat", options.projection_type));

        if (~options.overwrite_existing && isfile(output_filepath))
            if (options.verbose)
                fprintf("Skipping because output already exists: %s\n", output_filepath);
            end
            return;
        end
    end

    % ---------------------------------------------------------------------
    % STEP 3: Loop over all activity folders and load each activity's
    %         already-across-subjects SPD result
    % ---------------------------------------------------------------------
    for aa = 1:nActivities
        if (options.verbose)
            fprintf("Activity: %d/%d\n", aa, nActivities);
        end

        % Retrieve the current activity name/path.
        activityName = activityFiles(aa).name;
        activityPath = fullfile(activityFiles(aa).folder, activityName);

        % Search within this activity folder for the requested projection
        % type's across-subjects SPD result.
        %
        % Example pattern:
        %   *virtuallyFoveated_SPDResultsAcrossSubjects.mat
        %
        % or:
        %   *justProjection_SPDResultsAcrossSubjects.mat
        virtually_foveated_SPD_results = dir(fullfile(activityPath, sprintf('*%s_SPDResultsAcrossSubjects.mat', options.projection_type)));

        % If no matching file exists, warn and skip this activity.
        if isempty(virtually_foveated_SPD_results)
            fprintf('No SPDResultsAcrossSubjects file found for activity %s @ path: %s\n', activityName, activityPath);
            continue;
        end

        % Build the full path to the first matching SPD results file.
        %
        % We assume there should typically be exactly one relevant file.
        thisVideo = fullfile(virtually_foveated_SPD_results(1).folder, virtually_foveated_SPD_results(1).name);

        % Load the saved .mat file.
        loadedStruct = load(thisVideo);

        % The file is expected to contain a top-level struct called
        % activityData, and inside that a single field such as:
        %   activityData.virtuallyFoveatedActivityDataAcrossSubjects
        %
        % Since we want this function to be robust to exact field naming,
        % we grab whatever field exists inside activityData.
        loadedFieldNames = fieldnames(loadedStruct.activityData);
        dataFieldName = loadedFieldNames{1};

        % Extract the actual SPD result struct for this activity.
        individual_virtually_foveated_spd_results = loadedStruct.activityData.(dataFieldName);

        % -----------------------------------------------------------------
        % Save this activity's loaded result into stack arrays
        % -----------------------------------------------------------------
        % We stack each activity along the final dimension so we can average
        % across activities later.
        virtuallyFoveatedActivityDataAcrossActivities.exponentMaps(:,:,aa) = individual_virtually_foveated_spd_results.exponentMap;
        virtuallyFoveatedActivityDataAcrossActivities.varianceMaps(:,:,aa) = individual_virtually_foveated_spd_results.varianceMap;
        virtuallyFoveatedActivityDataAcrossActivities.spdByRegions(:,:,:,aa) = individual_virtually_foveated_spd_results.spdByRegion;
        virtuallyFoveatedActivityDataAcrossActivities.medianImages(:,:,aa) = individual_virtually_foveated_spd_results.medianImage;
        virtuallyFoveatedActivityDataAcrossActivities.frameDropVector{aa} = individual_virtually_foveated_spd_results.frameDropVector;

        % Frequency vector should be shared across all activities, so we
        % simply keep overwriting it with the same value.
        virtuallyFoveatedActivityDataAcrossActivities.frq = individual_virtually_foveated_spd_results.frq;

        % -----------------------------------------------------------------
        % Optionally also load justProjection data for combined plots
        % -----------------------------------------------------------------
        if (combine_figures)
            just_projection_SPD_results = dir(fullfile(activityPath, sprintf('*%s_SPDResultsAcrossSubjects.mat', "justProjection")));

            % If the justProjection comparison file is missing, warn but do
            % not stop the function.
            if isempty(just_projection_SPD_results)
                fprintf('No justProjection SPDResultsAcrossSubjects file found for activity %s @ path: %s\n', activityName, activityPath);
            else
                % Load the justProjection file.
                thisVideo = fullfile(just_projection_SPD_results(1).folder, just_projection_SPD_results(1).name);
                loadedJustProjectionStruct = load(thisVideo);

                % As above, dynamically retrieve the inner field name.
                loadedJustProjectionFieldNames = fieldnames(loadedJustProjectionStruct.activityData);
                justProjectionFieldName = loadedJustProjectionFieldNames{1};

                % Extract the justProjection SPD data struct.
                individual_just_projection_SPD_results = loadedJustProjectionStruct.activityData.(justProjectionFieldName);

                % Stack the justProjection results across activities too.
                justProjectionActivityDataAcrossActivities.exponentMaps(:,:,aa) = individual_just_projection_SPD_results.exponentMap;
                justProjectionActivityDataAcrossActivities.varianceMaps(:,:,aa) = individual_just_projection_SPD_results.varianceMap;
                justProjectionActivityDataAcrossActivities.spdByRegions(:,:,:,aa) = individual_just_projection_SPD_results.spdByRegion;
                justProjectionActivityDataAcrossActivities.medianImages(:,:,aa) = individual_just_projection_SPD_results.medianImage;
                justProjectionActivityDataAcrossActivities.frameDropVector{aa} = individual_just_projection_SPD_results.frameDropVector;
                justProjectionActivityDataAcrossActivities.frq = individual_just_projection_SPD_results.frq;
            end
        end
    end

    % ---------------------------------------------------------------------
    % STEP 4: Average the stacked results across activities
    % ---------------------------------------------------------------------
    % exponentMaps and varianceMaps are stacked in the 3rd dimension
    % medianImages are stacked in the 3rd dimension
    % spdByRegions is stacked in the 4th dimension
    %
    % We use 'omitmissing' to avoid NaNs / missing values contaminating the
    % average if some activities had missing data.
    virtuallyFoveatedAvgExponentMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossActivities.exponentMaps(:,:,1:nActivities), 3, 'omitmissing'));
    virtuallyFoveatedAvgVarianceMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossActivities.varianceMaps(:,:,1:nActivities), 3, 'omitmissing'));
    virtuallyFoveatedAvgSpdByRegion = squeeze(mean(virtuallyFoveatedActivityDataAcrossActivities.spdByRegions(:,:,:,1:nActivities), 4, 'omitmissing'));
    virtuallyFoveatedAvgMedianImage = squeeze(mean(virtuallyFoveatedActivityDataAcrossActivities.medianImages(:,:,1:nActivities), 3, 'omitmissing'));

    % Store the averaged results back into the output struct using the same
    % field names expected elsewhere in your pipeline.
    virtuallyFoveatedActivityDataAcrossActivities.exponentMap = virtuallyFoveatedAvgExponentMap;
    virtuallyFoveatedActivityDataAcrossActivities.varianceMap = virtuallyFoveatedAvgVarianceMap;
    virtuallyFoveatedActivityDataAcrossActivities.spdByRegion = virtuallyFoveatedAvgSpdByRegion;
    virtuallyFoveatedActivityDataAcrossActivities.medianImage = virtuallyFoveatedAvgMedianImage;

    % If we also loaded justProjection data, average it in the same way.
    if (combine_figures)
        justProjectionAvgExponentMap = squeeze(mean(justProjectionActivityDataAcrossActivities.exponentMaps(:,:,1:nActivities), 3, 'omitmissing'));
        justProjectionAvgVarianceMap = squeeze(mean(justProjectionActivityDataAcrossActivities.varianceMaps(:,:,1:nActivities), 3, 'omitmissing'));
        justProjectionAvgSpdByRegion = squeeze(mean(justProjectionActivityDataAcrossActivities.spdByRegions(:,:,:,1:nActivities), 4, 'omitmissing'));
        justProjectionAvgMedianImage = squeeze(mean(justProjectionActivityDataAcrossActivities.medianImages(:,:,1:nActivities), 3, 'omitmissing'));

        justProjectionActivityDataAcrossActivities.exponentMap = justProjectionAvgExponentMap;
        justProjectionActivityDataAcrossActivities.varianceMap = justProjectionAvgVarianceMap;
        justProjectionActivityDataAcrossActivities.spdByRegion = justProjectionAvgSpdByRegion;
        justProjectionActivityDataAcrossActivities.medianImage = justProjectionAvgMedianImage;
    end

    % ---------------------------------------------------------------------
    % STEP 5: Save the combined results and optionally export figures
    % ---------------------------------------------------------------------
    if (output_dir ~= "")
        % Create the output directory if needed.
        if (~isfolder(output_dir))
            mkdir(output_dir);
        end

        % Save results in a wrapper struct named activityData so the saved
        % file structure resembles the rest of your pipeline.
        activityData.virtuallyFoveatedActivityDataAcrossActivities = virtuallyFoveatedActivityDataAcrossActivities;
        activityData.justProjectionActivityDataAcrossActivities = justProjectionActivityDataAcrossActivities;

        save(output_filepath, 'activityData');

        % If requested, generate plots summarizing the across-activities
        % result and export them as PDFs.
        if (options.save_figures)

            % plotSPDs expects the main fields like exponentMap,
            % varianceMap, spdByRegion, etc. We remove activityNames before
            % plotting because that metadata field is not part of the SPD
            % data structure expected by plotSPDs.
            plotData = virtuallyFoveatedActivityDataAcrossActivities;
            if (isfield(plotData, 'activityNames'))
                plotData = rmfield(plotData, 'activityNames');
            end

            % Do the same cleanup for the optional justProjection data.
            justProjectionPlotData = justProjectionActivityDataAcrossActivities;
            if (combine_figures)
                justProjectionPlotData = rmfield(justProjectionPlotData, 'activityNames');
            end

            % Create the summary plots.
            [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs( ...
                plotData, ...
                "fovDegrees", fovDegrees, ...
                "exponent_clim", options.exponent_clim, ...
                "variance_clim", options.variance_clim, ...
                "spd_xlim", options.spd_xlim, ...
                "spd_ylim", options.spd_ylim, ...
                "justProjectionActivityData", justProjectionPlotData ...
            );

            % Build output filenames that explicitly say AcrossActivities.
            expPath = fullfile(output_dir, sprintf('allActivities_%s_exponentMapAcrossActivities.pdf', options.projection_type));
            varPath = fullfile(output_dir, sprintf('allActivities_%s_varianceMapAcrossActivities.pdf', options.projection_type));
            spdPath = fullfile(output_dir, sprintf('allActivities_%s_spdByRegionAcrossActivities.pdf', options.projection_type));

            % Export exponent map figure.
            if (~exist(expPath, 'file') || options.overwrite_existing)
                exportgraphics(exponentMapHandle, expPath, 'ContentType', 'vector');
            end

            % Export variance map figure.
            if (~exist(varPath, 'file') || options.overwrite_existing)
                exportgraphics(varianceMapHandle, varPath, 'ContentType', 'vector');
            end

            % Export SPD-by-region figure.
            if (~exist(spdPath, 'file') || options.overwrite_existing)
                exportgraphics(spdByRegionHandle, spdPath, 'ContentType', 'vector');
            end

            % Close figures so we do not accumulate open figure windows.
            close([exponentMapHandle varianceMapHandle spdByRegionHandle]);
        end
    end
end


function mustBeMemberString(value, allowed)
% mustBeMemberString
%
% Local custom validator for MATLAB argument blocks.
%
% This function is used instead of mustbeMember so it remains compatible
% with setups where built-in argument validators may not work properly with
% strings/chars.
%
% Inputs:
%   value   - value to validate
%   allowed - array/list of allowed values
%
% Behavior:
%   Converts both the input value and allowed values to strings, then checks
%   whether the input matches one of the permitted entries.

    valueStr = string(value);
    allowedStr = string(allowed);

    if ~any(valueStr == allowedStr)
        error("Value must be one of: %s. Got: %s", strjoin(allowedStr, ", "), valueStr);
    end
end