function activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
% Aggregate per-subject temporal SPD analysis results across activities
%
% Syntax:
%   activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir)
%   activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
%
% Description:
%   This function searches a directory of subject folders and aggregates
%   previously computed SPD analysis results across subjects for each
%   activity. Subject folders are expected to follow the naming convention
%   "FLIC_<subjectID>", and each subject folder may contain one or more
%   activity subfolders. Within each activity folder, the function looks
%   for a saved "*SPDResults.mat" file that contains the per-subject SPD
%   analysis outputs.
%
%   For each activity, the function collects the exponent maps, variance
%   (or intercept) maps, regional SPDs, median images, frame drop vectors,
%   and frequency support from all subjects with available data. These are
%   stored in a structured output indexed by activity name. The function
%   can optionally save the aggregated results to disk and export summary
%   figures showing the across-subject exponent map, variance map, and
%   regional SPD summaries.
%
% Inputs:
%   input_dir             - Char/string. Path to the directory containing
%                           subject folders named like "FLIC_2001",
%                           "FLIC_2002", etc.
%   output_dir            - Char/string. Path to the directory where
%                           aggregated results and optional figures should
%                           be saved. If empty, results are returned but
%                           not written to disk.
%
% Optional key/value pairs:
%   subjects              - Cell array. List of subject numeric IDs to
%                           include, e.g. {2001, 2003}. If empty, all
%                           valid subject folders in input_dir are used.
%   activities            - Cell array. List of activity folder names to
%                           process, e.g. {'walkIndoorFoveate'}. If empty,
%                           the function determines the union of all
%                           activity folders found across subjects.
%   verbose               - Logical. If true, print progress information
%                           while processing activities and subjects.
%   save_figures          - Logical. If true, generate and export summary
%                           figures for each activity to output_dir.
%   overwrite_existing    - Logical. If false, skip processing an activity
%                           when a saved across-subject results file for
%                           that activity already exists in output_dir.
%   fovDegrees            - Scalar. Horizontal or vertical field of view in
%                           degrees used when plotting the aggregated SPD
%                           maps.
%
% Outputs:
%   activityDataAcrossSubjects
%                         - Struct. Structure containing across-subject SPD
%                           analysis results for each activity. Each
%                           activity field may contain:
%                               .exponentMaps
%                               .interceptMaps
%                               .spdsByRegion
%                               .medianImage
%                               .frameDropVector
%                               .frq
%                           The structure also contains:
%                               .subjectIDs
%
% Examples:
%{
    input_dir = "/path/to/FLIC_processing";
    output_dir = "/path/to/SPD_across_subjects";

    activityDataAcrossSubjects = processSPDsAcrossSubjects( ...
        input_dir, ...
        output_dir, ...
        "subjects", {2001, 2003, 2005}, ...
        "activities", {'walkIndoorFoveate', 'taskOutdoor'}, ...
        "verbose", true, ...
        "save_figures", true ...
    );
%}
    arguments
        input_dir; 
        output_dir; 
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {};   % cell array of strings, e.g. {'walkIndoorFoveate'}
        options.verbose = false; 
        options.save_figures = false; 
        options.overwrite_existing = false; 
        options.projection_type {mustBeMember(options.projection_type, ["justProjection", "virtuallyFoveated"])} = "virtuallyFoveated";
        options.fovDegrees = 120; 

    end

    fovDegrees = options.fovDegrees; 

    % Initialize output
    activityDataAcrossSubjects = struct();

    % Find all subject folders matching ^FLIC_\d+$
    subjectFiles = dir(fullfile(input_dir, 'FLIC_*'));
    subjectFiles = subjectFiles([subjectFiles.isdir]);

    subjectNames = {subjectFiles.name};
    isValidSubject = ~cellfun('isempty', regexp(subjectNames, '^FLIC_\d+$', 'once'));
    subjectFiles = subjectFiles(isValidSubject);

    % If specific subjects were passed, filter to only those
    if ~isempty(options.subjects)
        requestedSubjectNames = cellfun(@(x) sprintf('FLIC_%d', x), ...
                                        options.subjects, ...
                                        'UniformOutput', false);
        keepSubjects = ismember({subjectFiles.name}, requestedSubjectNames);
        subjectFiles = subjectFiles(keepSubjects);
    end

    nSubjects = length(subjectFiles);

    % Save subjectIDs
    for ss = 1:nSubjects
        activityDataAcrossSubjects.subjectIDs{ss} = subjectFiles(ss).name;
    end

    % Determine activity list
    if ~isempty(options.activities)
        activityNames = options.activities;
    else
        activityNames = {};
        for ss = 1:nSubjects
            subjectPath = fullfile(subjectFiles(ss).folder, subjectFiles(ss).name);

            activityFiles = dir(subjectPath);
            activityFiles = activityFiles([activityFiles.isdir]);

            thisActivityNames = {activityFiles.name};
            isRealFolder = ~ismember(thisActivityNames, {'.', '..'});
            thisActivityNames = thisActivityNames(isRealFolder);

            activityNames = union(activityNames, thisActivityNames);
        end
    end

    % Loop over activities
    n_activities = numel(activityNames); 
    for aa = 1:length(activityNames)
        if(options.verbose)
            fprintf("Activity Name: %d/%d\n", aa, n_activities); 
        end     
        activityName = activityNames{aa};

        % Initialize activity field if needed
        if ~isfield(activityDataAcrossSubjects, activityName)
            activityDataAcrossSubjects.(activityName).exponentMaps = [];
            activityDataAcrossSubjects.(activityName).varianceMaps = [];
            activityDataAcrossSubjects.(activityName).spdsByRegion = [];
            activityDataAcrossSubjects.(activityName).medianImage = [];
            activityDataAcrossSubjects.(activityName).frameDropVector = {};
        end


        % Check to see if we can skip because we do not want to overwrite existing activities
        if(output_dir ~= "")
            output_filepath = fullfile(output_dir, sprintf("%s_%s_SPDResultsAcrossSubjects.mat", activityName, options.projection_type)); 
            if(~options.overwrite_existing && isfile(output_filepath))
                continue; 
            end
        end 

        % Loop over subjects
        for ss = 1:nSubjects
            if(options.verbose)
                fprintf("Subject: %d/%d\n", ss, nSubjects);
            end 

            subjectName = subjectFiles(ss).name;
            subjectPath = fullfile(subjectFiles(ss).folder, subjectName);
            activityPath = fullfile(subjectPath, activityName);

            % If activities were explicitly requested, error when missing
            if ~isfolder(activityPath)
                if ~isempty(options.activities)
                    error('Requested activity "%s" not found for subject %s.', ...
                          activityName, subjectName);
                else
                    fprintf('Activity %s not found for subject %s\n', activityName, subjectName);
                    continue;
                end
            end

            % Find the activityData struct for just this subject + activity pair
            frameDropVector = [];
            SPD_results = dir(fullfile(activityPath, sprintf('*%s_SPDResults.mat', options.projection_type)));

            if isempty(SPD_results)
                fprintf('No SPDResults file found for subject %s, activity %s @ path: %s\n', subjectName, activityName, activityPath);
                continue;
            end

            thisVideo = fullfile(SPD_results(1).folder, SPD_results(1).name);

            % Load in the per subject + activity SPD struct 
            individual_spd_results = load(thisVideo).(activityName); 
            

            % Save this subject + activity pair of data 
            activityDataAcrossSubjects.(activityName).exponentMaps(:,:,ss) = individual_spd_results.exponentMap; 
            activityDataAcrossSubjects.(activityName).varianceMaps(:,:,ss) = individual_spd_results.varianceMap; 
            activityDataAcrossSubjects.(activityName).spdsByRegion(:,:,:,ss) = individual_spd_results.spdByRegion; 
            activityDataAcrossSubjects.(activityName).medianImage(:,:,ss) = individual_spd_results.medianImage; 
            activityDataAcrossSubjects.(activityName).frameDropVector{ss} = individual_spd_results.frameDropVector; 
            activityDataAcrossSubjects.(activityName).frq = individual_spd_results.frq;
        end

        % Average across the subjects processed so far for this activity
        avgExponentMap = squeeze(mean(activityDataAcrossSubjects.(activityName).exponentMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgVarianceMap = squeeze(mean(activityDataAcrossSubjects.(activityName).varianceMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgSpdByRegion = squeeze(mean(activityDataAcrossSubjects.(activityName).spdsByRegion(:,:,:,1:nSubjects), 4, 'omitmissing'));

        % Define the image resolution in degrees
        degPerPix = fovDegrees / size(avgExponentMap, 1);

        % Output the results of this if desired
        if(output_dir ~= "")
            activityData = activityDataAcrossSubjects; 
            save(output_filepath, 'activityData');

            if(options.save_figures)
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityDataAcrossSubjects, "fovDegrees", options.fovDegrees); 
                
                exportgraphics(exponentMapHandle, fullfile(output_dir, sprintf('%s_%s_exponentMapAcrossSubjects.pdf', activityName, options.projection_type)), 'ContentType','vector');
                exportgraphics(varianceMapHandle, fullfile(output_dir, sprintf('%s_%s_varianceMapAcrossSubjects.pdf', activityName, options.projection_type)), 'ContentType','vector');
                exportgraphics(spdByRegionHandle, fullfile(output_dir, sprintf('%s_%s_spdByRegionAcrossSubjects.pdf', activityName, options.projection_type)), 'ContentType','vector');

                % Close all the figures after we saved them 
                close([exponentMapHandle varianceMapHandle spdByRegionHandle]); 

            end 

        end


    end
end