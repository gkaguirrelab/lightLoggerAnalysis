function activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
% Aggregate per-subject temporal SPD analysis results across activities
%
% Syntax:
%   activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir)
%   activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
%
% Description:
%   This function searches a directory of subject folders and aggregates
%   previously computed temporal SPD analysis results across subjects for
%   each activity. Subject folders are expected to follow the naming
%   convention "FLIC_<subjectID>", and each subject folder may contain one
%   or more activity subfolders. Within each activity folder, the function
%   looks for a saved "*SPDResults.mat" file that contains the per-subject
%   SPD analysis outputs.
%
%   For each activity, the function collects exponent maps, variance maps,
%   region-wise SPDs, median images, frame drop vectors, and frequency
%   support from all subjects with available data. These are stored in a
%   structured output indexed by activity name. The function can
%   optionally save the aggregated results to disk and export summary
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
%   fovDegrees            - Scalar. Field of view in degrees passed to
%                           plotSPDs when exporting summary figures.
%
% Outputs:
%   activityDataAcrossSubjects
%                         - Struct. Structure containing across-subject SPD
%                           analysis results for each activity. Each
%                           activity field may contain:
%                               .exponentMaps
%                               .varianceMaps
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
        options.overwrite_existing = false; 
        options.video_type {mustBeVideoType} = "virtuallyFoveated"; 
        options.save_figures = false; 
        options.fovDegrees = 120; 
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.color_mode {mustBeMember(options.color_mode, ["L+M+S", "L-M", "GRAY"])} = "L+M+S";
    end

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

    % Loop over subjects
    for ss = 1:nSubjects
        if options.verbose
            fprintf('Subject: %d/%d\n', ss, nSubjects);
        end

        subjectName = subjectFiles(ss).name;
        subjectPath = fullfile(subjectFiles(ss).folder, subjectName);

        % Find activity folders for this subject
        activityFiles = dir(subjectPath);
        activityFiles = activityFiles([activityFiles.isdir]);

        activityNames = {activityFiles.name};
        isRealFolder = ~ismember(activityNames, {'.', '..'});
        activityFiles = activityFiles(isRealFolder);

        availableActivityNames = {activityFiles.name};

        % If specific activities were requested, make sure they exist for this subject
        if ~isempty(options.activities)
            missingActivities = options.activities(~ismember(options.activities, availableActivityNames));
            if ~isempty(missingActivities)
                error('Requested activity "%s" not found for subject %s.', ...
                      missingActivities{1}, subjectName);
            end

            keepActivities = ismember(availableActivityNames, options.activities);
            activityFiles = activityFiles(keepActivities);
        end

        nActivities = numel(activityFiles);

        % Loop over activities for this subject
        for aa = 1:nActivities
            if options.verbose
                fprintf('Activity: %d/%d\n', aa, nActivities);
            end

            % Initialize activity data for this subject + activity 
            activityData = struct; 
            activityName = activityFiles(aa).name;
            activityPath = fullfile(activityFiles(aa).folder, activityName);

            % Check to see if we can skip because we do not want to overwrite existing activities
            if(output_dir ~= "")
                output_filepath = fullfile(output_dir, sprintf("%s_%s_%s_SPDResults.mat", subjectName, activityName, options.video_type)); 
                if(~options.overwrite_existing && isfile(output_filepath))
                    continue; 
                end
            end 

            % Find task video
            frameDropVector = [];
            task_video = dir(fullfile(activityPath, sprintf('*task_%s*.avi', options.video_type)));

            if isempty(task_video)
                fprintf('No task %s .avi file found for subject %s, activity %s\n', options.video_type, subjectName, activityName);
                continue;
            end

            thisVideo = fullfile(task_video(1).folder, task_video(1).name);

            % Call mapSPDs
            [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = mapSPDs(thisVideo, 'frameDropVector', ...
                                                                                                 frameDropVector,...
                                                                                                 'doPlot', false,...
                                                                                                 'color_mode', options.color_mode...
                                                                                                );
            activityData.(activityName).exponentMap = exponentMap; 
            activityData.(activityName).varianceMap = varianceMap; 
            activityData.(activityName).spdByRegion = spdByRegion; 
            activityData.(activityName).frq = frq;
            activityData.(activityName).medianImage = medianImage;
            activityData.(activityName).frameDropVector = frameDropVector;
                        
            % Output the results of this if desired
            if(output_dir ~= "")
                save(output_filepath, 'activityData');

                if(options.save_figures)
                    [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityData,... 
                                                                                         "fovDegrees", options.fovDegrees,...
                                                                                         "exponent_clim", options.exponent_clim,... 
                                                                                         "variance_clim", options.variance_clim, ...
                                                                                         "spd_xlim", options.spd_xlim, ...
                                                                                        "spd_ylim", options.spd_ylim ...
                                                                                    ); 
                    
                   % Build filepaths
                    expPath = fullfile(output_dir, sprintf('%s_%s_%s_exponentMap.pdf', subjectName, activityName, options.video_type));
                    varPath = fullfile(output_dir, sprintf('%s_%s_%s_varianceMap.pdf', subjectName, activityName, options.video_type));
                    spdPath = fullfile(output_dir, sprintf('%s_%s_%s_spdByRegion.pdf', subjectName, activityName, options.video_type));

                    % Exponent map
                    if (~exist(expPath, 'file') || options.overwrite_existing)
                        exportgraphics(exponentMapHandle, expPath, 'ContentType','vector');
                    end

                    % Variance map
                    if (~exist(varPath, 'file') || options.overwrite_existing)
                        exportgraphics(varianceMapHandle, varPath, 'ContentType','vector');
                    end

                    % SPD by region
                    if (~exist(spdPath, 'file') || options.overwrite_existing)
                        exportgraphics(spdByRegionHandle, spdPath, 'ContentType','vector');
                    end

                    % Close all the figures after we saved them 
                    close([exponentMapHandle varianceMapHandle spdByRegionHandle]); 
                end     
                
            end

        end

    end

end

function mustBeVideoType(value)
    % Accept both char and string
    valueStr = string(value);

    validOptions = ["justProjection", "virtuallyFoveated"];

    if ~any(valueStr == validOptions)
        error("Value must be 'justProjection' or 'virtuallyFoveated'. Got: %s", valueStr);
    end
end