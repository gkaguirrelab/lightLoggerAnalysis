function processSPDs(input_dir, output_dir, options)

    arguments
        input_dir; 
        output_dir; 
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {};   % cell array of strings, e.g. {'walkIndoorFoveate'}
        options.verbose = false;
        options.overwrite_existing = false; 
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
                output_filepath = fullfile(output_dir, sprintf("%s_%s_SPDResults.mat", subjectName, activityName)); 
                if(~options.overwrite_existing && isfile(output_filepath))
                    continue; 
                end
            end 

            % Find task video
            frameDropVector = [];
            task_video = dir(fullfile(activityPath, '*task*.avi'));

            if isempty(task_video)
                fprintf('No task .avi file found for subject %s, activity %s\n', subjectName, activityName);
                continue;
            end

            thisVideo = fullfile(task_video(1).folder, task_video(1).name);

            % Call mapSPDs
            [exponentMap, varianceMap, spdByRegion, frq, medianImage, frameDropVector] = mapSPDs(thisVideo, 'frameDropVector', frameDropVector);
            activityData.(activityName).exponentMap = exponentMap; 
            activityData.(activityName).varianceMap = varianceMap; 
            activityData.(activityName).spdByRegion = spdByRegion; 
            activityData.(activityName).frq = frq;
            activityData.(activityName).medianImage = medianImage;
            activityData.(activityName).frameDropVector = frameDropVector;
            
            % Output the results of this if desired
            if(output_dir ~= "")
                save(output_filepath, 'activityData');
            end
        end

    end

end