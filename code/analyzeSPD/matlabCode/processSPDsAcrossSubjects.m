function activityData = processSPDsAcrossSubjects(options)

    arguments
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {};   % cell array of strings, e.g. {'walkIndoorFoveate'}
        options.verbose = false; 
        options.save_at = ""; 
    end

    % First, let's get the location of where the virtually foveated videos live
    FLIC_processing_dir = getpref("lightLoggerAnalysis", "FLIC_processing_dir");

    % Initialize output
    activityData = struct();

    % Find all subject folders matching ^FLIC_\d+$
    subjectFiles = dir(fullfile(FLIC_processing_dir, 'FLIC_*'));
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
        activityData.subjectIDs{ss} = subjectFiles(ss).name;
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
        if(verbose)
            fprintf("Activity Name: %d\%d\n", aa, n_activities); 
        end     
        activityName = activityNames{aa};

        % Initialize activity field if needed
        if ~isfield(activityData, activityName)
            activityData.(activityName).exponentMaps = [];
            activityData.(activityName).interceptMaps = [];
            activityData.(activityName).spdsByRegion = [];
            activityData.(activityName).medianImage = [];
            activityData.(activityName).frameDropVector = {};
        end

        % Loop over subjects
        for ss = 1:nSubjects
            if(vebrose)
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

            % Find task video
            frameDropVector = [];
            task_video = dir(fullfile(activityPath, '*task*.avi'));

            if isempty(task_video)
                fprintf('No task .avi file found for subject %s, activity %s\n', subjectName, activityName);
                continue;
            end

            thisVideo = fullfile(task_video(1).folder, task_video(1).name);

            % Call mapSPDs
            [activityData.(activityName).exponentMaps(:,:,ss), ...
             activityData.(activityName).interceptMaps(:,:,ss), ...
             activityData.(activityName).spdsByRegion(:,:,:,ss), ...
             frq, ...
             activityData.(activityName).medianImage(:,:,ss), ...
             activityData.(activityName).frameDropVector{ss}] = ...
                mapSPDs(thisVideo, 'frameDropVector', frameDropVector);

            activityData.(activityName).frq = frq;
        end

        % Average across the subjects processed so far for this activity
        avgExponentMap = squeeze(mean(activityData.(activityName).exponentMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgVarianceMap = squeeze(mean(activityData.(activityName).interceptMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgSpdByRegion = squeeze(mean(activityData.(activityName).spdsByRegion(:,:,:,1:nSubjects), 4, 'omitmissing'));

        % Define the image resolution in degrees
        fovDegrees = 120;
        degPerPix = fovDegrees / size(avgExponentMap, 1);

        activityData.(activityName).fovDegrees = fovDegrees; 
        activityData.(activityName).degPerPix = degPerPix; 

        % 

    end
end