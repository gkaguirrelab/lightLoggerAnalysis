function activityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)

    arguments
        input_dir; 
        output_dir; 
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {};   % cell array of strings, e.g. {'walkIndoorFoveate'}
        options.verbose = false; 
        options.save_figures = false; 
        options.overwrite_existing = false; 
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
            fprintf("Activity Name: %d\%d\n", aa, n_activities); 
        end     
        activityName = activityNames{aa};

        % Initialize activity field if needed
        if ~isfield(activityDataAcrossSubjects, activityName)
            activityDataAcrossSubjects.(activityName).exponentMaps = [];
            activityDataAcrossSubjects.(activityName).interceptMaps = [];
            activityDataAcrossSubjects.(activityName).spdsByRegion = [];
            activityDataAcrossSubjects.(activityName).medianImage = [];
            activityDataAcrossSubjects.(activityName).frameDropVector = {};
        end


        % Check to see if we can skip because we do not want to overwrite existing activities
        if(output_dir ~= "")
            output_filepath = fullfile(output_dir, sprintf("%s_SPDResultsAcrossSubjects.mat", activityName)); 
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
            SPD_results = dir(fullfile(activityPath, '*SPDResults.mat'));

            if isempty(SPD_results)
                fprintf('No SPDResults file found for subject %s, activity %s\n', subjectName, activityName);
                continue;
            end

            thisVideo = fullfile(SPD_results(1).folder, SPD_results(1).name);

            % Load in the per subject + activity SPD struct 
            individual_spd_results = load(thisVideo).(activityName); 
            

            % Save this subject + activity pair of data 
            activityDataAcrossSubjects.(activityName).exponentMaps(:,:,ss) = individual_spd_results.exponentMap; 
            activityDataAcrossSubjects.(activityName).interceptMaps(:,:,ss) = individual_spd_results.interceptMap; 
            activityDataAcrossSubjects.(activityName).spdsByRegion(:,:,:,ss) = individual_spd_results.spdByRegion; 
            activityDataAcrossSubjects.(activityName).medianImage(:,:,ss) = individual_spd_results.medianImage; 
            activityDataAcrossSubjects.(activityName).frameDropVector{ss} = individual_spd_results.frameDropVector; 
            activityDataAcrossSubjects.(activityName).frq = individual_spd_results.frq;
        end

        % Average across the subjects processed so far for this activity
        avgExponentMap = squeeze(mean(activityDataAcrossSubjects.(activityName).exponentMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgVarianceMap = squeeze(mean(activityDataAcrossSubjects.(activityName).interceptMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        avgSpdByRegion = squeeze(mean(activityDataAcrossSubjects.(activityName).spdsByRegion(:,:,:,1:nSubjects), 4, 'omitmissing'));

        % Define the image resolution in degrees
        degPerPix = fovDegrees / size(avgExponentMap, 1);

        % Output the results of this if desired
        if(output_dir ~= "")
            save(output_filepath, 'activityDataAcrossSubjects');

            if(options.save_figures)
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(activityDataAcrossSubjects, "fovDegrees", options.fovDegrees); 
                
                exportgraphics(exponentMapHandle, fullfile(output_dir, sprintf('%s_exponentMapAcrossSubjects.pdf', activityName)), 'ContentType','vector');
                exportgraphics(varianceMapHandle, fullfile(output_dir, sprintf('%s_varianceMapAcrossSubjects.pdf', activityName)), 'ContentType','vector');
                exportgraphics(spdByRegionHandle, fullfile(output_dir, sprintf('%s_spdByRegionAcrossSubjects.pdf', activityName)), 'ContentType','vector');

                % Close all the figures after we saved them 
                close([exponentMapHandle varianceMapHandle spdByRegionHandle]); 

            end 

        end


    end
end