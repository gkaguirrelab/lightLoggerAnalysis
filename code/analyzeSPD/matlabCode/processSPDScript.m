function activityData = processSPDs(options)

    arguments
        options.subjects = {};     % cell array of ints, e.g. {2001, 2003}
        options.activities = {};   % cell array of strings, e.g. {'walkIndoorFoveate'}
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
    for aa = 1:length(activityNames)
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

        % Nan out areas beyond an elliptical field of view
        ellipseTransparentParams = [240, 240, 120000, .75, 0];
        p = ellipse_ex2im(ellipse_transparent2ex(ellipseTransparentParams));
        myEllipse = @(x,y) p(1).*x.^2 + p(2).*x.*y + p(3).*y.^2 + p(4).*x + p(5).*y + p(6);
        [X, Y] = meshgrid(1:480, 1:480);
        mask = double(myEllipse(X,Y) < 1e-9);
        avgExponentMap(mask == 0) = nan;
        avgVarianceMap(mask == 0) = nan;

        % Display the maps
        figure
        imagesc(-avgExponentMap);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        idxStarts = 1:12:(480 - 24 + 1);
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Average Exponent - %s', activityName))
        axis square
        colorbar

        figure
        imagesc(avgVarianceMap, [0.015 0.035]);
        hold on
        plot(240,240,'+k')
        for ii = [5, 10, 20, 40]
            viscircles([240 240], ii/degPerPix, 'Color', 'k', 'LineWidth', 1, 'EnhanceVisibility', false);
        end
        for ii = 1:39
            plot([idxStarts(ii), idxStarts(ii)], [1 480], '-', 'Color', [0.5 0.5 0.5]);
            plot([1 480], [idxStarts(ii), idxStarts(ii)], '-', 'Color', [0.5 0.5 0.5]);
        end
        title(sprintf('Average Contrast Variance - %s', activityName))
        axis square
        colorbar

        figure
        loglog(frq, squeeze(avgSpdByRegion(20,20,:)), '-k');
        hold on
        loglog(frq, squeeze(avgSpdByRegion(31,20,:)), '-r');
        plot([10^0 10^1.5], [10^-2 10^-5], ':k')
        legend({'center','periphery'});
        ylabel('Power [contrast^2/Hz]');
        xlabel('Frequency [log Hz]');
        title(sprintf('SPDs from the center and periphery - %s', activityName));
    end
end