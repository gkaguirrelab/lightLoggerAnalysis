%% batchPlotParticipantState.m
%
% Generate participant-state figures for multiple participants and
% activities using the local timeseriesData folder structure.
%
% Expected structure:
%
%   timeseriesData/
%       FLIC_0021/
%           chat/
%               data/
%                   imu.csv
%                   3d_eye_states.csv
%                   blinks.csv
%                   gaze.csv
%           dark/
%               data/
%                   imu.csv
%                   3d_eye_states.csv
%                   blinks.csv
%                   gaze.csv
%           read/
%               data/
%                   ...
%
% Each participant's dark recording is used to normalize all of that
% participant's non-dark activity recordings.

clear;
close all;
clc;


%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------

% Folder containing participant folders such as FLIC_0021, FLIC_1029, etc.
timeseriesRoot = ...
    '/Users/sophiamirabal/Downloads/timeseriesData';

% Parent folder in which figures will be saved.
outputRoot = ...
    '/Users/sophiamirabal/Downloads/participantStateFigures';


% Participants to process.
%
% Option 1: Explicitly specify participants:
% participantIDs = { ...
%     'FLIC_0021', ...
%     'FLIC_1029', ...
%     'FLIC_1044' ...
%     };

% Option 2: Automatically detect every FLIC participant folder.
% To use this instead, comment out participantIDs above and uncomment:

participantIDs = findParticipantFolders(timeseriesRoot);


% Activities to process.
activities = { ...
    'dark', ...
    'read', ...
    'chat', ...
    'walkIndoor', ...
    'walkOutdoor' ...
    'walkBiopond', ...
    'sitBiopond', ...   
    };

% The dark recording is used as a normalization reference, but is not
% included in the activity list above by default.


%% Plotting options

saveFigures = true;
summaryOnly = true;

winSizeSec = 5;

plotApertureSeparate = true;
normalizeApertureToDark = true;
apertureBaselinePercentile = 95;

includeLuminance = false;

% Usually false for a batch run because otherwise one additional dark
% diagnostic figure is created for every activity.
plotDarkDiagnostic = false;


%% ------------------------------------------------------------------------
% Batch processing
% -------------------------------------------------------------------------

if ~isfolder(timeseriesRoot)
    error('Timeseries root folder does not exist:\n%s', timeseriesRoot);
end

if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

nParticipants = numel(participantIDs);
nActivities = numel(activities);

% Record the outcome of every attempted participant/activity combination.
results = table( ...
    'Size', [nParticipants * nActivities, 4], ...
    'VariableTypes', {'string', 'string', 'string', 'string'}, ...
    'VariableNames', {'Participant', 'Activity', 'Status', 'Message'});

resultRow = 0;

fprintf('\n');
fprintf('============================================================\n');
fprintf('Beginning participant-state batch plotting\n');
fprintf('Participants: %d\n', nParticipants);
fprintf('Activities:   %d\n', nActivities);
fprintf('============================================================\n\n');


for p = 1:nParticipants

    participantID = string(participantIDs{p});

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('Participant %d of %d: %s\n', ...
        p, nParticipants, participantID);
    fprintf('------------------------------------------------------------\n');

    participantDir = fullfile(timeseriesRoot, participantID);

    % Locate this participant's dark reference data once. The same dark
    % recording is passed to every activity for this participant.
    darkDataDir = fullfile(participantDir, 'dark', 'data');

    darkEyeStateFile = fullfile(darkDataDir, '3d_eye_states.csv');
    darkGazeFile = fullfile(darkDataDir, 'gaze.csv');
    darkBlinkFile = fullfile(darkDataDir, 'blinks.csv');

    darkFiles = [ ...
        string(darkEyeStateFile), ...
        string(darkGazeFile), ...
        string(darkBlinkFile)];

    darkFilesExist = all(isfile(darkFiles));

    if normalizeApertureToDark && ~darkFilesExist

        warning( ...
            ['Skipping all activities for %s because one or more dark ', ...
             'normalization files are missing from:\n%s'], ...
            participantID, ...
            darkDataDir);

        for a = 1:nActivities

            resultRow = resultRow + 1;

            results.Participant(resultRow) = participantID;
            results.Activity(resultRow) = string(activities{a});
            results.Status(resultRow) = "Skipped";
            results.Message(resultRow) = ...
                "Missing dark-reference files";

        end

        continue

    end


    for a = 1:nActivities

        activity = string(activities{a});
        resultRow = resultRow + 1;

        fprintf('\n[%d/%d] %s — %s\n', ...
            a, nActivities, participantID, activity);

        activityDataDir = fullfile( ...
            participantDir, ...
            activity, ...
            'data');

        imuFile = fullfile(activityDataDir, 'imu.csv');
        eyeStateFile = fullfile(activityDataDir, '3d_eye_states.csv');
        blinkFile = fullfile(activityDataDir, 'blinks.csv');
        gazeFile = fullfile(activityDataDir, 'gaze.csv');

        requiredActivityFiles = [ ...
            string(imuFile), ...
            string(eyeStateFile), ...
            string(blinkFile), ...
            string(gazeFile)];

        missingActivityFiles = ...
            requiredActivityFiles(~isfile(requiredActivityFiles));

        results.Participant(resultRow) = participantID;
        results.Activity(resultRow) = activity;

        if ~isfolder(activityDataDir)

            warning( ...
                'Activity data folder does not exist. Skipping:\n%s', ...
                activityDataDir);

            results.Status(resultRow) = "Skipped";
            results.Message(resultRow) = ...
                "Activity data folder not found";

            continue

        elseif ~isempty(missingActivityFiles)

            warning( ...
                'Missing required file(s) for %s, %s:\n%s', ...
                participantID, ...
                activity, ...
                strjoin(missingActivityFiles, newline));

            results.Status(resultRow) = "Skipped";
            results.Message(resultRow) = ...
                "One or more activity CSV files are missing";

            continue

        end


        % Use a separate output folder for every participant and activity.
        % This prevents plotParticipantState's standard filenames from
        % overwriting figures from another session.
        activityOutputDir = fullfile( ...
            outputRoot, ...
            participantID, ...
            activity);

        if ~isfolder(activityOutputDir)
            mkdir(activityOutputDir);
        end

        figureTitle = sprintf( ...
            '%s: %s', ...
            participantID, ...
            formatActivityTitleForBatch(activity));

        try

            plotParticipantState( ...
                '', ...                           % raw_dir
                '', ...                           % processing_dir
                activityOutputDir, ...
                participantID, ...
                activity, ...
                figureTitle, ...
                imuFile, ...
                eyeStateFile, ...
                blinkFile, ...
                gazeFile, ...
                'save_figures', saveFigures, ...
                'winSizeSec', winSizeSec, ...
                'plot_aperture_separate', plotApertureSeparate, ...
                'normalize_aperture_to_dark', normalizeApertureToDark, ...
                'dark_eyeStateData', darkEyeStateFile, ...
                'dark_gazeData', darkGazeFile, ...
                'dark_blinkData', darkBlinkFile, ...
                'aperture_baseline_percentile', ...
                    apertureBaselinePercentile, ...
                'plot_dark_diagnostic', plotDarkDiagnostic, ...
                'include_luminance', includeLuminance, ...
                'summary_only', summaryOnly);
                

            results.Status(resultRow) = "Completed";
            results.Message(resultRow) = "";

            fprintf('Completed: %s — %s\n', participantID, activity);

            % If figures are not being saved, leave them open for
            % inspection. If they are being saved, the main function
            % already closes its exported figures.
            if saveFigures
                close all force
            end

        catch ME

            warning( ...
                'Could not process %s — %s:\n%s', ...
                participantID, ...
                activity, ...
                ME.message);

            results.Status(resultRow) = "Failed";
            results.Message(resultRow) = string(ME.message);

            % Prevent figures left behind by a failed iteration from
            % accumulating during the remainder of the batch.
            close all force

        end

    end

end


%% ------------------------------------------------------------------------
% Save and display batch summary
% -------------------------------------------------------------------------

summaryFile = fullfile( ...
    outputRoot, ...
    'participantState_batch_summary.csv');

writetable(results, summaryFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('Batch processing complete\n');
fprintf('Completed: %d\n', sum(results.Status == "Completed"));
fprintf('Skipped:   %d\n', sum(results.Status == "Skipped"));
fprintf('Failed:    %d\n', sum(results.Status == "Failed"));
fprintf('Summary:   %s\n', summaryFile);
fprintf('============================================================\n\n');

disp(results);


%% ------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function participantIDs = findParticipantFolders(timeseriesRoot)
% Automatically find folders whose names begin with FLIC_.

    directoryContents = dir(fullfile(timeseriesRoot, 'FLIC_*'));

    isParticipantFolder = [directoryContents.isdir];

    directoryContents = ...
        directoryContents(isParticipantFolder);

    participantIDs = {directoryContents.name};

    participantIDs = sort(participantIDs);

end


function formattedTitle = formatActivityTitleForBatch(activity)
% Convert camelCase activity names into readable display titles.

    activity = char(activity);

    % Add a space between a lowercase letter and a following uppercase
    % letter, for example walkOutdoor -> walk Outdoor.
    formattedTitle = regexprep( ...
        activity, ...
        '([a-z])([A-Z])', ...
        '$1 $2');

    formattedTitle = strtrim(formattedTitle);

    if ~isempty(formattedTitle)
        formattedTitle(1) = upper(formattedTitle(1));
    end

end