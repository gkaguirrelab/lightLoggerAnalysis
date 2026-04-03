function virtuallyFoveatedActivityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
% Aggregate per-subject temporal SPD analysis results across activities
%
% Syntax:
%   virtuallyFoveatedActivityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir)
%   virtuallyFoveatedActivityDataAcrossSubjects = processSPDsAcrossSubjects(input_dir, output_dir, options)
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
%   virtuallyFoveatedActivityDataAcrossSubjects
%                         - Struct. Structure containing across-subject SPD
%                           analysis results for each activity. Each
%                           activity field may contain:
%                               .exponentMaps
%                               .interceptMaps
%                               .spdByRegion
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

    virtuallyFoveatedActivityDataAcrossSubjects = processSPDsAcrossSubjects( ...
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
        options.exponent_clim = false
        options.variance_clim = false
        options.spd_xlim = false
        options.spd_ylim = false
        options.combine_figures = false; 
        options.n_participants = 1; 
        options.across_subject_deviation = 0; 
    end
    
    fovDegrees = options.fovDegrees; 
    combine_figures = options.combine_figures; 

    % Initialize output
    virtuallyFoveatedActivityDataAcrossSubjects = struct();
    justProjectionActivityDataAcrossSubjects = false; 
    if(combine_figures)
        justProjectionActivityDataAcrossSubjects = struct(); 
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

    % Save subjectIDs
    for ss = 1:nSubjects
        virtuallyFoveatedActivityDataAcrossSubjects.subjectIDs{ss} = subjectFiles(ss).name;
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
        if ~isfield(virtuallyFoveatedActivityDataAcrossSubjects, activityName)
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).exponentMaps = [];
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).varianceMaps = [];
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).spdByRegion = [];
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).medianImage = [];
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).frameDropVector = {};
        end

        if(combine_figures && ~isfield(justProjectionActivityDataAcrossSubjects, activityName))
            justProjectionActivityDataAcrossSubjects.(activityName).exponentMaps = [];
            justProjectionActivityDataAcrossSubjects.(activityName).varianceMaps = [];
            justProjectionActivityDataAcrossSubjects.(activityName).spdByRegion = [];
            justProjectionActivityDataAcrossSubjects.(activityName).medianImage = [];
            justProjectionActivityDataAcrossSubjects.(activityName).frameDropVector = {};
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
            virtually_foveated_SPD_results = dir(fullfile(activityPath, sprintf('*%s_SPDResults.mat', options.projection_type)));

            if isempty(virtually_foveated_SPD_results)
                fprintf('No SPDResults file found for subject %s, activity %s @ path: %s\n', subjectName, activityName, activityPath);
                continue;
            end

            thisVideo = fullfile(virtually_foveated_SPD_results(1).folder, virtually_foveated_SPD_results(1).name);

            % Load in the per subject + activity SPD struct 
            individual_virtually_foveated_spd_results = load(thisVideo).activityData.(activityName); 
            
            % Save this subject + activity pair of data 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).exponentMaps(:,:,ss) = individual_virtually_foveated_spd_results.exponentMap; 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).varianceMaps(:,:,ss) = individual_virtually_foveated_spd_results.varianceMap; 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).spdByRegions(:,:,:,ss) = individual_virtually_foveated_spd_results.spdByRegion; 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).medianImages(:,:,ss) = individual_virtually_foveated_spd_results.medianImage; 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).frameDropVector{ss} = individual_virtually_foveated_spd_results.frameDropVector; 
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).frq = individual_virtually_foveated_spd_results.frq;

            if(combine_figures)
                % Load in the per subject + activity SPD struct for the non foveated 
                just_projection_SPD_results = dir(fullfile(activityPath, sprintf('*%s_SPDResults.mat', "justProjection")));
                thisVideo = fullfile(just_projection_SPD_results(1).folder, just_projection_SPD_results(1).name);
                individual_just_projection_SPD_results = load(thisVideo).activityData.(activityName); 

                justProjectionActivityDataAcrossSubjects.(activityName).exponentMaps(:,:,ss) = individual_just_projection_SPD_results.exponentMap; 
                justProjectionActivityDataAcrossSubjects.(activityName).varianceMaps(:,:,ss) = individual_just_projection_SPD_results.varianceMap; 
                justProjectionActivityDataAcrossSubjects.(activityName).spdByRegions(:,:,:,ss) = individual_just_projection_SPD_results.spdByRegion; 
                justProjectionActivityDataAcrossSubjects.(activityName).medianImages(:,:,ss) = individual_just_projection_SPD_results.medianImage; 
                justProjectionActivityDataAcrossSubjects.(activityName).frameDropVector{ss} = individual_just_projection_SPD_results.frameDropVector; 
                justProjectionActivityDataAcrossSubjects.(activityName).frq = individual_just_projection_SPD_results.frq;

            end 

        end

        % Average across the subjects processed so far for this activity
        virtuallyFoveatedAvgExponentMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.(activityName).exponentMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        virtuallyFoveatedAvgVarianceMap = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.(activityName).varianceMaps(:,:,1:nSubjects), 3, 'omitmissing'));
        virtuallyFoveatedAvgSpdByRegion = squeeze(mean(virtuallyFoveatedActivityDataAcrossSubjects.(activityName).spdByRegions(:,:,:,1:nSubjects), 4, 'omitmissing'));

        % Store the averages
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).exponentMap = virtuallyFoveatedAvgExponentMap; 
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).varianceMap = virtuallyFoveatedAvgVarianceMap;
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).spdByRegion = virtuallyFoveatedAvgSpdByRegion;  

        % Store the averages for the just projection struct if desired 
        if(combine_figures)
            justProjectionAvgExponentMap = squeeze(mean(  justProjectionActivityDataAcrossSubjects.(activityName).exponentMaps(:,:,1:nSubjects), 3, 'omitmissing'));
            justProjectionAvgVarianceMap = squeeze(mean(  justProjectionActivityDataAcrossSubjects.(activityName).varianceMaps(:,:,1:nSubjects), 3, 'omitmissing'));
            justProjectionAvgSpdByRegion = squeeze(mean(  justProjectionActivityDataAcrossSubjects.(activityName).spdByRegions(:,:,:,1:nSubjects), 4, 'omitmissing'));

            justProjectionActivityDataAcrossSubjects.(activityName).exponentMap = justProjectionAvgExponentMap; 
            justProjectionActivityDataAcrossSubjects.(activityName).varianceMap = justProjectionAvgVarianceMap;
            justProjectionActivityDataAcrossSubjects.(activityName).spdByRegion = justProjectionAvgSpdByRegion;  
        end 

        % Also, importantly, we need to compute the standard deviation of the regionAverages means
        % ACROSS subjects 
        for ii = 1:nSubjects
            
            [~. ~. ~, subject_region_averages] = plotSPDs("return_region_averages_only", true)
            virtuallyFoveatedActivityDataAcrossSubjects.(activityName).regionAveragesAcrossSubjects{ii} = regionAverages.virtuallyFoveated; 

            if(combine_figures)
                justProjectionActivityDataAcrossSubjects.(activityName).regionAveragesAcrossSubjects{ii} = regionAverages.justProjection; 
            end 
        end 

        % Importantly, save the subjects used to calculate this 
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).subjects = options.subjects; 
        
        % Then calculate standard deviation
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).acrossSubjectSTD.center = std(regionAverages.virtuallyFoveated.center);
        virtuallyFoveatedActivityDataAcrossSubjects.(activityName).acrossSubjectSTD.periphery = std(regionAverages.virtuallyFoveated.periphery); 
        if(combine_figures)
            justProjectionActivityDataAcrossSubjects.(activityName).subjects = options.subjects; 

            justProjectionActivityDataAcrossSubjects.(activityName).acrossSubjectSTD.center = std(regionAverages.justProjection.center);
            justProjectionActivityDataAcrossSubjects.(activityName).acrossSubjectSTD.periphery = std(regionAverages.justProjection.center); 

        end 


        % Output the results of this if desired
        if(output_dir ~= "")
            % Remove this field we no longer need to keep track of subjects
            virtuallyFoveatedActivityDataAcrossSubjects = rmfield(virtuallyFoveatedActivityDataAcrossSubjects, 'subjectIDs'); 
            
            % Save the activity data
            activityData = virtuallyFoveatedActivityDataAcrossSubjects; 
            save(output_filepath, 'activityData');

            if(options.save_figures)
                assert(isstruct(virtuallyFoveatedActivityDataAcrossSubjects)); 
                         
                [exponentMapHandle, varianceMapHandle, spdByRegionHandle] = plotSPDs(virtuallyFoveatedActivityDataAcrossSubjects, ...
                                                                                     "fovDegrees", options.fovDegrees,...
                                                                                     "exponent_clim", options.exponent_clim,... 
                                                                                     "variance_clim", options.variance_clim, ...
                                                                                     "spd_xlim", options.spd_xlim, ...
                                                                                     "spd_ylim", options.spd_ylim, ...
                                                                                     "justProjectionActivityData", justProjectionActivityDataAcrossSubjects,...
                                                                                     "num_participants", options.n_participants...
                                                                                     ); 
                

                % Build filepaths
                expPath = fullfile(output_dir, sprintf('%s_%s_exponentMapAcrossSubjects.pdf', activityName, options.projection_type));
                varPath = fullfile(output_dir, sprintf('%s_%s_varianceMapAcrossSubjects.pdf', activityName, options.projection_type));
                spdPath = fullfile(output_dir, sprintf('%s_%s_spdByRegionAcrossSubjects.pdf', activityName, options.projection_type));

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