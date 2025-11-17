function [sceneGeometry,p, gazeOffset] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, options)
% Estimate eye and scene geometry for a gaze calibration measurement
%
% Syntax:
%  estimateSceneGeometry(perimeterFile, frameSet, gazeTargets)
%
% Description:
%   The appearance of the eye in a camera image is influenced by several
%   parameters, including biometric properties of the eye and the position
%   of the camera. This function estimates these parameters using a set of
%   observations of the pupil perimeter made during fixation upon known
%   targets at known visual angle locations. The analysis assumes that the
%   camera is fixed in position relative to the head.
%
%   We perform a set of searches with an ever expanding set of parameters.
%
% Inputs:
%	perimeterFile         - Char vector. Full path to a file that contains
%	                        boundary points on the perimeter of the pupil
%	                        for each frame of a video.
%   frameSet              - A 1xm vector that specifies the m frame indices
%                           (indexed from 1) which identify the set of
%                           frames from the video to guide the search.
%   gazeTargets           - A 2xm matrix that provides the positions, in
%                           degrees of visual angle, of fixation targets
%                           that correspond to each of the frames.
%
% Optional key/value pairs:
%  'setupArgs'            - Cell array. These are key-value pairs to be
%                           used when generating the sceneGeometry.
%  'confidenceThreshold'  - Numeric scalar. Only pupile perimeter points
%                           with a confidence value above this threshold
%                           will be used in the search.
%  'x0'                   - 1x11 numeric array. Can pass a starting point
%                           for the search.
%
% Outputs
%   sceneGeometry         - The resulting sceneGeometry structure.
%   p                     - 1x11 array. The parameters at the solution.
%   gazeOffset            - 1x2 array that has the offset (in degrees)
%                           between the primary position of the eye and the
%                           [0 0] coordinate of the gaze target space.
%
% Examples:
%{
    % Define some properties of the eye and of the scene that will be fixed
    % for the scene search
    sceneArgs = {
        'cameraGlintSourceRelative',[3.6;-3.2;6.9], ...
        'intrinsicCameraMatrix',[561.471804, 0.0, 200; 0.0, 562.494105, 200; 0.0, 0.0, 1.0],...
        'sensorResolution',[400 400],...
        'radialDistortionVector',[0 0]};

    % Add the args for this particular observer
    observerArgs = {'sphericalAmetropia',-1.00,'spectacleLens',[-1.00,-0.25,0]};

    % Combine the two argument sets
    setupArgs = [sceneArgs observerArgs];

    % More parameters
    subjectID = 'FLIC_2005';
    confidenceThreshold = 0.8;

    % Set up the paths
    dropboxBaseDir = getpref('lightLoggerAnalysis','dropboxBaseDir');
    saveFolders = [dropboxBaseDir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/'];
    perimeterFile = [saveFolders, subjectID, '_gazeCal_perimeter.mat'];

    % Identify the frames and gazes
    frameSet = [10349; 12362; 12783; 13153; 13518; 13962];
    gazeTargets = [
             0         0
             0   15.0000
             0  -15.0000
       15.0000         0
      -15.0000         0
        7.5000    7.5000];

    % The x0 guess
    x0 = [-28.6484   -7.3094   51.0564   24.3158    0.5042   12.1706    0.9918 0.9927   18.8754   49.3395   40.5355];

    % Indicate that we want to use the parpool. Can only make use of
    % nWorkers <= the total number of gaze targets
    nWorkers = 6;

    % Run the routine
    [sceneGeometry,p,gazeOffset] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'glintFileName',glintFileName, 'setupArgs', setupArgs, 'x0', x0, 'nWorkers', nWorkers);

%}

arguments
    perimeterFile char
    frameSet double {mustBeVector}
    gazeTargets double {mustBeMatrix}
    options.setupArgs cell {mustBeVector} = {}
    options.confidenceThreshold double {mustBeScalarOrEmpty} = 0.75
    options.x0 double = [];
    options.glintFileName char = '';
    options.paramSearchSets = {};
    options.verbosity = 'stage'; % 'none','stage','iter';
    options.nWorkers double {mustBeNumeric} = [];
end

% Extract the optional arguments
setupArgs = options.setupArgs;
confidenceThreshold = options.confidenceThreshold;
nWorkers = options.nWorkers;

% Set up the parallel pool
if isempty(nWorkers)
    nWorkers = 1;
end
switch options.verbosity
    case 'none'
        nWorkers = startParpool( nWorkers, false );
    otherwise
        nWorkers = startParpool( nWorkers, true );
end

% Load the perimeter file
load(perimeterFile,'perimeter');

% Use just the selected frames
perimeter = perimeter.data(frameSet);

% If the glintFileName is defined, load the glints
if ~isempty(options.glintFileName)
    load(options.glintFileName,'glintData')
    glintData.X = glintData.X(frameSet);
    glintData.Y = glintData.Y(frameSet);
else
    glintData.X = nan(size(frameSet));
    glintData.Y = nan(size(frameSet));
end

% Remove any perimeters that have no points that are above the confidence
% threshold
goodIdx = cellfun(@(x) any(x.confidence > confidenceThreshold), perimeter);
perimeter = perimeter(goodIdx);
glintData.X = glintData.X(goodIdx); glintData.Y = glintData.Y(goodIdx);

% Create the sceneGeometry starting point
sceneGeometry = createSceneGeometry(setupArgs{:});

% Adjust the visual angles of the gaze targets to account for magnification
% / minification produced by any spectacle lens
if isfield(sceneGeometry.refraction.retinaToCamera.magnification,'spectacle')
    gazeTargets = gazeTargets * sceneGeometry.refraction.retinaToCamera.magnification.spectacle;
end

% Define some verbosity options
switch options.verbosity
    case {'none','stage'}
        optimset.Display = 'off';
    case {'iter'}
        optimset.Display = 'iter';
    otherwise
        error('Not a valid verbosity setting');
end

% Define a nested objective function that computes error in estimation of
% the gaze target angular positions. We discount the error in the mean eye
% pose. We nest the function to support having the pBest and fValBest
% available to the main function.
pBest = [];
fValBest = Inf;

    function angleError = calcAngleError(gazeTargets,eyePoses,p)
        eyePoses = eyePoses(:,1:2);
        meanEyePoses = mean(eyePoses,'omitmissing');
        eyePoses = eyePoses - meanEyePoses;
        meanGazeTargets = mean(gazeTargets);
        gazeTargets = gazeTargets - meanGazeTargets;
        idx = ~isnan(eyePoses);
        angleError = norm(gazeTargets(idx)-eyePoses(idx)) + (norm(meanEyePoses)/10)^2;
        % Store the best solution seen
        if angleError < fValBest
            fValBest = angleError;
            pBest = p;
        end
    end


% Define an objective that minimizes mismatch between targets and eye
% rotations, and uses updates in camera position
myNewScene = @(p) updateSceneGeometry(sceneGeometry,p,setupArgs);
myEyePoses = @(p) estimateEyePoses(perimeter,glintData,myNewScene(p),confidenceThreshold);
myObj = @(p) calcAngleError(gazeTargets,myEyePoses(p),p);

% Define the X0 and bounds on the search for camera position and eye
% parameters.
%       x, y, z camear translation
%       azi, ele, tor camera rotation
%       eye rotation centers (common, differential)
%       cornea (axial length, k1, k2)
lb = [-50 -15 40 10 -5 -10 0.8 0.8 10 30 30];
ub = [-20 -05 100 40  5  20 1.2 1.2 20 60 60];

% Use passed or built-in x0
if ~isempty(options.x0)
    x0 = options.x0;
else
    x0 = [-25 -10 50 25  0  10 1.0 1.0 14 44 45];
end

% Define a progressive search strategy
if isempty(options.paramSearchSets)
    paramSearchSets = {1:3,4:6,1:6,7:8,9:11,1:11};
else
    paramSearchSets = options.paramSearchSets;
end

% Report the search start point
if ~strcmp(options.verbosity,'none')
    fprintf('at start fval = %2.2f\n',myObj(x0));
end

% Search across the paramsets
for ss = 1:length(paramSearchSets)

    % Free and lock the appropriate parameters
    paramSet = paramSearchSets{ss};
    thisLB = x0; thisUB = x0;
    thisLB(paramSet) = lb(paramSet); thisUB(paramSet) = ub(paramSet);

    % Announce
    if ~strcmp(options.verbosity,'none')
        fprintf('search %d of %d...',ss,length(paramSearchSets));
    end

    % Search
    [p,fVal] = bads(myObj,x0, thisLB, thisUB, thisLB, thisUB, [], optimset);

    % Check to make sure that the returned solution is the best
    if fVal > fValBest
        p = pBest;
        fVal = fValBest;
        % Announce
        if ~strcmp(options.verbosity,'none')
            fprintf('fval = %2.2f (using pBest)\n',fVal);
        end
    else
        % Announce
        if ~strcmp(options.verbosity,'none')
            fprintf('fval = %2.2f\n',fVal);
        end
    end

    % Update the sceneGeometry at the solution
    sceneGeometry = updateSceneGeometry(sceneGeometry,p,setupArgs);

    % Update the x0
    x0 = p;

end


%% Plot
% Show the visual angle results
eyePoses = estimateEyePoses(perimeter,glintData,sceneGeometry,confidenceThreshold);
figure
gazeOffset = mean(eyePoses(:,1:2))-mean(gazeTargets(:,1:2));
plot(gazeTargets(:,1),gazeTargets(:,2),'o'); hold on
plot(eyePoses(:,1)-gazeOffset(1),eyePoses(:,2)-gazeOffset(2),'x')
axis equal
xlim([-30 30]); ylim([-30 30]);
xlabel('azimuth [deg]'); ylabel('elevation [deg]');
title(sprintf('Gaze offset %2.1f, %2.1f [azi, ele]',gazeOffset));

% Render the eyes and the pupil perimeters
modelEyeLabelNames = {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'  'glint_01'};
modelEyePlotColors = {'.w' 'ob' '-g' '.y' 'or'};

for ii = 1:length(perimeter)
    renderEyePose(eyePoses(ii,:),sceneGeometry,...
        'modelEyeLabelNames',modelEyeLabelNames,'modelEyePlotColors',modelEyePlotColors);
    hold on
    conf = perimeter{ii}.confidence;
    goodIdx = conf > confidenceThreshold;
    Xp = perimeter{ii}.Xp(goodIdx);
    Yp = perimeter{ii}.Yp(goodIdx);
    plot(Xp,Yp,'xg')
    plot(glintData.X(ii),glintData.Y(ii),'xr')
end


end % main function


%% LOCAL FUNCTIONS

function sceneGeometry = updateSceneGeometry(sceneGeometry,p,setupArgs)

% A persistent variable so we can detect what has changed in p
persistent pLast
if isempty(pLast)
    pLast = zeros(size(p));
end

% Camera translation and rotation
sceneGeometry.cameraPosition.translation = p(1:3)';
sceneGeometry.cameraPosition.rotation = p(4:6);

% Eye rotation center
sceneGeometry.eye.meta.rotationCenterScalers = p(7:8);
sceneGeometry.eye.rotationCenters = human.rotationCenters( sceneGeometry.eye );

% Cornea properties
sceneGeometry.eye.meta.corneaAxialRadius = p(9);
sceneGeometry.eye.meta.kvals = p(10:11);
sceneGeometry.eye.cornea = human.cornea( sceneGeometry.eye );

if any(pLast(9:11) ~= p(9:11))
    % Update the glint optical system if cornea has changes
    sceneGeometry.refraction.glint.opticalSystem = ...
        assembleOpticalSystem( sceneGeometry.eye, 'surfaceSetName', 'glint', 'skipMagCalc', true, setupArgs{:});
    % Update the stopToMedium optical system with the cornea
    sceneGeometry.refraction.stopToMedium.opticalSystem = ...
        assembleOpticalSystem( sceneGeometry.eye, 'surfaceSetName', 'stopToMedium', 'skipMagCalc', true, setupArgs{:});
end

% Update pLast
pLast = p;

end


% Fit the set of frames for this scene geometry and return the fits and errors
function [eyePoses,ellipseRMSEs] = estimateEyePoses(perimeter,glintData,sceneGeometry,confidenceThreshold)

% Define the return variables
eyePoses = nan(length(perimeter),4);
ellipseRMSEs = nan(length(perimeter),1);

% Loop through the perimeters
parfor ii = 1:length(perimeter)
    conf = perimeter{ii}.confidence;
    goodIdx = conf > confidenceThreshold;
    Xp = perimeter{ii}.Xp(goodIdx);
    Yp = perimeter{ii}.Yp(goodIdx);
    glintCoord = [glintData.X(ii) glintData.Y(ii)];
    % Make sure we have Xp points
    if ~isempty(Xp)
        [eyePoses(ii,:),~,ellipseRMSEs(ii)] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry,...
            'glintTol',5,...
            'cameraTransX0',[0;0;0],...
            'cameraTransBounds', [0;0;0]);
    end
end

end


