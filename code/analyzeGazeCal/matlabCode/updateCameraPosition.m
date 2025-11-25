function [sceneGeometry,p, gazeOffset] = updateCameraPosition(sceneGeometryFile, perimeterFile, frameSet, gazeTargets, options)
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
    % More parameters
    subjectID = 'FLIC_2002';
    confidenceThreshold = 0.8;

    % Set up the paths
    dropboxBaseDir = getpref('lightLoggerAnalysis','dropboxBaseDir');
    sceneFolder = [dropboxBaseDir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/gazeCalibration/temporalFrequency/'];
    sceneGeometryFile = [sceneFolder, subjectID, '_gazeCal_SceneGeometry.mat'];
    perimeterFolder = [dropboxBaseDir, '/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/', subjectID, '/walkIndoor/temporalFrequency/'];
    perimeterFile = [perimeterFolder, subjectID, '_walkIndoor_tf_perimeter_contrast-1x5.mat'];

    % Identify the frames and gazes
    frameSet = [14253; 14763; 15153; 15586; 16013];
    gazeTargets = [
    2.9181   17.2401
   11.4868   17.9062
   -9.8265   15.3623
   -3.7883   32.2315
    1.8061    7.9323];

    % Indicate that we want to use the parpool. Can only make use of
    % nWorkers <= the total number of gaze targets
    nWorkers = 5;

    % Run the routine
    [sceneGeometry,p,gazeOffset] = updateCameraPosition(sceneGeometryFile, perimeterFile, frameSet, gazeTargets, 'nWorkers', nWorkers);

%}

arguments
    sceneGeometryFile char
    perimeterFile char
    frameSet double {mustBeVector}
    gazeTargets double {mustBeMatrix}
    options.confidenceThreshold double {mustBeScalarOrEmpty} = 0.75
    options.x0 double = [];
    options.glintFileName char = '';
    options.verbosity = 'iter'; % 'none','stage','iter';
    options.nWorkers double {mustBeNumeric} = [];
end

% Load the perimeter file
load(sceneGeometryFile,'sceneGeometry');

% Load the perimeter file
load(perimeterFile,'perimeter');

% Extract the optional arguments
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
myNewScene = @(p) updateSceneGeometry(sceneGeometry,p);
myEyePoses = @(p) estimateEyePoses(perimeter,glintData,myNewScene(p),confidenceThreshold);
myObj = @(p) calcAngleError(gazeTargets,myEyePoses(p),p);

% Define the X0 and bounds on the search for camera position and eye
% parameters.
%       x, y, z camear translation
%       azi, ele, tor camera rotation
%       eye rotation centers (common, differential)
%       cornea (axial length, k1, k2)
lb = [-10 -10 -10 -10 -10 -10];
ub = [ 10  10  10  10  10  10];

% Use passed or built-in x0
if ~isempty(options.x0)
    x0 = options.x0;
else
    x0 = [0 0 0 0 0 0];
end

% Report the search start point
if ~strcmp(options.verbosity,'none')
    fprintf('at start fval = %2.2f\n',myObj(x0));
end

    % Search
    [p,fVal] = bads(myObj,x0, lb, ub, lb, ub, [], optimset);

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
    sceneGeometry = updateSceneGeometry(sceneGeometry,p);


%% Plot
% Show the visual angle results
eyePoses = estimateEyePoses(perimeter,glintData,sceneGeometry,confidenceThreshold);
figure
gazeOffset = mean(eyePoses(:,1:2))-mean(gazeTargets(:,1:2));
plot(gazeTargets(:,1),gazeTargets(:,2),'o'); hold on
plot(eyePoses(:,1)-gazeOffset(1),eyePoses(:,2)-gazeOffset(2),'x')
axis equal
xlim([-40 40]); ylim([-40 40]);
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

function sceneGeometry = updateSceneGeometry(sceneGeometry,p)

% Camera translation and rotation
sceneGeometry.cameraPosition.translation = sceneGeometry.cameraPosition.translation + p(1:3)';
sceneGeometry.cameraPosition.rotation = sceneGeometry.cameraPosition.rotation + p(4:6);

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


