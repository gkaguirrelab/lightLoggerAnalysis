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
%   p                     - The parameters at the solution.
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
    observerArgs = {'sphericalAmetropia',2.25,'spectacleLens',[2.25,2,80]};

    % Combine the two argument sets
    setupArgs = [sceneArgs observerArgs];

    % Define the input variables for this particular gaze cal video
    perimeterFile = '/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Flic Experimenter/FLIC_data/lightLogger/Processing/FLIC_200X_gazeCalibration_session1_perimeter.mat';
    

fullFrameSet = [10242
       10551
       11038
       11424
       11759
       12138
       12504
       13054
       13361
       13770
       14200
       14508
       14849
       15318
       15695
       16081
       16515
       16896
       17225
       17762
       18113
       18501
       18905
       19282
       19672
       20063
       20368
       20751
       21277
       21585
       22065
       22388
       22830
       23188]';

fullGazeTargets = [-1, 1].*[ ...
            0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
            0, 15; 0, -15; -15, 0; 15, 0;...
            -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
            0, 10; 0, -7.5; -7.5, 0; 7.5, 0;...
            0, 0; -15, 15; -15, -15; 15, 15; 15, -15; ...
            0, 15; 0, -15; -15, 0; 15, 0;...
            -7.5, 7.5; -7.5, -7.5; 7.5, 7.5; 7.5, -7.5; ...
            0, 10; 0, -7.5; -7.5, 0; 7.5, 0]; % where gaze targets are
            %taken from runGazeCalibration.m (in this original file,
            %negative x = left, negative y = down. In gka model eye,
            %negative x = right, negative y = down. So we flip the sign of x).

smallSetIdx = [1, 6, 7, 8, 9];

    %frameSet = fullFrameSet(smallSetIdx);
    %gazeTargets = fullGazeTargets(smallSetIdx,:);

frameSet = fullFrameSet;
    gazeTargets = fullGazeTargets;

    % This is the x0, in case we want to pass that
    x0 = [-29.9355  -10.3699   52.3664   24.2541    2.1374   15.5410    0.9895    1.0024   17.6172   43.0417   41.0665];

    % Run the routine
    [sceneGeometry,p,gazeOffset] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs);

    % Search again starting from the prior search result
    [sceneGeometry,p] = estimateSceneGeometry(perimeterFile, frameSet, gazeTargets, 'setupArgs', setupArgs, 'x0', p);

%}

arguments
    perimeterFile char
    frameSet double {mustBeVector}
    gazeTargets double {mustBeMatrix}
    options.setupArgs cell {mustBeVector} = {}
    options.confidenceThreshold double {mustBeScalarOrEmpty} = 0.75
    options.x0 double = [];
    options.paramSearchSets = {};
    options.verbosity = 'stage'; % 'none','stage','iter';
end

% Extract the optional arguments
setupArgs = options.setupArgs;
confidenceThreshold = options.confidenceThreshold;

% Load the perimeter file
load(perimeterFile,'perimeter');

% Use just the selected frames
perimeter = perimeter.data(frameSet);

% Remove any perimeters that have no points that are above the confidence
% threshold
goodIdx = cellfun(@(x) any(x.confidence > confidenceThreshold), perimeter);
perimeter = perimeter(goodIdx);

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

% Define an objective that minimizes mismatch between targets and eye
% rotations, and uses updates in camera position
myNewScene = @(p) updateSceneGeometry(sceneGeometry,p,setupArgs);
myEyePoses = @(p) estimateEyePoses(perimeter,myNewScene(p),confidenceThreshold);
myObj = @(p) calcAngleError(gazeTargets,myEyePoses(p));

% Define the X0 and bounds on the search for camera position and eye
% parameters.
%       x, y, z camear translation
%       azi, ele, tor camera rotation
%       eye rotation centers (common, differential)
%       cornea (axial length, k1, k2)
lb = [-50 -15 40 10 -5 -10 0.8 0.8 10 40 40];
ub = [-20 -05 100 40  5  20 1.2 1.2 20 50 50];

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

    % Announce
    if ~strcmp(options.verbosity,'none')
        fprintf('fval = %2.2f\n',fVal);
    end

    % Update the sceneGeometry at the solution
    sceneGeometry = updateSceneGeometry(sceneGeometry,p,setupArgs);

    % Update the x0
    x0 = p;

end


%% Plot
% Show the visual angle results
eyePoses = estimateEyePoses(perimeter,sceneGeometry,confidenceThreshold);
figure
gazeOffset = mean(eyePoses(:,1:2))-mean(gazeTargets(:,1:2));
plot(gazeTargets(:,1),gazeTargets(:,2),'o'); hold on
plot(eyePoses(:,1)-gazeOffset(1),eyePoses(:,2)-gazeOffset(2),'x')
axis equal
xlim([-30 30]); ylim([-30 30]);
xlabel('azimuth [deg]'); ylabel('elevation [deg]');
title(sprintf('Gaze offset %2.1f, %2.1f [azi, ele]',gazeOffset));

% Render the eyes and the pupil perimeters
for ii = 1:length(perimeter)
    renderEyePose(eyePoses(ii,:),sceneGeometry);
    hold on
    conf = perimeter{ii}.confidence;
    goodIdx = conf > confidenceThreshold;
    Xp = perimeter{ii}.Xp(goodIdx);
    Yp = perimeter{ii}.Yp(goodIdx);
    plot(Xp,Yp,'*k')
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
function [eyePoses,ellipseRMSEs] = estimateEyePoses(perimeter,sceneGeometry,confidenceThreshold)

% For now, we are not using the glint
glintCoord = [];

% Define the return variables
eyePoses = nan(length(perimeter),4);
ellipseRMSEs = nan(length(perimeter),1);

% Loop through the perimeters
for ii = 1:length(perimeter)
    conf = perimeter{ii}.confidence;
    goodIdx = conf > confidenceThreshold;
    Xp = perimeter{ii}.Xp(goodIdx);
    Yp = perimeter{ii}.Yp(goodIdx);
    % Make sure we have Xp points
    if ~isempty(Xp)
        [eyePoses(ii,:),~,ellipseRMSEs(ii)] = eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry,...
            'cameraTransX0',[0;0;0],...
            'cameraTransBounds', [0;0;0]);
    end
    drawnow
end
end


% An objective function that computes error in estimation of the gaze
% target angular positions. We discount the error in the mean eye pose.
function angleError = calcAngleError(gazeTargets,eyePoses)
eyePoses = eyePoses(:,1:2);
meanEyePoses = mean(eyePoses,'omitmissing');
eyePoses = eyePoses - meanEyePoses;
meanGazeTargets = mean(gazeTargets);
gazeTargets = gazeTargets - meanGazeTargets;
idx = ~isnan(eyePoses);
angleError = norm(gazeTargets(idx)-eyePoses(idx)) + norm(meanEyePoses)/10;
end
