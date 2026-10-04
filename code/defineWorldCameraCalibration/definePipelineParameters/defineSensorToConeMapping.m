% defineSensorToConeMapping.m

% This script uses ISETbio to define a variable that supports the
% conversion of integrated radiance values from the camera sensors to cone
% isomerization rates for an observer. The current implementation assumes
% the lens density typical for a 25 year-old observer.

% Hard code some options
options.age  = 25;
options.basePupilDiameterMm = 3.0;

% Load the IMX219 sensitivity functions on the first pass
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');
wlsSensor = T.wls;
channelNames = {'red','green','blue'};

% Extract the sensitivity vectors into an N-wave x 3 array
camSens = [T.(channelNames{1}), T.(channelNames{2}), T.(channelNames{3})];

% Instantiate ISETbio lens directly
lens = Lens('wave', wlsSensor);

% Apply the standard CIE 2006 age-dependent scaling factor to the lens density.
if options.age <= 60
    ageScalar = 1 + 0.02 * (options.age - 32);
else
    ageScalar = 1.56 + 0.0667 * (options.age - 60);
end
lens.density = ageScalar;

% Pre-calculate constants for photon conversion
h = 6.626e-34; % Planck's constant (J*s)
c = 2.998e8;   % Speed of light (m/s)
energyToQuanta = (wlsSensor .* 1e-9) ./ (h * c); % Quanta per Joule

% Eye geometry constants for calculating retinal irradiance
focalLengthM = 0.017;
pupilAreaM2 = pi * (options.basePupilDiameterMm / 2000)^2;

% Precompute the pseudo-inverse of the camera sensitivities
pinvCamSens = pinv(camSens');

% Create a 1D grid of eccentricities out to the edges of the visual field
eccGrid = 0:0.25:120;
transformTable = zeros(length(eccGrid), 9);
densityTable = zeros(length(eccGrid), 3); % NEW: Store L, M, S cone densities

% Loop over scalar eccentricities to build the LUT
for ii = 1:length(eccGrid)
    ecc = eccGrid(ii);

    % Clamp eccentricity safely within the 29-degree radius of the 58-deg FOV lattice
    safeEcc = min(ecc, 28);

    % Instantiate minimal cMosaic to retrieve macular pigment and outer segment properties
    % Evaluate on the VERTICAL meridian [0 safeEcc] instead of horizontal 
    % to bypass the optic disk (blind spot) located at 12-18 degrees horizontally.
    suppressedText = evalc('cm = cMosaic(''eccentricityDegs'', [0 safeEcc], ''sizeDegs'', [0.2 0.2], ''wave'', wlsSensor);');

    % Clamp eccentricity for the anatomical query to 60 degrees (approx 18mm) 
    % to stay safely within the bounds of the empirical dataset.
    safeEccAperture = min(ecc, 60);
    eccMeters = safeEccAperture * (300 * 1e-6);
    
    % Dynamically calculate the scaled cone outer segment aperture area.
    [~, apertureDiameterMeters] = coneSizeReadData('eccentricity', eccMeters, 'angle', 0);
    coneApertureM2 = pi * (apertureDiameterMeters / 2)^2;
    
    % Extract total cone density ---
    coneDensitySqMm = coneDensityReadData('eccentricity', eccMeters, 'angle', 0);

    % Convert cones/mm^2 to cones/degree^2 (assuming ~300 microns/degree)
    mmPerDeg = 300 * 1e-3;
    coneDensitySqDeg = coneDensitySqMm * (mmPerDeg^2);

    % Assume a fixed, 10% S-cone fraction
    sFraction = 0.1;
    
    % Distribute the remaining fraction to L and M cones (using standard 2:1 ratio)
    lmFraction = 1.0 - sFraction;
    lFraction = lmFraction * (0.6 / 0.9);
    mFraction = lmFraction * (0.3 / 0.9);

    % Store the densities
    densityTable(ii, 1) = coneDensitySqDeg * lFraction; % L-cone
    densityTable(ii, 2) = coneDensitySqDeg * mFraction; % M-cone
    densityTable(ii, 3) = coneDensitySqDeg * sFraction; % S-cone

    % Recalculate radiometric scalar for this eccentricity
    radiometricScalar = (pupilAreaM2 / focalLengthM^2) * coneApertureM2;

    % Base LMS absorptance (inherently includes eccentricity-scaled optical density)
    baseLMS = diag(lens.transmittance) * cm.pigment.absorptance;
    baseLMS = diag(energyToQuanta) * baseLMS .* radiometricScalar;

    % Effective LMS sensitivities at this eccentricity 
    effectiveLMS = diag(cm.macular.transmittance) * baseLMS;

    % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
    T_mat = effectiveLMS' * pinvCamSens;

    % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
    transformTable(ii, :) = T_mat(:)';

    assert(~any(isnan(T_mat(:)')))

end

coneMapVar = struct();
coneMapVar.eccGrid = eccGrid;
coneMapVar.transformTable = transformTable;
coneMapVar.densityTable = densityTable;
coneMapVar.basePupilDiameterMm = options.basePupilDiameterMm;
coneMapVar.observerAge = options.age;

% Save the coneMapVar in the "derived" directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radianceToConeRateSupport.mat');
readme = ['Created by defineSensorToConeMapping.\n'...
    'coneMapVar -- a structure with values needed for conversion of radiance to cone isomerization rate.\n'];
save(saveFileName,'readme','coneMapVar');