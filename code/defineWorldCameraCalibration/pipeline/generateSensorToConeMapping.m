function coneMapVar = generateSensorToConeMapping(options)
% GENERATESENSORTOCONEMAPPING Creates a 1D mapping LUT for converting
% camera RGB radiance to LMS isomerization rates based on eccentricity.
% Relies upon ISETbio tools for this.
%
% Name-Value Arguments:
%   age                 - Observer age in years (default: 56)
%   basePupilDiameterMm - Canonical pupil diameter in mm (default: 3.0)

arguments
    options.age (1,1) double = 56
    options.basePupilDiameterMm (1,1) double = 3.0
end

persistent T wlsSensor channelNames camSens

% Load the IMX219 sensitivity functions on the first pass
if isempty(camSens)
    dataFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'data',...
        'IMX219_spectralSensitivity.mat');
    load(dataFileName,'T');
    wlsSensor = T.wls;
    channelNames = {'red','green','blue'};

    % Extract the sensitivity vectors into an N-wave x 3 array
    camSens = [T.(channelNames{1}), T.(channelNames{2}), T.(channelNames{3})];
end

% Instantiate ISETbio optical and receptor components directly
lens = Lens('wave', wlsSensor);

% Apply the standard CIE 2006 age-dependent scaling factor to the lens
% density. The baseline density curve is standardized for a 32-year-old
% observer.
if options.age <= 60
    ageScalar = 1 + 0.02 * (options.age - 32);
else
    ageScalar = 1.56 + 0.0667 * (options.age - 60);
end

% In ISETbio, the Lens density property is a scalar multiplier 
% that scales the internal unitDensity spectrum.
lens.density = ageScalar;
pigment = photoPigment('wave', wlsSensor);

% Pre-calculate constants for photon conversion
h = 6.626e-34; % Planck's constant (J*s)
c = 2.998e8;   % Speed of light (m/s)
energyToQuanta = (wlsSensor .* 1e-9) ./ (h * c); % Quanta per Joule

% Eye geometry constants for calculating retinal irradiance
focalLengthM = 0.017;
pupilAreaM2 = pi * (options.basePupilDiameterMm / 2000)^2;

% Extract the cone outer segment aperture area directly
coneApertureM2 = pigment.pdArea;

% Base spectral sensitivity (Lens Transmittance * Cone Absorptance)
radiometricScalar = (pupilAreaM2 / focalLengthM^2) * coneApertureM2;

baseLMS = diag(lens.transmittance) * pigment.absorptance;
baseLMS = diag(energyToQuanta) * baseLMS .* radiometricScalar;

% Precompute the pseudo-inverse of the camera sensitivities
pinvCamSens = pinv(camSens');

% Create a 1D grid of eccentricities out to the edges of the visual field
eccGrid = 0:0.25:120;
transformTable = zeros(length(eccGrid), 9);

% Loop over scalar eccentricities to build the LUT
for ii = 1:length(eccGrid)
    ecc = eccGrid(ii);

    % Instantiate fresh object inside the loop to clear internal
    % transmittance caches
    mac = Macular('wave', wlsSensor);

    % Use the eccDensity method to obtain the properly scaled optical
    % density for this eccentricity
    mac.density = mac.eccDensity(ecc);

    % Effective LMS sensitivities at this eccentricity
    effectiveLMS = diag(mac.transmittance) * baseLMS;

    % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
    T_mat = effectiveLMS' * pinvCamSens;

    % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
    transformTable(ii, :) = T_mat(:)';
end

coneMapVar = struct();
coneMapVar.eccGrid = eccGrid;
coneMapVar.transformTable = transformTable;
coneMapVar.basePupilDiameterMm = options.basePupilDiameterMm;
coneMapVar.observerAge = options.age;

end