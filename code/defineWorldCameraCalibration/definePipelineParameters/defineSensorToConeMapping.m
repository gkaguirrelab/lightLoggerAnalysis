% defineSensorToConeMapping.m

% This script uses ISETbio to define a variable that supports the
% conversion of integrated radiance values from the camera sensors to cone
% isomerization rates for an observer. We loop over a range of ages (and
% thus lens density parameters) and pupil diameters to precompute optics.

% Housekeeping
clear

% Hard code some options
options.ageRange  = [18,50];
options.basePupilDiameterMm = 3.0;
options.pupilDiameterMmRange = 1.5:0.5:8.0; % Range for pre-computing PSFs

% Load the IMX219 sensitivity functions
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

% Instantiate continuous Macular and PhotoPigment objects.
% Bypassing discrete cMosaic spatial generation eliminates lattice noise
% and drastically improves execution speed.
suppressedText = evalc('dummyCm = cMosaic(''sizeDegs'', [0.1 0.1], ''wave'', wlsSensor);');
macularObj = dummyCm.macular;
pigmentObj = dummyCm.pigment;

% Cache the baseline foveal densities to scale inside the loop
baseMacularDensity = macularObj.density;
basePigmentOpticalDensity = pigmentObj.opticalDensity;

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
densityTable = zeros(length(eccGrid), 3); % Store L, M, S cone densities

% Preallocate a struct array with empty fields for the cone mapping
coneMapVar(options.ageRange(2)) = struct(...
    'eccGrid', [], 'transformTable', [],...
    'densityTable', [], 'basePupilDiameterMm', [],...
    'observerAge', []);

% Report we are starting age loop
fprintf('Generating %d ISETbio age profiles', range(options.ageRange))

% Loop over observer ages
for age = options.ageRange(1):options.ageRange(2)
    fprintf('.');

    % Apply the standard CIE 2006 age-dependent scaling factor to the lens density
    if age <= 60
        ageScalar = 1 + 0.02 * (age - 32);
    else
        ageScalar = 1.56 + 0.0667 * (age - 60);
    end
    lens.density = ageScalar;

    % Loop over scalar eccentricities to build the LUT
    for ii = 1:length(eccGrid)
        ecc = eccGrid(ii);

        % Clamp eccentricity for the anatomical query to 60 degrees
        safeEccAperture = min(ecc, 60);
        eccMeters = safeEccAperture * (300 * 1e-6);

        % Calculate the scaled cone outer segment aperture area
        [~, apertureDiameterMeters] = coneSizeReadData('eccentricity', eccMeters, 'angle', 0);
        coneApertureM2 = pi * (apertureDiameterMeters / 2)^2;

        % Extract total cone density
        coneDensitySqMm = coneDensityReadData('eccentricity', eccMeters, 'angle', 0);
        mmPerDeg = 300 * 1e-3;
        coneDensitySqDeg = coneDensitySqMm * (mmPerDeg^2);

        % Assume a fixed, 10% S-cone fraction
        sFraction = 0.1;
        lmFraction = 1.0 - sFraction;
        lFraction = lmFraction * (0.6 / 0.9);
        mFraction = lmFraction * (0.3 / 0.9);

        % Store the densities
        densityTable(ii, 1) = coneDensitySqDeg * lFraction;
        densityTable(ii, 2) = coneDensitySqDeg * mFraction;
        densityTable(ii, 3) = coneDensitySqDeg * sFraction;

        % Scale macular pigment density
        macularObj.density = baseMacularDensity .* exp(-ecc / 2.0);
        pigmentObj.opticalDensity = basePigmentOpticalDensity;

        % Recalculate radiometric scalar for this eccentricity
        radiometricScalar = (pupilAreaM2 / focalLengthM^2) * coneApertureM2;

        % Base LMS absorptance and sensitivities
        baseLMS = diag(lens.transmittance) * pigmentObj.absorptance;
        baseLMS = diag(energyToQuanta) * baseLMS .* radiometricScalar;
        effectiveLMS = diag(macularObj.transmittance) * baseLMS;

        % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
        T_mat = effectiveLMS' * pinvCamSens;

        % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
        transformTable(ii, :) = T_mat(:)';
        assert(~any(isnan(T_mat(:)')))
    end

    % Store the results
    coneMapVar(age).eccGrid = eccGrid;
    coneMapVar(age).transformTable = transformTable;
    coneMapVar(age).densityTable = densityTable;
    coneMapVar(age).basePupilDiameterMm = options.basePupilDiameterMm;
    coneMapVar(age).observerAge = age;
end
fprintf('done\n');

% Precompute pre-receptoral optics (PSFs) for the specified pupil range
fprintf('Generating ISETbio pre-receptoral optics for pupil sizes');
opticsSupport = struct();
peakWavelengths = [562, 530, 430]; % L, M, S peaks

for p = 1:length(options.pupilDiameterMmRange)
    fprintf('.');
    pupilMm = options.pupilDiameterMmRange(p);

    opticsSupport(p).pupilDiameterMm = pupilMm;

    for c = 1:3
        % Instantiate a clean wavefront object for EACH wavelength 
        % independently to bypass ISETbio multi-wavelength bugs.
        wvf = wvfCreate();
        wvf = wvfSet(wvf, 'calc wavelengths', peakWavelengths(c));

        % The measured pupil size must remain at 8.0mm to preserve the correct 
        % scaling of Thibos higher-order biological aberrations.
        wvf = wvfSet(wvf, 'measured pupil diameter', 8.0);
        wvf = wvfSet(wvf, 'calc pupil diameter', pupilMm);

        % Provide sufficient spatial resolution
        wvf = wvfSet(wvf, 'spatial samples', 401);

        % Correct spherical refractive error (defocus)
        wvf = wvfSet(wvf, 'zcoeffs', 0, {'defocus'});

        % Compute the wavefront and PSF for this single wavelength
        wvf = wvfCompute(wvf);

        % Extract the 2D PSF and spatial sampling resolution using 'min' for arc minutes
        opticsSupport(p).psf{c} = wvfGet(wvf, 'psf', peakWavelengths(c));
        psfSupport = wvfGet(wvf, 'psf spatial samples', 'min', peakWavelengths(c));
        opticsSupport(p).psfSpacingArcMin(c) = abs(psfSupport(2) - psfSupport(1));
    end
end
fprintf('done\n');

% Save both variables in the derived directory
saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radianceToConeRateSupport.mat');
readme = ['Created by defineSensorToConeMapping.\n'...
    'coneMapVar -- a structure with values needed for conversion of radiance to cone isomerization rate.\n'...
    'opticsSupport -- a structure containing pre-computed PSFs for a range of pupil diameters.\n'];
save(saveFileName,'readme','coneMapVar','opticsSupport');