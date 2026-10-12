% defineSensorToConeMappingConeMosaic.m

% This script uses ISETbio to define a variable that supports the
% conversion of integrated radiance values from the camera sensors to cone
% isomerization rates for an observer. We loop over a range of ages (and
% thus lens density parameters) and pupil diameters to precompute optics.

% Housekeeping
clear

% Hard code some options
% Meridian: Select between {'Upper', 'Lower', 'Nasal', 'Temporal'}
options.meridian= 'Nasal';
options.visualizeEmployedConeMosaics = false;
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


% Precompute the pseudo-inverse of the camera sensitivities
pinvCamSens = pinv(camSens');

% Create a 1D grid of eccentricities out to the edges of the visual field
% eccGrid = 0:0.25:120;

eccGrid = 0:1:20;

% Preallocate a struct array with empty fields for the cone mapping
coneMapVar(options.ageRange(2)) = struct(...
    'wavelengthSupportNM', [], ...
    'eccGrid', [], ...
    'transformTable', [], ...
    'effectiveLMSTable', [], ...
    'densityTable', [], ...
    'basePupilDiameterMm', [],...
    'observerAge', []);

% Set this to true to employ the @coneMosaic based derivation of the
% effectiveLMS. Set to false to use the previous, first-principles based
% approach.
useCMosaicMethod = true;

% Loop over observer ages
for age = options.ageRange(1):options.ageRange(1)
    fprintf('.');

    % Apply the standard CIE 2006 age-dependent scaling factor to the lens density
    if age <= 60
        ageScalar = 1 + 0.02 * (age - 32);
    else
        ageScalar = 1.56 + 0.0667 * (age - 60);
    end
    lensDensity = ageScalar;

    if (useCMosaicMethod)
        theConeMapVar = generateEffectiveLMSusingCmosaicMethod(wlsSensor, eccGrid, lensDensity, pinvCamSens, options, age);
    else
        theConeMapVar = generateEffectiveLMSusingOldMethod(wlsSensor, eccGrid, lensDensity, pinvCamSens, options, age);
    end

    coneMapVar(age) = theConeMapVar;
end % age


%
% Have not yet worked on the optics part
%

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

if (useCMosaicMethod)
    % Save both variables in the derived directory
    saveFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        sprintf('radianceToConeRateSupportCMosaicBased%sMeridian.mat', options.meridian));
else
    % Save both variables in the derived directory
    saveFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radianceToConeRateSupport.mat');
end


readme = ['Created by defineSensorToConeMappingConeMosaic.\n'...
    'coneMapVar -- a structure with values needed for conversion of radiance to cone isomerization rate.\n'...
    'opticsSupport -- a structure containing pre-computed PSFs for a range of pupil diameters.\n'];
save(saveFileName,'readme','coneMapVar','opticsSupport');

fprintf('Saved data to %s\n', saveFileName);

%
% SUPPORTING ROUTINES
%

function theConeMapVar = generateEffectiveLMSusingCmosaicMethod(wlsSensor, eccGrid, lensDensity, pinvCamSens, options, age)

    transformTable = zeros(length(eccGrid), 9);
    densityTable = zeros(length(eccGrid), 3);

    % Loop over scalar eccentricities to build the LUT
    for ii = 1:length(eccGrid)
        fprintf('Probing cone mosaic at ecc: %2.1f degs (%s meridian)\n', eccGrid(ii), options.meridian);
        ecc = eccGrid(ii);
        
        safeEccAperture = min(ecc, 60);
        eccMeters = safeEccAperture * (300 * 1e-6);

        % Obtain an estimate of the cone aperture at the examined eccentricity
        % in order to define an appropriate size for the employed cone mosaic 
        [~, apertureDiameterMeters] = coneSizeReadData('eccentricity', eccMeters, 'angle', 0);

        % Determine cone mosaic size (degs) based on current eccentricity
        if (ecc == 0)
            coneMosaicSizeDegs = 0.25;
        else
            coneMosaicSizeDegs = -0.5 + 100*(apertureDiameterMeters*1e6/300);
        end
       
        if (strcmp(options.meridian, 'Nasal') && (abs(ecc-15.5)<= 3))
            % For nasal eccentricities (ecc > 0), around 15.5, where the optic disk is located,
            % enlarge the cone mosaic size to ensure that we have some cones in the generated cone mosaic
            % to compute the effectiveLMS
            coneMosaicSizeDegs = 4;
        end

        switch (options.meridian)
            case 'Nasal'
                eccentricityDegs = [ecc 0];
                
            case 'Temporal'
                eccentricityDegs = [-ecc 0];  

            case 'Upper'
                eccentricityDegs = [0 ecc];

            case 'Lower'
                eccentricityDegs = [0 -ecc];
            otherwise
                error('Uknown meridian: ''%s''.', options.meridian)
        end

        % Generate cone mosaic at the desired eccentricity
        theConeMosaic = cMosaic(...
            'eccentricityDegs', eccentricityDegs, ...
            'sizeDegs', coneMosaicSizeDegs*[1 1], ...
            'wave', wlsSensor, ...
            'sourceLatticeSizeDegs', 64, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', struct('smoothLocalVariations', true));

        if (options.visualizeEmployedConeMosaics)
            hFig = figure(8000+ii);
            ax = subplot(1,1,1);
            theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax);
            drawnow;
        end

        % Compute and store the L,M,S cone densities
        coneDensitiesConesPerSquaredDegree = computeConeDensity(theConeMosaic);
        densityTable(ii, 1) = coneDensitiesConesPerSquaredDegree(1);
        densityTable(ii, 2) = coneDensitiesConesPerSquaredDegree(2);
        densityTable(ii, 3) = coneDensitiesConesPerSquaredDegree(3);

        % Generate diffraction-limited optics
        micronsPerDegree = theConeMosaic.micronsPerDegree;
        opticalDefocusDiopters = 0;
        imagingPeakWavelength = 550;

        [theDiffractionLimitedOptics, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = ...
            diffractionLimitedOptics(options.basePupilDiameterMm, wlsSensor, imagingPeakWavelength, ...
                micronsPerDegree, opticalDefocusDiopters, ...
                'wavefrontSpatialSamples', 101, ...
                'noLCA', true);

        % Apply the desired lens density
        theOptics = oiGet(theDiffractionLimitedOptics, 'optics');
        theLens = opticsGet(theOptics, 'lens');
    
        % Lens with desired density
        theUpdatedLens = theLens;
        theUpdatedLens.density = lensDensity;
        theUpdatedOptics = opticsSet(theOptics, 'lens', theUpdatedLens);
        theDiffractionLimitedOptics = oiSet(theDiffractionLimitedOptics, 'optics', theUpdatedOptics);
    
        % Determine the cones we will sample for computing the effectiveLMS
        %reportLMSconePositions = false;
        %[theLconeIndex, theMconeIndex, theSconeIndex] = ...
        %    determineSampledLMSconeIndices(theConeMosaic, reportLMSconePositions);

        % Stimulus params
        pixelsNum = 5; % uniform scene, so we dont need many pixels
        theEnergyWattsPerSteradianPerM2PerNM = 1.0 * 1E-3;  % 1.0 mWatts/Sr/M2/nm

        % Preallocate memory
        quantalExcitationRateSpectraFromCMosaic = zeros(numel(wlsSensor),3);

        % Probe the cone mosaic with monochromatic stimuli to extract the effectiveLMS 
        parfor iWave = 1:numel(wlsSensor)

            monochromaticBandsNM = wlsSensor(iWave);
            theScene = generateMonochromaticScene(...
                max(theConeMosaic.sizeDegs), pixelsNum, wlsSensor, ...
                theEnergyWattsPerSteradianPerM2PerNM, monochromaticBandsNM);
    
            % Compute average LMS excitation spectra using the employed coneMosaic
            coneAbsorptionsPhotonsPerConePerSecFromCmosaicSingleCone = ...
                computeExcitationRatesUsingTheConeMosaic(...
                    theScene, theDiffractionLimitedOptics, ...
                    theConeMosaic);

            % Divide the spectra by the stimulus energy
            for iCone = 1:3
                quantalExcitationRateSpectraFromCMosaic(iWave,iCone) = ...
                    mean(coneAbsorptionsPhotonsPerConePerSecFromCmosaicSingleCone(:,iCone),1) / ...
                    theEnergyWattsPerSteradianPerM2PerNM;
            end
        end % parfor iWave
    

        effectiveLMS = quantalExcitationRateSpectraFromCMosaic;
    

        % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
        T_mat = effectiveLMS' * pinvCamSens;
    
        % Store the fundamentals
        effectiveLMSTable(ii,:,:) = effectiveLMS;

        % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
        transformTable(ii, :) = T_mat(:)';
        assert(~any(isnan(T_mat(:)')))
    end % ii

    % Store the results
    theConeMapVar.wavelengthSupportNM = wlsSensor; 
    theConeMapVar.effectiveLMSTable = effectiveLMSTable;
    theConeMapVar.eccGrid = eccGrid;
    theConeMapVar.transformTable = transformTable;
    theConeMapVar.densityTable = densityTable;
    theConeMapVar.basePupilDiameterMm = options.basePupilDiameterMm;
    theConeMapVar.observerAge = age;

end

%
% Method to compute the effectiveLMS by interrogating cone mosaics
%

function coneAbsorptionsPhotonsPerConePerSec = computeExcitationRatesUsingTheConeMosaic(...
    theScene, theOptics, theConeMosaic)

    % Compute the retinal image of the scene
    theRetinalImage = oiCompute(theOptics, theScene);

    % Compute the L-, M-, and S-cone exciations to the monochromatic scene
    theConeExcitations = theConeMosaic.compute(...
            theRetinalImage, ...
            'nTrials', 1);
   % assert(theConeMosaic.coneTypes(theLconeIndex) == cMosaic.LCONE_ID, 'did not find an L-cone');
   % assert(theConeMosaic.coneTypes(theMconeIndex) == cMosaic.MCONE_ID, 'did not find an M-cone');
   % assert(theConeMosaic.coneTypes(theSconeIndex) == cMosaic.SCONE_ID, 'did not find an S-cone');

    % From counts to rates
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 1) = ...
        mean(theConeExcitations(1,1,theConeMosaic.lConeIndices),3) / theConeMosaic.integrationTime;
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 2) = ...
        mean(theConeExcitations(1,1,theConeMosaic.mConeIndices),3) / theConeMosaic.integrationTime;
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 3) = ...
        mean(theConeExcitations(1,1,theConeMosaic.sConeIndices),3) / theConeMosaic.integrationTime;
end

%
% Method to determine the cone densities from the examined cone mosaics.
%
function coneDensitiesConesPerSquaredDegree = computeConeDensity(theConeMosaic)
    
    distancesFromConeMosaicCenter = sqrt(sum(bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, mean(theConeMosaic.coneRFpositionsDegs,1)).^2,2));
    
    theRadiusDegs = 0.5*min(theConeMosaic.sizeDegs);
    theConeIndices = find(distancesFromConeMosaicCenter <= theRadiusDegs);
    
    idx = find(theConeMosaic.coneTypes(theConeIndices) == cMosaic.LCONE_ID);
    theLconeIndices = theConeIndices(idx);

    idx = find(theConeMosaic.coneTypes(theConeIndices) == cMosaic.MCONE_ID);
    theMconeIndices = theConeIndices(idx);

    idx = find(theConeMosaic.coneTypes(theConeIndices) == cMosaic.SCONE_ID);
    theSconeIndices = theConeIndices(idx);

    fprintf('Will compute density based on %d (out of %d) cones (%d+%d+%d) within a radius of %2.2fdegs\n', ...
        numel(theConeIndices), theConeMosaic.conesNum, ...
        numel(theLconeIndices), numel(theMconeIndices), numel(theSconeIndices), theRadiusDegs);

    theArea = pi * theRadiusDegs^2;
    
    coneDensitiesConesPerSquaredDegree(1) = numel(theLconeIndices)/theArea;
    coneDensitiesConesPerSquaredDegree(2) = numel(theMconeIndices)/theArea;
    coneDensitiesConesPerSquaredDegree(3) = numel(theSconeIndices)/theArea;
end

%
% Method to determine the center-most L, M, and S-cone in a given conem mosaic
%
function [theLconeIndex, theMconeIndex, theSconeIndex] = determineSampledLMSconeIndices(theConeMosaic, reportLMSconePositions)
    theLconeIndices = theConeMosaic.lConeIndices;
    radialEccs = sqrt(sum(theConeMosaic.coneRFpositionsDegs(theLconeIndices,:).^2,2));
    [~, idx] = sort(radialEccs);
    theLconeIndex = theLconeIndices(idx(1));

    theMconeIndices = theConeMosaic.mConeIndices;
    radialEccs = sqrt(sum(theConeMosaic.coneRFpositionsDegs(theMconeIndices,:).^2,2));
    [~, idx] = sort(radialEccs);
    theMconeIndex = theMconeIndices(idx(1));

    theSconeIndices = theConeMosaic.sConeIndices;
    radialEccs = sqrt(sum(theConeMosaic.coneRFpositionsDegs(theSconeIndices,:).^2,2));
    [~, idx] = sort(radialEccs);
    theSconeIndex = theSconeIndices(idx(1));

    if (reportLMSconePositions)
        fprintf('L-cone at (%2.2f, %2.2f) microns\n',...
            theConeMosaic.coneRFpositionsMicrons(theLconeIndex,1), ...
            theConeMosaic.coneRFpositionsMicrons(theLconeIndex,2));
    
        fprintf('M-cone at (%2.2f, %2.2f) microns\n',...
            theConeMosaic.coneRFpositionsMicrons(theMconeIndex,1), ...
            theConeMosaic.coneRFpositionsMicrons(theMconeIndex,2));
    
        fprintf('S-cone at (%2.2f, %2.2f) microns\n',...
            theConeMosaic.coneRFpositionsMicrons(theSconeIndex,1), ...
            theConeMosaic.coneRFpositionsMicrons(theSconeIndex,2));
    end

end

%
%  Method to compute a monochromatic scene with which we extract the cone
%  fundamentals directly from the cone mosaic
%
function theScene = generateMonochromaticScene(fovDegs, pixelsNum, ...
    spectralSupportNanoMeters, theEnergyWattsPerSteradianPerM2PenNM, targetWavelengths)
    % Generate equivalent ISETBio scene
    theScene  = sceneCreate('uniform bb', pixelsNum, 5000, spectralSupportNanoMeters);
    
    theScene = sceneSet(theScene,'fov', fovDegs);

    % Rewrite the scene with the desired monochromatic spectrum
    % allowing photons only at the current monomchromatic wavelength
    energyWithinBand = theEnergyWattsPerSteradianPerM2PenNM;
    energy = sceneGet(theScene,'energy')*0;

    for iMonoChromaticBand = 1:numel(spectralSupportNanoMeters)
        if (ismember(spectralSupportNanoMeters(iMonoChromaticBand), targetWavelengths))
            energy(:,:,iMonoChromaticBand) = energyWithinBand * ones(pixelsNum, pixelsNum);
        end
    end
    theScene = sceneSet(theScene,'energy',energy);
end



function theConeMapVar = generateEffectiveLMSusingOldMethod(wlsSensor, eccGrid, lensDensity, pinvCamSens, options, age)

    % Instantiate ISETbio lens directly
    lens = Lens('wave', wlsSensor);
    lens.density = lensDensity;

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

    transformTable = zeros(length(eccGrid), 9);
    densityTable = zeros(length(eccGrid), 3); % Store L, M, S cone densities

    % Report we are starting age loop
    fprintf('Generating %d ISETbio age profiles', range(options.ageRange))

    % Loop over scalar eccentricities to build the LUT
    for ii = 1:length(eccGrid)
        fprintf('Analyzing ecc: %f\n', eccGrid(ii));
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

        if (1==2)
            figure(1); clf;
            plot(wlsSensor, effectiveLMS(:,1), 'r-');
            hold on;
            plot(wlsSensor, effectiveLMS(:,2), 'g-');
            plot(wlsSensor, effectiveLMS(:,3), 'b-');
            
            save('effectiveLMSoldmethod.mat', 'effectiveLMS', 'wlsSensor');
            disp('here')
            pause;
        end

        % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
        T_mat = effectiveLMS' * pinvCamSens;

        % Store the fundamentals
        effectiveLMSTable(ii,:,:) = effectiveLMS;

        % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
        transformTable(ii, :) = T_mat(:)';
        assert(~any(isnan(T_mat(:)')))
    end

    
    % Store the results
    theConeMapVar.wavelengthSupportNM = wlsSensor; 
    theConeMapVar.effectiveLMSTable = effectiveLMSTable;
    theConeMapVar.eccGrid = eccGrid;
    theConeMapVar.transformTable = transformTable;
    theConeMapVar.densityTable = densityTable;
    theConeMapVar.basePupilDiameterMm = options.basePupilDiameterMm;
    theConeMapVar.observerAge = age;
end

