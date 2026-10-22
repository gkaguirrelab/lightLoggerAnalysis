% defineSensorToConeMappingConeMosaicWithOptics.m

% This script uses ISETbio to define a variable that supports the
% conversion of integrated radiance values from the camera sensors to cone
% isomerization rates for an observer. We loop over a range of ages (and
% thus lens density parameters) and pupil diameters to precompute optics.

% Housekeeping
clear

% Hard code some options
options.optics.zernikeDataBase = 'Polans2015';

options.optics.subjectID = 5;

% 1D (defocus, oblique & vertical astigmatism) optics optimization 
options.optics.StrehlRatioOptimizationParams = struct(...
       'defocusDioptersRange', -3:0.25:3, ...
       'obliqueAstigmatismDioptersRange', nan, ...
       'verticalAstigmatismDioptersRange', nan ...
   );

% 3D (defocus, oblique & vertical astigmatism) optics optimization 
% options.optics.StrehlRatioOptimizationParams = struct(...
%         'defocusDioptersRange', -2:0.25:2, ...
%         'obliqueAstigmatismDioptersRange', -2:0.25:2, ...
%         'verticalAstigmatismDioptersRange', -2:0.25:2 ...
%     );


% 601 (or more) spatial samples are required to support the 
% extent of PSFs at extreme eccentricities (> 30 degs), 
% especially at short wavelengths
options.optics.wavefrontSpatialSamples = 601;
options.ageRange  = 18;
options.basePupilDiameterMm = 3.0;
options.pupilDiameterMmRange = 3.0; %1.5:0.5:4.0; % 1.5:0.5:8.0; % Range for pre-computing PSFs

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

    % 2D eccentricity grid
    deltaEccGridDeg = 10;
    [X,Y] = meshgrid(-40:deltaEccGridDeg:40, -20:deltaEccGridDeg:20);
    eccGrid = [X(:) Y(:)];

    % Add some extreme eccs now
    eccGrid(size(eccGrid,1)+1,:) = [-40  -10];
    eccGrid(size(eccGrid,1)+1,:) = [-40   0];
    eccGrid(size(eccGrid,1)+1,:) = [-40   10];

    eccGrid(size(eccGrid,1)+1,:) = [40  -10];
    eccGrid(size(eccGrid,1)+1,:) = [40   0];
    eccGrid(size(eccGrid,1)+1,:) = [40   10];

    eccGrid = eccGrid';

% To speed up things lets subsample the wavelenth axis;
deltaWls = 10;
wlsSensorSubsample = wlsSensor(1):deltaWls:wlsSensor(end);

% Compute defocus and astigmatism for maximizing the foveal Strehl ratio
opticsSupport = computeOptimalStrehlRatioDefocusAndAstigmatismForFovealOptics(options, wlsSensorSubsample);

% Preallocate a struct array with empty fields for the cone mapping
coneMapVar(options.ageRange(1)) = struct(...
    'theEmployedConeMosaics', [], ...
    'supportWavelength', [], ...
    'eccGrid', [], ...
    'transformTable', [], ...
    'effectiveLMSTable', [], ...
    'densityTable', [], ...
    'basePupilDiameterMm', [],...
    'observerAge', []);

% Loop over observer ages
for age = options.ageRange(1):options.ageRange(1)

    % Apply the standard CIE 2006 age-dependent scaling factor to the lens density
    if age <= 60
        ageScalar = 1 + 0.02 * (age - 32);
    else
        ageScalar = 1.56 + 0.0667 * (age - 60);
    end
    lensDensity = ageScalar;

        theConeMapVar = generateEffectiveLMSusingCmosaicMethod(wlsSensor, deltaWls, eccGrid, lensDensity, pinvCamSens, options, age);
    coneMapVar(age) = theConeMapVar;
end


% Compute cone-weighted PSFs for one particular lens age effectiveLMSTable
theConeMapVar = coneMapVar(options.ageRange(1));
    
    opticsSupport.thePSFs = generateOpticsUsingCmosaicMethod(...
            theConeMapVar.theEmployedConeMosaics, theConeMapVar.supportWavelength, theConeMapVar.effectiveLMSTable, options, opticsSupport.theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters);
    % Save both variables in the derived directory
    saveFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'radianceToConeSupport.mat');

readme = ['Created by defineSensorToConeMappingConeMosaic.\n'...
    'coneMapVar -- a structure with values needed for conversion of radiance to cone isomerization rate.\n'...
    'opticsSupport -- a structure containing pre-computed PSFs for a range of pupil diameters.\n'];
save(saveFileName,'readme','coneMapVar','opticsSupport', '-v7.3');

fprintf('Saved data to %s\n', saveFileName);

%
% SUPPORTING ROUTINES
%


function opticsSupport = computeOptimalStrehlRatioDefocusAndAstigmatismForFovealOptics(options, wlsSensorSubsample)
    opticsParams.zernikeDataBase = options.optics.zernikeDataBase;
    opticsParams.subjectID = options.optics.subjectID;
    opticsParams.StrehlRatioOptimizationParams = options.optics.StrehlRatioOptimizationParams;

    opticsParams.wavefrontSpatialSamples = options.optics.wavefrontSpatialSamples;
    opticsParams.pupilDiameterMM = 3.0;
    opticsParams.zeroCenterPSF = false;
    opticsParams.noLCA = false;
        
    % Generate foveal cone mosaic
    eccentricityDegs = [0 0];
    coneMosaicSizeDegs = [0.2 0.2];
    theFovealConeMosaic = cMosaic(...
            'eccentricityDegs', eccentricityDegs, ...
            'sizeDegs', coneMosaicSizeDegs, ...
            'wave', wlsSensorSubsample, ...
            'sourceLatticeSizeDegs', 64, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', struct('smoothLocalVariations', true));

    % Search over the 3D search grid to detemine the defocus and
    % astigmatism values that are needed maximize the Strehl ratio
    % We only do this in the fovea 
    psfUpsampleFactor = [];
    useParfor = true;
    visualizeStrehlRatioOptimization = true;
    contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection = true;

    fprintf('Optimizing foveal optics.\n');
            
    [~, ~, zCoeffsMicrons] = theFovealConeMosaic.oiEnsembleGenerate(theFovealConeMosaic.eccentricityDegs, ...
        'zernikeDataBase', opticsParams.zernikeDataBase, ...
        'subjectID', opticsParams.subjectID, ...
        'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
        'subtractCentralRefraction', false, ...
        'subtractCentralRefraction3D', false, ...
        'zeroCenterPSF', ~true, ...
        'withZeroedPistonAndTiltZernikeCoefficients', ~true, ...
        'flipPSFUpsideDown', false, ...
        'noLCA', opticsParams.noLCA, ...
        'refractiveErrorDiopters', 0, ...
        'obliqueAndVerticalAstigmatismErrorsDiopters', [0 0], ...
        'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples);


    % The Artal/Polans Z-coeff indices start from 0, not 1
    theDefocusZcoeffIndex = wvfOSAIndexToVectorIndex('defocus')-1;
    theObliqueAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('oblique_astigmatism')-1;
    theVerticalAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('vertical_astigmatism')-1;

    % From microns to diopters
    defocusDioptersCenter = wvfDefocusMicronsToDiopters(zCoeffsMicrons(theDefocusZcoeffIndex), opticsParams.pupilDiameterMM);
    obliqueAstigmatismDioptersCenter = wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theObliqueAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM);
    verticalAstigmatismDioptersCenter = wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theVerticalAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM);

    % Round to closest quarter diopter
    defocusDioptersCenter = sign(defocusDioptersCenter) * 0.25*round(4*abs(defocusDioptersCenter));
    obliqueAstigmatismDioptersCenter = sign(obliqueAstigmatismDioptersCenter) * 0.25*round(4*abs(obliqueAstigmatismDioptersCenter));
    verticalAstigmatismDioptersCenter = sign (verticalAstigmatismDioptersCenter) * 0.25*round(4*abs(verticalAstigmatismDioptersCenter));

    examinedDefocusDioptersRange = opticsParams.StrehlRatioOptimizationParams.defocusDioptersRange;
    examinedObliqueAstigmatismDioptersRange = opticsParams.StrehlRatioOptimizationParams.obliqueAstigmatismDioptersRange;
    examinedVerticalAstigmatismDioptersRange = opticsParams.StrehlRatioOptimizationParams.verticalAstigmatismDioptersRange;
       
    % Center the 3D search around (0,0,0) defocus and astigmatism Z coeffs

    examinedDefocusDiopters = examinedDefocusDioptersRange - defocusDioptersCenter;
    if (isnan(examinedObliqueAstigmatismDioptersRange))
        examinedObliqueAstigmatismDiopters = 0;
    else
        examinedObliqueAstigmatismDiopters = examinedObliqueAstigmatismDioptersRange - obliqueAstigmatismDioptersCenter;
    end
    if (isnan(examinedVerticalAstigmatismDioptersRange))
        examinedVerticalAstigmatismDiopters = 0;
    else
        examinedVerticalAstigmatismDiopters = examinedVerticalAstigmatismDioptersRange -verticalAstigmatismDioptersCenter;
    end

    % 3D search grid
    [X,Y,Z] = ndgrid(...
            examinedDefocusDiopters, ...
            examinedObliqueAstigmatismDiopters, ...
            examinedVerticalAstigmatismDiopters);


    examinedDefocusAndAstigmatismDioptersGrid = [...
            reshape(X, [numel(X) 1]) reshape(Y, [numel(X) 1]) reshape(Z, [numel(X) 1]) ];

    % Go !
    [theOptics, thePSF, theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters, ...
        theOptimalStrehlRatio, StrehlRatioAsAFunctionOfDefocusAndAstigmatism] = ...
          RGCMosaicConstructor.helper.optics.optimized3DStrehlRatioPSF(...
			    examinedDefocusAndAstigmatismDioptersGrid, ...
			    theFovealConeMosaic, theFovealConeMosaic.eccentricityDegs, opticsParams, ...
                opticsParams.wavefrontSpatialSamples, psfUpsampleFactor, ...
			    visualizeStrehlRatioOptimization, contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection, ...
                examinedDefocusDiopters, examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, ...
                useParfor);

    opticsSupport.pupilDiameterMmRange = options.pupilDiameterMmRange;
    opticsSupport.theBasePSF = thePSF;
    opticsSupport.theOptimalStrehlRatio = theOptimalStrehlRatio;
    opticsSupport.theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters = theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters;
    opticsSupport.StrehlRatioAsAFunctionOfDefocusAndAstigmatism = StrehlRatioAsAFunctionOfDefocusAndAstigmatism;

end

function thePSFs = generateOpticsUsingCmosaicMethod(theEmployedConeMosaics, ...
    theEffectiveLMSTableWavelengthSupport, theEffectiveLMSTable, ...
    options, theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters)

    opticsParams.zernikeDataBase = options.optics.zernikeDataBase;
    opticsParams.subjectID = options.optics.subjectID;
    opticsParams.StrehlRatioOptimizationParams = options.optics.StrehlRatioOptimizationParams;

    opticsParams.wavefrontSpatialSamples = options.optics.wavefrontSpatialSamples;
    opticsParams.pupilDiameterMM = 3.0;
    opticsParams.zeroCenterPSF = false;
    opticsParams.noLCA = false;

    % Preallocate memory for thePSFs
    thePSFs = cell(numel(theEmployedConeMosaics), length(options.pupilDiameterMmRange));

    psfUpsampleFactor = [];
    useParfor = true;
    visualizeStrehlRatioOptimization = ~true;
    contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection = ~true;
    
    for ii = 1:numel(theEmployedConeMosaics)

        theConeMosaic = theEmployedConeMosaics{ii};
        theConeFundamentalsAtCurrentEccentricity = squeeze(theEffectiveLMSTable(ii,:,:));

        examinedDefocusDiopters = theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters(1);
        examinedObliqueAstigmatismDiopters = theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters(2);
        examinedVerticalAstigmatismDiopters = theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters(3);

        % 3D search grid
        [X,Y,Z] = ndgrid(...
                    examinedDefocusDiopters, ...
                    examinedObliqueAstigmatismDiopters, ...
                    examinedVerticalAstigmatismDiopters);
        
        examinedDefocusAndAstigmatismDioptersGrid = [...
                    reshape(X, [numel(X) 1]) reshape(Y, [numel(X) 1]) reshape(Z, [numel(X) 1]) ];


        % Interpolated cone fundamentals: same for all pupil sizes, only
        % computed for the first pupil size
        theInterpolatedConeFundamentals = [];

        for iPupilSize = 1:length(options.pupilDiameterMmRange)
            opticsParams.pupilDiameterMM = options.pupilDiameterMmRange(iPupilSize);

            % Compute the full spectra PSF
            [~, thePSF] = RGCMosaicConstructor.helper.optics.optimized3DStrehlRatioPSF(...
					        examinedDefocusAndAstigmatismDioptersGrid, ...
					        theConeMosaic, theConeMosaic.eccentricityDegs, opticsParams, ...
                            opticsParams.wavefrontSpatialSamples, psfUpsampleFactor, ...
					        visualizeStrehlRatioOptimization, contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection, ...
                            examinedDefocusDiopters, examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters, ...
                            useParfor);
   
            if (isempty(theInterpolatedConeFundamentals))
                theInterpolatedConeFundamentals = nan(numel(thePSF.supportWavelength), 3);
            end

            % We do not store the full-spectra PSFs (this would required >60 GB data), 
            % we only store the L/M/S-cone--weighted PSFs, with
            % the weights taken from theEffectiveLMSTable at the current
            % eccentricity and current lens age
            for iConeType = 1:3
                [thePSF.coneWeightedData(iConeType, :,:), theInterpolatedConeFundamentals(:,iConeType)] = fullSpectralPSFtoConeWeightedPSF(...
                    theConeFundamentalsAtCurrentEccentricity(:,iConeType), theEffectiveLMSTableWavelengthSupport, ...
                    theInterpolatedConeFundamentals(:,iConeType), thePSF.data, thePSF.supportWavelength);
            end

            % Zero the full spectral PSF to save memory
            thePSF.data = [];
            thePSF.supportWavelength = [];

            % Save the cone weighted data
            thePSFs{ii,iPupilSize} = thePSF;
        end % for iPupilSize

    end % for ii
end

function [theConePSF, theInterpolatedConeFundamental] = fullSpectralPSFtoConeWeightedPSF(...
    theConeFundamental, theConeFundamentalWavelengthSupport, theInterpolatedConeFundamental, ...
    thePSFdata, thePSFWavelengthSupport)
    
    if (any(isnan(theInterpolatedConeFundamental(:))))
        % Interpolate cone fundamental to the spectral support of the PSF
        theInterpolatedConeFundamental = interp1(theConeFundamentalWavelengthSupport, ...
            theConeFundamental, thePSFWavelengthSupport);
    end
    
    % Reshape the full-spectra PSF for efficient computation of theConePSF
    wavelengthsNum = size(thePSFdata,3);
    rowsNum = size(thePSFdata,1);
    colsNum = size(thePSFdata,2);
    thePSFdata = reshape(thePSFdata, [rowsNum*colsNum wavelengthsNum]);

    % Compute theConePSF
    theConePSF = reshape(thePSFdata * theInterpolatedConeFundamental, [rowsNum colsNum]);
end


function theConeMapVar = generateEffectiveLMSusingCmosaicMethod(wlsSensor, deltaWls, eccGrid, lensDensity, pinvCamSens, options, age)

    theEmployedConeMosaics = cell(1,length(eccGrid));
    transformTable = zeros(length(eccGrid), 9);
    densityTable = zeros(length(eccGrid), 3);

    % Loop over scalar eccentricities to build the LUT
    for ii = 1:length(eccGrid)
        theEccDegs = eccGrid(:,ii);
        
            % 2D grid
            fprintf('Cone mosaic at [%2.1f,%2.1f] degs...', eccGrid(1,ii), eccGrid(2,ii));

            radialEccDegs = norm(theEccDegs);
            safeEccApertureDegs = min(radialEccDegs, 60);
            radialEccMeters = safeEccApertureDegs * (300 * 1e-6);
    
            % Obtain an estimate of the cone aperture at the examined eccentricity
            % in order to define an appropriate size for the employed cone mosaic 
            [~, apertureDiameterMeters] = coneSizeReadData('eccentricity', radialEccMeters, 'angle', 0);
    
            % Determine cone mosaic size (degs) based on current eccentricity
            if (radialEccDegs < 0.1)
                coneMosaicSizeDegs = 0.25;
            else
                coneMosaicSizeDegs = -0.5 + 100*(apertureDiameterMeters*1e6/300);
            end

            angularEccDegs = atan2d(theEccDegs(2), theEccDegs(1));
            if (abs(angularEccDegs) < 15) && (abs(radialEccDegs-15.5)<= 3)
                % For nasal eccentricities around 15.5, where the optic disk is located,
                % enlarge the cone mosaic size to ensure that we have some cones in the generated cone mosaic
                % to compute the effectiveLMS
                coneMosaicSizeDegs = 4;
            end

            eccentricityDegs = [theEccDegs(1) theEccDegs(2)];
        
        % To speed up things lets subsample the wavelength axis;
        wlsSensorSubsample = wlsSensor(1):deltaWls:wlsSensor(end);
        
        % Generate cone mosaic at the desired eccentricity
        theConeMosaic = cMosaic(...
            'eccentricityDegs', eccentricityDegs, ...
            'sizeDegs', coneMosaicSizeDegs*[1 1], ...
            'wave', wlsSensorSubsample, ...
            'sourceLatticeSizeDegs', 64, ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'coneApertureModifiers', struct('smoothLocalVariations', true));

        % Compute and store the L,M,S cone densities
        coneDensitiesConesPerSquaredDegree = computeConeDensity(theConeMosaic);
        densityTable(ii, 1) = coneDensitiesConesPerSquaredDegree(1);
        densityTable(ii, 2) = coneDensitiesConesPerSquaredDegree(2);
        densityTable(ii, 3) = coneDensitiesConesPerSquaredDegree(3);

        % Generate diffraction-limited optics
        micronsPerDegree = theConeMosaic.micronsPerDegree;
        opticalDefocusDiopters = 0;
        imagingPeakWavelength = 550;

        theDiffractionLimitedOptics = diffractionLimitedOptics(...
            options.basePupilDiameterMm, wlsSensorSubsample, imagingPeakWavelength, ...
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
    
        % Stimulus params
        pixelsNum = 5; % uniform scene, so we dont need many pixels
        theEnergyWattsPerSteradianPerM2PerNM = 1.0 * 1E-3;  % 1.0 mWatts/Sr/M2/nm

        % Preallocate memory
        quantalExcitationRateSpectraFromCMosaic = zeros(numel(wlsSensorSubsample),3);

        % Probe the cone mosaic with monochromatic stimuli to extract the effectiveLMS 
        parfor iWave = 1:numel(wlsSensorSubsample)

            monochromaticBandsNM = wlsSensorSubsample(iWave);
            theScene = generateMonochromaticScene(...
                max(theConeMosaic.sizeDegs), pixelsNum, wlsSensorSubsample, ...
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


        % Back to original sampling
        effectiveLMS = zeros(numel(wlsSensor),3);
        for iCone = 1:3
            effectiveLMS(:,iCone) = ...
                interp1(wlsSensorSubsample, quantalExcitationRateSpectraFromCMosaic(:,iCone), wlsSensor);
        end

        % Calculate 3x3 transformation matrix from Camera RGB to LMS rates
        T_mat = effectiveLMS' * pinvCamSens;
    
        % Preallocate memory for the fundamentals
        if (ii==1)
            effectiveLMSTable = zeros(...
                length(eccGrid), ...
	            size(effectiveLMS,1), ...
	            size(effectiveLMS,2) );
        end

        % Store the fundamentals
        effectiveLMSTable(ii,:,:) = effectiveLMS;

        % Flatten the 3x3 matrix (column-major) into a 1x9 row for the LUT
        transformTable(ii, :) = T_mat(:)';
        assert(~any(isnan(T_mat(:)')))

        % Save the employed cone mosaics
        theEmployedConeMosaics{ii} = theConeMosaic;
    end % ii


    % Store the results
    theConeMapVar.theEmployedConeMosaics = theEmployedConeMosaics;
    theConeMapVar.supportWavelength = wlsSensor; 
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

    theArea = pi * theRadiusDegs^2;
    
    coneDensitiesConesPerSquaredDegree(1) = numel(theLconeIndices)/theArea;
    coneDensitiesConesPerSquaredDegree(2) = numel(theMconeIndices)/theArea;
    coneDensitiesConesPerSquaredDegree(3) = numel(theSconeIndices)/theArea;
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

