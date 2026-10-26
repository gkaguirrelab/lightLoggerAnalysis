% defineSensorToConeMapping.m
%
% This script uses ISETbio to define a variable that supports the
% conversion of integrated radiance values from the camera sensors to cone
% isomerization rates for an observer. We loop over a range of ages (and
% thus lens density parameters) and pupil diameters to precompute optics.

%% 1. CONFIGURATION
clear;

% Optics & Subject Parameters
options.optics.zernikeDataBase = 'Polans2015';
options.optics.subjectID = 2;
options.optics.wavefrontSpatialSamples = 601; % 601+ supports extreme eccs (>30 deg)
options.ageRange = 18;
options.basePupilDiameterMm = 3.0;
options.pupilDiameterMmRange = 3.0; % 1.5:0.5:8.0 for full pre-computation

% 1D (defocus, oblique & vertical astigmatism) optics optimization
options.optics.StrehlRatioOptimizationParams = struct(...
    'defocusDioptersRange', -3:0.25:3, ...
    'obliqueAstigmatismDioptersRange', nan, ...
    'verticalAstigmatismDioptersRange', nan ...
);

% Grid & Sampling Parameters
% Grid & Sampling Parameters
options.numLogSamplesX = 7; % Number of positive X-axis samples (excluding 0)
options.numLogSamplesY = 4; % Number of Y-axis samples per hemisphere (excluding 0)
options.deltaWls = 10; % sampling of the spectra
options.imagingPeakWavelength = 550;
options.pixelsNum = 5; % Uniform scene pixel count
options.energyWattsPerSteradian = 1.0 * 1E-3; % 1.0 mWatts/Sr/M2/nm

%% 2. INITIALIZATION & GRID SETUP

% Load the IMX219 sensitivity functions
dataFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'data',...
    'IMX219_spectralSensitivity.mat');
load(dataFileName,'T');
wlsSensor = T.wls;
channelNames = {'red','green','blue'};

% Extract the sensitivity vectors into an N-wave x 3 array and pseudo-inverse
camSens = [T.(channelNames{1}), T.(channelNames{2}), T.(channelNames{3})];
pinvCamSens = pinv(camSens');

% Subsample wavelength axis to improve execution time
wlsSensorSubsample = wlsSensor(1):options.deltaWls:wlsSensor(end);

% 2D eccentricity grid construction (Log-spaced, restricted to X >= 0)
% Total locations = (1 + options.numLogSamplesX) * (1 + 2 * options.numLogSamplesY)

% Create log-spaced X coordinates (0 to 40)
x_pos = logspace(log10(1), log10(40), options.numLogSamplesX);
x_vals = [0, x_pos];

% Create log-spaced Y coordinates (-20 to 20)
y_pos = logspace(log10(1), log10(20), options.numLogSamplesY);
y_vals = [-fliplr(y_pos), 0, y_pos];

% Generate meshgrid and format as 2xN array
[X, Y] = meshgrid(x_vals, y_vals);
eccGrid = [X(:), Y(:)]';


%% 3. OPTICS & CONE MAPPING EXECUTION

% Compute defocus and astigmatism for maximizing the foveal Strehl ratio
opticsSupport = computeOptimalStrehlRatioDefocusAndAstigmatismForFovealOptics(options, wlsSensorSubsample);

% Preallocate a struct array for the cone mapping
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
for age = options.ageRange(1):options.ageRange(end)

    if age <= 60
        lensDensity = 1 + 0.02 * (age - 32);
    else
        lensDensity = 1.56 + 0.0667 * (age - 60);
    end

    % Call the unified routine
    [theConeMapVar, thePSFs] = computeGridMappingAndOptics(...
        wlsSensor, wlsSensorSubsample, eccGrid, lensDensity, pinvCamSens, options, age, ...
        opticsSupport.theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters);

    coneMapVar(age) = theConeMapVar;

    if age == options.ageRange(1)
        opticsSupport.thePSFs = thePSFs;
    end
end


%% 4. EXPORT DERIVED VARIABLES

saveFileName = fullfile(...
    tbLocateProjectSilent('lightLoggerAnalysis'),...
    'derived',...
    'radianceToConeSupport.mat');

readme = ['Created by defineSensorToConeMapping.\n'...
    'coneMapVar -- a structure with values needed for conversion of radiance to cone isomerization rate.\n'...
    'opticsSupport -- a structure containing pre-computed PSFs for a range of pupil diameters.\n'];
save(saveFileName,'readme','coneMapVar','opticsSupport', '-v7.3');
fprintf('Saved data to %s\n', saveFileName);


%% SUPPORTING ROUTINES

function opticsSupport = computeOptimalStrehlRatioDefocusAndAstigmatismForFovealOptics(options, wlsSensorSubsample)
    opticsParams.zernikeDataBase = options.optics.zernikeDataBase;
    opticsParams.subjectID = options.optics.subjectID;
    opticsParams.StrehlRatioOptimizationParams = options.optics.StrehlRatioOptimizationParams;
    
    opticsParams.wavefrontSpatialSamples = options.optics.wavefrontSpatialSamples;
    opticsParams.pupilDiameterMM = options.basePupilDiameterMm;
    opticsParams.zeroCenterPSF = false;
    opticsParams.noLCA = false;
    
    theFovealConeMosaic = cMosaic(...
        'eccentricityDegs', [0 0], ...
        'sizeDegs', [0.2 0.2], ...
        'wave', wlsSensorSubsample, ...
        'sourceLatticeSizeDegs', 64, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'coneApertureModifiers', struct('smoothLocalVariations', true));
    
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
    
    theDefocusZcoeffIndex = wvfOSAIndexToVectorIndex('defocus')-1;
    theObliqueAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('oblique_astigmatism')-1;
    theVerticalAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('vertical_astigmatism')-1;
    
    defocusDioptersCenter = sign(wvfDefocusMicronsToDiopters(zCoeffsMicrons(theDefocusZcoeffIndex), opticsParams.pupilDiameterMM)) * 0.25*round(4*abs(wvfDefocusMicronsToDiopters(zCoeffsMicrons(theDefocusZcoeffIndex), opticsParams.pupilDiameterMM)));
    obliqueAstigmatismDioptersCenter = sign(wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theObliqueAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM)) * 0.25*round(4*abs(wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theObliqueAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM)));
    verticalAstigmatismDioptersCenter = sign(wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theVerticalAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM)) * 0.25*round(4*abs(wvfAstigmatismMicronsToDiopters(zCoeffsMicrons(theVerticalAstigmatismZcoeffIndex), opticsParams.pupilDiameterMM)));
    
    examinedDefocusDiopters = opticsParams.StrehlRatioOptimizationParams.defocusDioptersRange - defocusDioptersCenter;
    
    if (isnan(opticsParams.StrehlRatioOptimizationParams.obliqueAstigmatismDioptersRange))
        examinedObliqueAstigmatismDiopters = 0;
    else
        examinedObliqueAstigmatismDiopters = opticsParams.StrehlRatioOptimizationParams.obliqueAstigmatismDioptersRange - obliqueAstigmatismDioptersCenter;
    end
    
    if (isnan(opticsParams.StrehlRatioOptimizationParams.verticalAstigmatismDioptersRange))
        examinedVerticalAstigmatismDiopters = 0;
    else
        examinedVerticalAstigmatismDiopters = opticsParams.StrehlRatioOptimizationParams.verticalAstigmatismDioptersRange - verticalAstigmatismDioptersCenter;
    end
    
    [X,Y,Z] = ndgrid(examinedDefocusDiopters, examinedObliqueAstigmatismDiopters, examinedVerticalAstigmatismDiopters);
    examinedDefocusAndAstigmatismDioptersGrid = [reshape(X, [numel(X) 1]) reshape(Y, [numel(X) 1]) reshape(Z, [numel(X) 1]) ];
    
    [~, thePSF, theFovealOptimalStrehlRatioDefocusAndAstigmatismDiopters, theOptimalStrehlRatio, StrehlRatioAsAFunctionOfDefocusAndAstigmatism] = ...
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


function [theConeMapVar, gridPSFs] = computeGridMappingAndOptics(wlsSensor, wlsSensorSubsample, eccGrid, lensDensity, pinvCamSens, options, age, fovealOpticsDiopters)
    
    numEccs = size(eccGrid, 2);
    numPupils = length(options.pupilDiameterMmRange);
    
    % Preallocate arrays for parfor slicing
    theEmployedConeMosaics = cell(1, numEccs);
    transformTable = zeros(numEccs, 9);
    densityTable = zeros(numEccs, 3);
    effectiveLMSTable = zeros(numEccs, numel(wlsSensor), 3);
    gridPSFs = cell(numEccs, numPupils);
    
    % Hoist base optics using muted helper
    dummyMosaic = cMosaic('eccentricityDegs', [0 0], 'sizeDegs', [0.1 0.1]);
    micronsPerDegree = dummyMosaic.micronsPerDegree;
    theBaseOptics = generateBaseOpticsMuted(options.basePupilDiameterMm, wlsSensorSubsample, options.imagingPeakWavelength, micronsPerDegree);
    
    theOptics = oiGet(theBaseOptics, 'optics');
    theLens = opticsGet(theOptics, 'lens');
    theLens.density = lensDensity;
    theOptics = opticsSet(theOptics, 'lens', theLens);
    theBaseOptics = oiSet(theBaseOptics, 'optics', theOptics);
    
    examinedDefocusAndAstigmatismDioptersGrid = [fovealOpticsDiopters(1), fovealOpticsDiopters(2), fovealOpticsDiopters(3)];
    
    opticsParams.zernikeDataBase = options.optics.zernikeDataBase;
    opticsParams.subjectID = options.optics.subjectID;
    opticsParams.StrehlRatioOptimizationParams = options.optics.StrehlRatioOptimizationParams;
    opticsParams.wavefrontSpatialSamples = options.optics.wavefrontSpatialSamples;
    opticsParams.zeroCenterPSF = false;
    opticsParams.noLCA = false;
    
    fprintf('Processing %d grid locations in parallel...\n', numEccs);
    
    parfor ii = 1:numEccs
        
        theEccDegs = eccGrid(:,ii);
        radialEccDegs = norm(theEccDegs);
        safeEccApertureDegs = min(radialEccDegs, 60);
        radialEccMeters = safeEccApertureDegs * (300 * 1e-6);
        isOffCenter = (theEccDegs(1) > 0);
        
        [~, apertureDiameterMeters] = coneSizeReadData('eccentricity', radialEccMeters, 'angle', 0);
        
        if (radialEccDegs < 0.1)
            coneMosaicSizeDegs = 0.25;
        else
            coneMosaicSizeDegs = -0.5 + 100*(apertureDiameterMeters*1e6/300);
        end
        
        angularEccDegs = atan2d(theEccDegs(2), theEccDegs(1));
        if (abs(angularEccDegs) < 15) && (abs(radialEccDegs-15.5)<= 3)
            coneMosaicSizeDegs = 4;
        end
        
        % Format as 1x2 row vector for cMosaic compatibility
        eccentricityDegs = [theEccDegs(1), theEccDegs(2)];

        % 1 & 2. Create Mosaics, Compute Densities, and Extract LMS
        theConeMosaic = generateConeMosaicMuted(eccentricityDegs, coneMosaicSizeDegs, wlsSensorSubsample);
        effectiveLMS = computeLMSForMosaic(theConeMosaic, wlsSensorSubsample, wlsSensor, options, theBaseOptics);

        if isOffCenter
            % Generate the mirrored temporal mosaic (explicitly a 1x2 row vector)
            temporalEccDegs = [-theEccDegs(1), theEccDegs(2)];
            temporalMosaic = generateConeMosaicMuted(temporalEccDegs, coneMosaicSizeDegs, wlsSensorSubsample);

            % Average Nasal and Temporal Densities
            densN = computeConeDensity(theConeMosaic);
            densT = computeConeDensity(temporalMosaic);
            densityTable(ii, :) = (densN(1:3) + densT(1:3)) / 2;

            % Extract Temporal LMS and Average
            temporalLMS = computeLMSForMosaic(temporalMosaic, wlsSensorSubsample, wlsSensor, options, theBaseOptics);
            effectiveLMS = (effectiveLMS + temporalLMS) / 2;
        else
            % Center meridian (no mirroring needed)
            coneDensities = computeConeDensity(theConeMosaic);
            densityTable(ii, :) = coneDensities(1:3);
        end
        
        T_mat = effectiveLMS' * pinvCamSens;
        effectiveLMSTable(ii,:,:) = effectiveLMS;
        transformTable(ii, :) = T_mat(:)';
        theEmployedConeMosaics{ii} = theConeMosaic;
        
        % 3. OPTICS GENERATION
        useParforInner = false; 
        psfUpsampleFactor = [];
        visualizeStrehlRatioOptimization = false;
        contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection = false;
        localOpticsParams = opticsParams;
        
        for iPupilSize = 1:numPupils
            localOpticsParams.pupilDiameterMM = options.pupilDiameterMmRange(iPupilSize);

            % Simulate PSF optimization ONCE to save time
            thePSF = generatePSFMuted(examinedDefocusAndAstigmatismDioptersGrid, ...
                theConeMosaic, localOpticsParams, psfUpsampleFactor, ...
                visualizeStrehlRatioOptimization, contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection, ...
                fovealOpticsDiopters, useParforInner);

            theInterpolatedConeFundamentals = nan(numel(thePSF.supportWavelength), 3);
            for iConeType = 1:3
                [thePSF.coneWeightedData(iConeType, :,:), theInterpolatedConeFundamentals(:,iConeType)] = fullSpectralPSFtoConeWeightedPSF(...
                    effectiveLMS(:,iConeType), wlsSensor, ...
                    theInterpolatedConeFundamentals(:,iConeType), thePSF.data, thePSF.supportWavelength);
            end
            
            % Symmetrize off-center PSFs to capture the binocular mean without a second Strehl optimization
            if isOffCenter
                for iConeType = 1:3
                    orig = squeeze(thePSF.coneWeightedData(iConeType, :, :));
                    thePSF.coneWeightedData(iConeType, :, :) = (orig + fliplr(orig)) / 2;
                end
            end

            thePSF.data = [];
            thePSF.supportWavelength = [];
            gridPSFs{ii, iPupilSize} = thePSF;
        end
        
        fprintf('Completed grid location %d of %d\n', ii, numEccs);
    end 
    
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


function [theConePSF, theInterpolatedConeFundamental] = fullSpectralPSFtoConeWeightedPSF(...
    theConeFundamental, theConeFundamentalWavelengthSupport, theInterpolatedConeFundamental, ...
    thePSFdata, thePSFWavelengthSupport)

    if (any(isnan(theInterpolatedConeFundamental(:))))
        theInterpolatedConeFundamental = interp1(theConeFundamentalWavelengthSupport, ...
            theConeFundamental, thePSFWavelengthSupport);
    end
    
    wavelengthsNum = size(thePSFdata,3);
    rowsNum = size(thePSFdata,1);
    colsNum = size(thePSFdata,2);
    thePSFdata = reshape(thePSFdata, [rowsNum*colsNum wavelengthsNum]);
    
    theConePSF = reshape(thePSFdata * theInterpolatedConeFundamental, [rowsNum colsNum]);
end


function coneAbsorptionsPhotonsPerConePerSec = computeExcitationRatesUsingTheConeMosaic(...
    theScene, theOptics, theConeMosaic)
    theRetinalImage = oiCompute(theOptics, theScene);
    [~, theConeExcitations] = evalc("theConeMosaic.compute(theRetinalImage, 'nTrials', 1)");
    
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 1) = mean(theConeExcitations(1,1,theConeMosaic.lConeIndices),3) / theConeMosaic.integrationTime;
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 2) = mean(theConeExcitations(1,1,theConeMosaic.mConeIndices),3) / theConeMosaic.integrationTime;
    coneAbsorptionsPhotonsPerConePerSec(1, 1, 3) = mean(theConeExcitations(1,1,theConeMosaic.sConeIndices),3) / theConeMosaic.integrationTime;
end


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


function theScene = generateMonochromaticScene(fovDegs, pixelsNum, ...
    spectralSupportNanoMeters, theEnergyWattsPerSteradianPerM2PenNM, targetWavelengths)
    
    theScene  = sceneCreate('uniform bb', pixelsNum, 5000, spectralSupportNanoMeters);
    theScene = sceneSet(theScene,'fov', fovDegs);
    energyWithinBand = theEnergyWattsPerSteradianPerM2PenNM;
    energy = sceneGet(theScene,'energy')*0;
    
    for iMonoChromaticBand = 1:numel(spectralSupportNanoMeters)
        if (ismember(spectralSupportNanoMeters(iMonoChromaticBand), targetWavelengths))
            energy(:,:,iMonoChromaticBand) = energyWithinBand * ones(pixelsNum, pixelsNum);
        end
    end
    theScene = sceneSet(theScene,'energy',energy);
end


% --- PARFOR TRANSPARENCY & MODULARITY HELPERS ---

function effectiveLMS = computeLMSForMosaic(theConeMosaic, wlsSensorSubsample, wlsSensor, options, theBaseOptics)
    quantalExcitationRateSpectraFromCMosaic = zeros(numel(wlsSensorSubsample),3);
    for iWave = 1:numel(wlsSensorSubsample)
        monochromaticBandsNM = wlsSensorSubsample(iWave);
        theScene = generateMonochromaticScene(...
            max(theConeMosaic.sizeDegs), options.pixelsNum, wlsSensorSubsample, ...
            options.energyWattsPerSteradian, monochromaticBandsNM);
        
        coneAbsorptionsPhotonsPerConePerSecFromCmosaicSingleCone = ...
            computeExcitationRatesUsingTheConeMosaic(theScene, theBaseOptics, theConeMosaic);
        
        for iCone = 1:3
            quantalExcitationRateSpectraFromCMosaic(iWave,iCone) = ...
                mean(coneAbsorptionsPhotonsPerConePerSecFromCmosaicSingleCone(:,iCone),1) / ...
                options.energyWattsPerSteradian;
        end
    end 
    
    effectiveLMS = zeros(numel(wlsSensor),3);
    for iCone = 1:3
        effectiveLMS(:,iCone) = interp1(wlsSensorSubsample, quantalExcitationRateSpectraFromCMosaic(:,iCone), wlsSensor);
    end
end

function theBaseOptics = generateBaseOpticsMuted(basePupilDiameterMm, wlsSensorSubsample, imagingPeakWavelength, micronsPerDegree)
    cmdStr = "diffractionLimitedOptics(" + ...
        "basePupilDiameterMm, wlsSensorSubsample, imagingPeakWavelength, " + ...
        "micronsPerDegree, 0, 'wavefrontSpatialSamples', 101, 'noLCA', true)";
    [~, theBaseOptics] = evalc(cmdStr);
end

function theConeMosaic = generateConeMosaicMuted(eccentricityDegs, coneMosaicSizeDegs, wlsSensorSubsample)
    cmdStr = "cMosaic(" + ...
        "'eccentricityDegs', eccentricityDegs, " + ...
        "'sizeDegs', coneMosaicSizeDegs*[1 1], " + ...
        "'wave', wlsSensorSubsample, " + ...
        "'sourceLatticeSizeDegs', 64, " + ...
        "'rodIntrusionAdjustedConeAperture', true, " + ...
        "'coneApertureModifiers', struct('smoothLocalVariations', true))";
    [~, theConeMosaic] = evalc(cmdStr);
end

function thePSF = generatePSFMuted(examinedDefocusAndAstigmatismDioptersGrid, theConeMosaic, localOpticsParams, psfUpsampleFactor, visualizeStrehlRatioOptimization, contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection, fovealOpticsDiopters, useParforInner)
    cmdStr = ['RGCMosaicConstructor.helper.optics.optimized3DStrehlRatioPSF(' ...
        'examinedDefocusAndAstigmatismDioptersGrid, ' ...
        'theConeMosaic, theConeMosaic.eccentricityDegs, localOpticsParams, ' ...
        'localOpticsParams.wavefrontSpatialSamples, psfUpsampleFactor, ' ...
        'visualizeStrehlRatioOptimization, contrastMaxStrehlRatioPSFtoAsMeasuredAndCentralCorrection, ' ...
        'fovealOpticsDiopters(1), fovealOpticsDiopters(2), fovealOpticsDiopters(3), ' ...
        'useParforInner)'];
    [~, ~, thePSF] = evalc(cmdStr);
end