function inspectGeneratedData()

    whichMeridian = 'Temporal';
    rootDir = '/Volumes/SSDdisk/MATLAB/projects/lightLoggerAnalysis/derived';
    load(fullfile(rootDir, sprintf('radianceToConeRateSupportCMosaicBased%sMeridianWithOptimizedOptics.mat', whichMeridian)), ...
        'coneMapVar', 'opticsSupport');

    age = 18;

    figure(33);
    coneMapVar(age).eccGrid
    for iEcc = 1:numel(coneMapVar(age).eccGrid)
        thePSF = opticsSupport.thePSFs{iEcc};
        subplot(4,5,iEcc)
        thePSF.supportWavelength
        pause
        imagesc(squeeze(thePSF.data(:,:,100)))
    end
end
