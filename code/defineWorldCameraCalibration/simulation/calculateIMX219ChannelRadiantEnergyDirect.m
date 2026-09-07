function channelEnergy = calculateIMX219ChannelRadiantEnergyDirect(cameraRadianceMap)
% CALCULATECHANNELRADIANTENERGY Computes the total radiant flux/energy 
% in each channel (R, G, B) from a camera radiance map, weighting each pixel 
% by its solid angle (steradians) and accounting for the Bayer pattern.

persistent bayerIdx pixelSolidAngles channelNames
if isempty(bayerIdx)
    bayerPattern = "BGGR";
    channelNames = {'red', 'green', 'blue'};
    [rows, cols] = size(cameraRadianceMap);
    
    % Obtain Bayer channel pixel indices
    [bayerIdx{1}, bayerIdx{2}, bayerIdx{3}] = returnBayerIndices(cameraRadianceMap, bayerPattern);
    
    % Load precomputed pixel solid angles from the derived file
    projectRoot = tbLocateProjectSilent('lightLoggerAnalysis');
    deltaPath = fullfile(projectRoot, 'derived', 'deltaSteradians.mat');
    assert(isfile(deltaPath), ...
        'calculateChannelRadiantEnergy:MissingDeltaSteradians', ...
        'Could not find derived deltaSteradians.mat file. Please run defineDeltaSteradians first.');
    
    load(deltaPath, 'deltaSteradians');
    pixelSolidAngles = deltaSteradians;
end

channelEnergy = struct();
for cc = 1:3
    % Isolate pixels belonging to this specific Bayer channel
    channelMask = false(size(cameraRadianceMap));
    channelMask(bayerIdx{cc}) = true;
    
    % Integrate radiance over the solid angle subtended by each channel's pixels
    channelPixels = cameraRadianceMap .* channelMask;
    channelEnergy.(channelNames{cc}) = sum(channelPixels(:) .* pixelSolidAngles(:));
end

end