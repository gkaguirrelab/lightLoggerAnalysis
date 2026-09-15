function cameraScoreHistory = simulateAGCRadiance(radianceLevels, initialCameraScore, speedSetting)
% SIMULATEAGCRADIANCE Simulates adaptive gain control behavior.
%
% Inputs:
%   radianceLevels      - 1D Array. Vector of absolute light levels (radiance) for each frame.
%   initialCameraScore  - Double. Starting combined camera sensitivity score.
%   speedSetting        - Double. The base speed setting parameter. Defaults to 0.95.
%
% Outputs:
%   cameraScoreHistory  - 1D Array. Combined camera sensitivity score at each frame.

    % Handle default arguments
    if nargin < 3 || isempty(speedSetting)
        speedSetting = 0.95;
    end

    % Define hardware bounds and operational targets
    signalTarget = 127;
    exposureRange = [37, 8333];
    aGainRange = [1, 10.666];
    dGainRange = [1, 10]; 

    % Decompose the initial camera score into separate AGC settings
    % using the provided cameraScoreToAGCSettings function.
    AGCSettings = cameraScoreToAGCSettings(initialCameraScore);
    
    currentExposure = AGCSettings.exposure; %[cite: 4]
    currentAGain = AGCSettings.Again;       %[cite: 4]
    currentDGain = AGCSettings.Dgain;       %[cite: 4]

    % Initialize history tracking array
    numFrames = length(radianceLevels);
    cameraScoreHistory = zeros(1, numFrames);

    % --- Simulation Loop ---
    for k = 1:numFrames
        % 1. Record current sensitivity score
        thisCameraScore = currentAGain * currentDGain * currentExposure;
        cameraScoreHistory(k) = thisCameraScore;

        % 2. Calculate observed signal with sensor clipping limits
        sRaw = radianceLevels(k) * thisCameraScore; 
        sK = max(0, min(255, sRaw)); 

        % 3. Calculate initial correction factor
        correction = 1 + (signalTarget - sK) / signalTarget;

        % 4. Determine dynamic speed scaling
        speed = speedSetting;
        if sK == 0 || sK == 255
            speed = speedSetting^3;
        elseif abs(correction - 1) < 0.25
            speed = speedSetting^2;
        end

        % 5. Apply speed smoothing
        correction = 1 + ((1 - speed) * (correction - 1));

        % 6. Apply Hardware Clamping and Priority Rules (Exposure -> AGain -> DGain)
        if correction > 1
            if currentExposure < exposureRange(2)
                currentExposure = max(exposureRange(1), min(currentExposure * correction, exposureRange(2)));
            elseif currentAGain < aGainRange(2)
                currentAGain = max(aGainRange(1), min(currentAGain * correction, aGainRange(2)));
            else
                currentDGain = max(dGainRange(1), min(currentDGain * correction, dGainRange(2)));
            end
        elseif correction < 1
            if currentDGain > dGainRange(1)
                currentDGain = max(dGainRange(1), min(currentDGain * correction, dGainRange(2)));
            elseif currentAGain > aGainRange(1)
                currentAGain = max(aGainRange(1), min(currentAGain * correction, aGainRange(2)));
            else
                currentExposure = max(exposureRange(1), min(currentExposure * correction, exposureRange(2)));
            end
        end
    end

end