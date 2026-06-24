function correctedFrame = applyFieldingFunction(inputFrame, flatField)

% applyFieldingFunction corrects one image/frame using a flat-field map.
%
% inputFrame = frame to correct
% flatField  = normalized flat-field image, same size as inputFrame
% correctedFrame = brightness-corrected frame

inputFrame = double(inputFrame);
flatField = double(flatField);

% Normalize flat field to peak = 1, in case it is not already
flatFieldNorm = flatField ./ max(flatField(:), [], 'omitnan');

% Avoid division problems
flatFieldNorm(flatFieldNorm <= 0) = NaN;

% Correct image
correctedFrame = inputFrame ./ flatFieldNorm;
end 