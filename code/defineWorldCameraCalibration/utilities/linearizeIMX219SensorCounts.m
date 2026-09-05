function yLinear = linearizeIMX219SensorCounts(y, n)
% This function implements a correction of the soft saturating non-
% linearity that is explored in the routine defineFullWellCapacityEffect.
% The function takes raw sensor counts from the IMX219 chip, along with the
% empirically measured exponent parameter, and returns the unbounded 
% linearized value (representing the true photon-equivalent portion).

% The "dark signal" value is the empirically measured sensor value reported
% by the camera under conditions of zero true photon capture. We load this
% value into a persistent variable.
persistent darkSignal
if isempty(darkSignal)
    paramFileName = fullfile(...
        tbLocateProjectSilent('lightLoggerAnalysis'),...
        'derived',...
        'darkSignal.mat');
    load(paramFileName,'darkSignal');
end
Smin = darkSignal;

% Define the fixed asymptotic value
Smax = 2^8-1 - Smin;

% Force all image pixels to above the min value
y(y<Smin) = Smin;

% Subtract the baseline offset
yPrime = y - Smin;

% Calculate and return the unclipped, unbounded linear portion (a * x)
% Note: Input values corresponding precisely to maximum full-well capacity 
% will natively evaluate to Inf, preserving physical magnitude assumptions 
% for saturated pixels.
yLinear = yPrime ./ (1 - (yPrime ./ Smax).^n).^(1./n);

end