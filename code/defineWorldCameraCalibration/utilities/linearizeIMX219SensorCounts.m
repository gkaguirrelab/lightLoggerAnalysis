function yLinear = linearizeIMX219SensorCounts(y, n)
% This function implements a correction of the soft saturating non-
% linearity that is explored in the routine defineFullWellCapacityEffect.
% The function takes raw sensor counts from the IMX219 chip, along with the
% empirically measured exponent parameter, and returns the linearized
% value. The routine has hard-coded the minimum dark value and the
% bit-depth of the camera. The linearized values are returned as floats.

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

% Calculate the unclipped, linear portion (a * x)
ax = yPrime ./ (1 - (yPrime ./ Smax).^n).^(1./n);

% Calculate the maximum possible value of ax
yMax = Smax - 1;
aMax = yMax ./ (1 - (yMax ./ Smax).^n).^(1./n);

% We scale the linearized y so that it maps the 0-254 original range to
% 0-254 on the output side. Any input values that were larger than 254 get
% mapped to 255.
yLinear = (ax ./ aMax) .* (2^8-2);
yLinear(isinf(yLinear)) = (2^8-1);

end
