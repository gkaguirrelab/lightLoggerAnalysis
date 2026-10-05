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

% Definitions
Smin = darkSignal;
Smax = 2^8-1 - Smin;

% First steps
y = double(y) + 0.375;
yPrime = y - darkSignal;

% Use a non-negative version strictly for the asymptotic gain nonlinearity 
% to avoid complex numbers from fractional powers of negative numbers, 
% while allowing yPrime to retain signed values.
yPrimeNonlinear = max(0, yPrime);
asymptoticGain = 1 ./ (1 - (yPrimeNonlinear ./ Smax).^n).^(1./n);

% Apply asymptotic gain to the signed linear signal
yLinear = yPrime .* asymptoticGain;


end