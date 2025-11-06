function [r2,amplitude,phase,fit] = fourierRegression( signal, tSecs, fHz )
% Fourier regression using a sin and cosine at the specified frequency
%
% Syntax:
%   [r2,amplitude,phase,fit] = fourierRegression( signal, tSecs, fHz )
%
% Description:
%   This function performs a least-squares regression of the input signal 
%   onto sine and cosine basis functions at the specified frequency fHz. It
%   computes the amplitude and phase of the best-fit sinusoid, reconstructs 
%   the fitted waveform, and returns the coefficient of determination (r2) 
%   as a measure of accuracy.
% Inputs:
%   signal                - Numeric vector. The measured signal values.
%   tSecs                 - Numeric vector. Time stamps in seconds 
%                           corresponding to each sample in the signal.
%   fHz                   - Scalar. The target frequency in Hertz at which
%                           to fit the sinusoid.
%
% Outputs:
%   r2                    - Scalar. The coefficient (0-1) of determination 
%                           between the fitted sinusoid and the original
%                           signal.
%   amplitude             - Scalar. The amplitude of the fitted sinusoid.  
%   phase                 - Scalar. The phase offset of the fitted sinusoid
%                           in radians(?)
%   fit                   - Numeric vector. The reconstructed sinusoidal
%                           fit values aligned at tSecs.
%
% Examples:
%{
    [r2, amplitude, phase, fit] = fourierRegression(signal, tSecs, 10)
    
    figure
    plot(signal, '.r')
    hold on
    plot(fit,'b')
%}


arguments
    signal (1,:) {mustBeFloat}
    tSecs (1,:) {mustBeNumeric}
    fHz (1,1) {mustBeNumeric}
end

% Set up the regression matrix
X = [];
X(:,1) = sin(  tSecs./(1/fHz).*2*pi );
X(:,2) = cos(  tSecs./(1/fHz).*2*pi );

% Perform the fit
b = X\signal';

% Derive some results
fit = (X * b)';

% Adjust the mean offset of the fit to match the signal
fit = fit - mean(fit) + mean(signal);

% Calculate and return the amplitude, phase, and r2
amplitude  = norm(b);
phase = -atan(b(2)/b(1));
r2 = corr(fit',signal')^2;

end