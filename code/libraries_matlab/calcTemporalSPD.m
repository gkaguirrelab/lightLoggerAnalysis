function [spd,frq] = calcTemporalSPD( v, fps, options )
% Return the temporal spectral power density for a world camera video
%
% Syntax:
%   [spd,frq] = calcTemporalSPD( video, t )
%
% Description:
%   This routine accepts 640x480 "raw" format videos from the light logger.
%   The duration of the passed video can vary. The routine returns the
%   spectral power density of the video across the time domain. By default,
%   this will be the temporal SPD for the spatial average value for each
%   frame. Optional arguments adjust this behavior.
%
% Inputs:
%   v                     - tx480x640 array of 8 bit unsigned integers.
%                           This is a "chunk" of the world camera video
%   fps                   - Scalar. Sampling rate of the recording (frames
%                           per second). Defaults to 200.
%
% Optional key/value pairs:
%  'lineResolution'       - Logical. If set to true, the SPD is calculated
%                           using each vertical line of the sensor,
%                           providing for temporal resolution = fps * 480.
%  'byChannel'            - Logical. Return the SPD separately for the R,
%                           G, and B sensor channels.
% Outputs:
%   spd                   - 1xf. The spectral power density in units of
%                           contrast^2/Hz.
%   frq                   - 1xf. The frequencies (in Hz) for the spd values
%
% Examples:
%{
%}

arguments
    v (:,480,640) {mustBeNumeric}
    fps (1,1) {mustBeScalarOrEmpty} = 200
    options.applyFieldCorrection (1,1) logical = false
    options.lineResolution (1,1) logical = false
    options.spatialChannel (1,1) string = 'RGB';
end

% Combine data across space. We modify v to either nan components we do
% not wish to include, or replace elements with derived values
switch options.spatialChannel
    case 'RGB'
        % no change
    otherwise
        error('Not a defined channel setting')
end

% If we wish to treat each line of the video as a separate time point, then
% we reformat v to have the dimensions 640 x 1 x [t * 480]
if options.lineResolution
    % Need to implement this step properly
    % v = reshape(v',1,numel(v));    
    fps = fps * 480;
end

% Take the mean of each frame, accounting for nan values. This is the
% signal
signal = squeeze(mean(mean(v,2,"omitmissing"),3,"omitmissing"));

% Convert to contrast units
signal = (signal - mean(signal))/mean(signal);

% FFT of the signal in contrast units
[frq,amp] = simpleFFT(signal,fps);

% Discard the zero and nyquit frequencies
frq = frq(2:end-1); amp = amp(2:end-1);

% Covert from amplitude (contrast) to spectral power density
% (contrast^2/Hz)
spd=(amp.^2)./frq;

end

