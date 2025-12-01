function [spd, frq] = calcTemporalSPD(v, fps, options)
% Return the temporal spectral power density for a video
%
% Syntax:
%   [spd, frq] = calcTemporalSPD(video, t)
%t
% Description:
%   This routine accepts "raw" format videos from the light logger. The
%   duration of the passed video can vary. The routine returns the spectral
%   power density of the video across the time domain. By default, this
%   will be the temporal SPD for the spatial average value for each frame.
%   Optional arguments adjust this behavior.
%
% Inputs:
%   v                     - txVertxHoriz array of 8 bit unsigned integers.
%                           This is a "chunk" of the world camera video
%   fps                   - Scalar. Sampling rate of the recording (frames
%                           per second). Defaults to 120.
%
% Optional key/value pairs:
%  'lineResolution'       - Logical. If set to true, the SPD is calculated
%                           using each vertical line of the sensor,
%                           providing for temporal resolution = fps * 480.
%  'byChannel'            - Logical. Return the SPD separately for the R,
%                           G, and B sensor channels.
%
% Outputs:
%   spd                   - 1xf. The spectral power density in units of
%                           contrast^2/Hz.
%   frq                   - 1xf. The frequencies (in Hz) for the spd values
%
% Examples:
%{
    fps = 120;
    num_frames = 30 * fps;
    v = rand(num_frames, 480, 640)*255;
    [spd, frq] = calcTemporalSPD(v, fps);
    plotSPD(spd, frq);
%}

arguments
    v (:,480,480) {mustBeNumeric}
    fps (1,1) {mustBeScalarOrEmpty} = 120
    options.lineResolution (1,1) logical = true
    options.regionMatrix (:,:) {mustBeNumeric} = ones(480,480);
    options.applyFieldCorrection (1,1) logical = false
    options.byChannel (1,1) logical = false
    options.postreceptoralChannel {mustBeMember(options.postreceptoralChannel,{'LM', 'L-M', 'S'})} = 'LM'
    options.camera (1,:) char {mustBeMember(options.camera, {'standard', 'imx219'})} = 'imx219'
    options.nan_threshold = 0.1;
end

% Convert video data to double and get dimensions
v = double(v);
[nFrames, nRows, nCols] = size(v);

% Reshape 3D video into 2D matrix
if options.lineResolution
    % Sum across columns and Transpose to get nFrames x nRows matrix
    v_RowsFrames = squeeze(sum(v, 3))';

    % Reshape into a vector with length nFrames * nRows
    signal = v_RowsFrames(:);

    % Convert to contrast units relative to mean
    signal = (signal - mean(signal)) / mean(signal);

    % Modify fps to reflect rolling shutter temporal resolution
    fps = fps * nRows;

else
    % Ensure we are selecting the right mask
    mask = options.regionMatrix ~= 0;

    % Flatten spatial dims into a single pixel dimension
    v_flat = reshape(v, [nFrames, nRows * nCols]);      % nFrames x (nRows*nCols)

    % Select only pixels in the region
    regionPixels = v_flat(:, mask(:));

    % Mean the pixels per image per frame to achieve the signal
    signal = mean(regionPixels, 2, 'omitmissing');

    % Convert to contrast units
    signal = (signal - mean(signal, 'omitmissing')) / mean(signal, 'omitmissing');

end

if(numel(signal(:)) == 0)
    error("Signal is empty");
end

% Find the percent nan of the signal. If it is large, simply return nan now
if( numel(find(isnan(signal(:)))) / numel(signal(:)) >= options.nan_threshold )
    frq = nan; spd = nan;
    return;
end

% Otherwise, obtain the PSD of the signal
[frq, spd] = simplePSD(double(signal), fps);


end

