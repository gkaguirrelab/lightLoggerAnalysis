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
    fps = 200;
    num_frames = 30 * fps;
    v = rand(num_frames, 480, 640)*255;
    [spd, frq] = calcTemporalSPD(v, fps);
    plotSPD(spd, frq);
%}

arguments
    v (:,480,640) {mustBeNumeric}
    fps (1,1) {mustBeScalarOrEmpty} = 200
    options.applyFieldCorrection (1,1) logical = false
    options.lineResolution (1,1) logical = false
    options.spatialChannel (1,1) string = {'R', 'B', 'G'};
end

% Get dimensions from the input video data (v)
[nFrames, nRows, nCols] = size(v);

% Combine data across space. We modify v to either nan components we do
% not wish to include, or replace elements with derived values
switch options.spatialChannel
    case {'R', 'B', 'G'}
        % Create 2D binary mask across all frames where 1 is a pixel of the
        % desired color, everything else is 0
        [~, idxMatrix_mat_2D] = returnPixelIdx(options.spatialChannel, 'nRows', nRows, 'nCols', nCols);
        
        % Create a 3D binary mask by replicating 2D mask across the time
        % dimension. Ensures mask and 'v' are the same size.
        idxMatrix_mat_3D = repmat(idxMatrix_mat_2D, [nFrames, 1, 1]);
        
        % Apply v-compatible mask and set '0' pixel values to NaN
        v(idxMatrix_mat_3D ~= 1) = NaN;

        % ASSERT FUNCTION FOR NaN COUNT     ----(might not need this)
        assert(sum(isnan(v), 'all') == sum(idxMatrix_mat_3D == 0, 'all'), ...
        'Assertion failed: Number of NaNs in v (%d) does not match expected masked pixels (%d) for %s channel.');
        % ASSERT FUNCTION FOR NaN VALUES
        assert(all(isnan(v(idxMatrix_mat_3D == 0))) && (all(~isnan(v(idxMatrix_mat_3D ~= 0)))), ...
            'Assertion failed: Incorrect NaN pattern after spatial channel masking for selected channel.');

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

% PSD of the signal in units of contrast^2/Hz
[frq, spd] = simplePSD(signal, fps);

end

