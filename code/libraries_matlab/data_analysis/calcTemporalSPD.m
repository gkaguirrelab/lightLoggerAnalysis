function [spd,frq] = calcTemporalSPD( v, fps, options )
% Return the temporal spectral power density for a world camera video
%
% Syntax:
%   [spd,frq] = calcTemporalSPD( video, t )
%t
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
%                           per second). Defaults to 120.
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
    fps = 120;
    num_frames = 30 * fps;
    v = rand(num_frames, 480, 640)*255;
    [spd, frq] = calcTemporalSPD(v, fps);
    plotSPD(spd, frq);
%}

arguments
    v (:,480,640) {mustBeNumeric}
    fps (1,1) {mustBeScalarOrEmpty} = 120
    options.lineResolution (1,1) logical = true
    options.regionMatrix (:,:) {mustBeNumeric} = ones(480,640);
    options.applyFieldCorrection (1,1) logical = false
    options.byChannel (1,1) logical = false
    options.postreceptoralChannel {mustBeMember(options.postreceptoralChannel,{'LM', 'L-M', 'S'})} = 'LM'
    options.camera (1,:) char {mustBeMember(options.camera, {'standard', 'imx219'})} = 'imx219'
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
    mat2D_V = reshape(v, nFrames, nRows * nCols);

    % Prepare output
    channels = {'R','G','B'};
    rgbSignal = nan(nFrames, numel(channels));
    
    for cc = 1:numel(channels)
    
        % Get 2D Bayer Mask for this channel
        [~, mask2D] = returnPixelIdx(channels{cc}, 'nRows', nRows, 'nCols', nCols);
    
        % Flatten and apply regionMatrix mask with logical 'AND'
        maskVec = mask2D(:) & options.regionMatrix(:);
    
        % Compute average per frame
        rgbSignal(:, cc) = mean(mat2D_V(:, maskVec), 2);
    
    end
    
    % Convert RGB --> LMS contrast relative to background
    background_RGB = mean(rgbSignal, 1);
    modulation_RGB = rgbSignal - background_RGB;
    modulation_LMS = cameraToCones(modulation_RGB, options.camera);
    background_LMS = cameraToCones(background_RGB, options.camera);
    
    lmsSignal = modulation_LMS ./ background_LMS;
    
    % Select a post-receptoral channel
    switch options.postreceptoralChannel
        case {'LM'}
            signal = (lmsSignal(:,1)+lmsSignal(:,2))/2;
        case {'L-M'}
            signal = (lmsSignal(:,1)-lmsSignal(:,2));
        case {'S'}
            signal = ((lmsSignal(:,3)-lmsSignal(:,1))+lmsSignal(:,2))/2;
    end

end

% PSD of the signal
[frq, spd] = simplePSD(signal, fps);

end

