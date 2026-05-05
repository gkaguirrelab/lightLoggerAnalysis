function converted = rgb2lms(frame, T_receptors, T_camera, options)
% 
% 
% Description: 
%   TODO 
% Inputs: 
%   TODO
% 
% Output:
%   TODO
% 
% Table that contains the average spectral sensitivity functions of the
% camera sensors in commercial mobile phones, reported in:
%   Tominaga S, Nishi S, Ohtera R. Measurement and estimation of spectral
%   sensitivity functions for mobile phone cameras. Sensors. 2021 Jul
%   22;21(15):4985.
% 
% Examples:
%{
    TODO
%}

    arguments 
        frame 
        T_receptors = []
        T_camera = []
        options.camera {mustBeMember(options.camera, ["imx219", "standard"])} = "imx219"
    end 

    % Input cannot be 0s or 255s as we do not have a valid conversion 
    % for these values to LMS
    assert(~any(frame(:) <= 0));
    assert(~any(frame(:) >= 255)); 

    if(isempty(T_receptors) || isempty(T_camera))
        [T_receptors, T_camera] = generate_LMS_transformation_info(options.camera); 
    end

    % Let's get the size of the input frame 
    [nRows, nCols, nChannels] = size(frame);
    assert(nChannels == 3, 'Input frame must have size nRows x nCols x 3.');

    % Create the spectrum implied by the rgb camera weights, and then project
    % that on the receptors. ensure convert to double for floating point math 
    frame_flat = double(reshape(frame, [], 3));     
    converted_flat = (T_receptors * (frame_flat * T_camera )' )';

    % Then we need to reshape back to the shape of the image
    converted = reshape(converted_flat, nRows, nCols, 3);
    end 