function converted = rgb2lms(frame, T_receptors, T_camera, options)
% Convert an RGB camera frame into cone-response LMS coordinates.
%
% Syntax:
%   converted = rgb2lms(frame, T_receptors, T_camera, options)
%
% Description:
%   This function interprets each RGB pixel as weights on the camera's
%   spectral sensitivity basis, reconstructs the implied spectral signal,
%   and projects that signal onto human L, M, and S cone sensitivity
%   functions. The conversion is applied framewise in matrix form after
%   reshaping the image into a list of pixels, then the result is restored
%   to image layout. Inputs at or beyond the 0 and 255 bounds are rejected
%   because the downstream LMS analyses expect those values to have been
%   masked or excluded beforehand.
%
% Inputs:
%   frame                    - Numeric image array of size
%                              `[nRows, nCols, 3]` containing RGB camera
%                              values.
%   T_receptors              - Optional 3-by-N matrix of receptor
%                              sensitivities. If omitted, it is loaded
%                              from `generate_LMS_transformation_info`.
%   T_camera                 - Optional 3-by-N matrix describing the
%                              camera spectral sensitivities used to map
%                              RGB values into spectral space.
%   options                  - Name/value options. `camera` selects which
%                              built-in camera sensitivity model to use
%                              when the transform matrices are omitted.
%
% Outputs:
%   converted                - Numeric image array of size
%                              `[nRows, nCols, 3]` containing the L, M,
%                              and S responses for each pixel.
%
% Examples:
%{
    converted = rgb2lms(frame, [], [], "camera", "imx219");
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
