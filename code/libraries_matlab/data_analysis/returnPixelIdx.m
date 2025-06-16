function [idxMatrix_flat, idxMatrix_mat] = returnPixelIdx(pixelClass, options)
% Generates 1D and 2D pixel location masks for a specific RBG value. 
%
% Syntax:l
%   idxMatrix = returnPixelIdx(pixelClass, options)
% 
% Description;
%   Calculates the Bayer filter mask (1 for specified color, 0 otherwise)
%   for the frame dimensions found in the specified chunk.
%
%   This mask represents the physical locations of the R, G, or B pixels
%   on the camera sensor. This mask is static and does not change based on
%   the content of the video frame, only its dimensions.
% 
% Inputs:
%   pixelClass          - String ('R', 'G', or 'B') specifying the color to identify.
%
% Optional key/value pairs:
%   'nRows'             - Scalar integer. The number of rows (height) of the camera frame.
%                         Defaults to 480 if not provided.
% 
%   'nCols'             - Scalar integer. The number of columns (width) of the camera frame.
%                         Defaults to 640 if not provided.
%
% Output:
%   idxMatrix_flat      1D column vector (mask).
%                       - '1' at indices where the specified color pixel is located, and '0' elsewhere.
%                       - Length is (nRows * nCols).
% 
%   idxMatrix_mat       2D matrix (mask).
%                       - '1' at indices where the specified color pixel is located, and '0' elsewhere.
%                       - Dimensions are (nRows * nCols).
%
% Examples:
%{
    ---- To Use With Chunks Data ----
    1. ** create necessary file paths **
    2. ** parse Chunks as normal **

    [~, actual_nRows, actual_nCols] = size(chunks{i}.W.v)
    [X_mask_flat, X_mask_mat] = returnPixelIdx('X', 'nRows', actual_nRows, 'nCols', actual_nCols);
    
    OR (recommended)

    Use plotMeanPixel.m script to express and visualize mean signal of each frame in
    the chunk for all three color channels.
%}

arguments
    pixelClass (1,1) string {mustBeMember(pixelClass, {'R', 'G', 'B'})}
    options.nRows (1,1) {mustBePositive, mustBeInteger} = 480
    options.nCols (1,1) {mustBePositive, mustBeInteger} = 640
end

nRows = options.nRows;
nCols = options.nCols;

assert(mod(nRows, 2) == 0 && mod(nCols, 2) == 0, ...
    'Assertion failed: Frame dimensions must be even for consistent Bayer pattern generation.');

% Initialize 2D Bayer mask
temp_idxMatrix = zeros(nRows, nCols, 'uint8');

switch pixelClass
    % Mark positions for RED pixels
    case 'R'
        for rr = 1:nRows
            for cc = 1:nCols
                if mod(rr - 1, 2) ~= 0 && mod(cc - 1, 2) ~= 0
                   temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
        assert((sum(temp_idxMatrix, 'all')) == (nRows * nCols / 2)), ...
            'Assertion failed: Incorrect count of active pixels (1s) for RED channel';
    % Mark positions for GREEN pixels 
    case 'G'
        for rr = 1:nRows
            for cc = 1:nCols
                if (mod(rr - 1, 2) == 0 && mod(cc - 1, 2) ~= 0) || ...
                   (mod(rr - 1, 2) ~= 0 && mod(cc - 1, 2) == 0)
                    temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
        assert((sum(temp_idxMatrix, 'all')) == (nRows * nCols / 4)), ...
            'Assertion failed: Incorrect count of active pixels (1s) for RED channel';
    % Mark positions for BLUE pixels 
    case 'B'
        for rr = 1:nRows
            for cc = 1:nCols
                if mod(rr - 1, 2) == 0 && mod(cc - 1, 2) == 0
                    temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
        assert((sum(temp_idxMatrix, 'all')) == (nRows * nCols / 2)), ...
            'Assertion failed: Incorrect count of active pixels (1s) for BLUE channel';
    otherwise
            error("Invalid pixelClass. Use 'R', 'G', or 'B'.");
end

% Flatten the 2D mask into a 1D column vector
idxMatrix_flat = temp_idxMatrix(:);
% Preserve 2D mask 
idxMatrix_mat = temp_idxMatrix;
