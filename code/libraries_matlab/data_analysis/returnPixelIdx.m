function [idxMatrix_flat, idxMatrix_mat] = returnPixelIdx(pixelClass, options)
% Generates 1D pixel location mask for a specific and RBG value. 
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
%   idxMatrix     Single 1D column vector (a mask).
%                 - '1' at indices where the specified color pixel is located, and '0' elsewhere.
%                 - Length is (nRows * nCols).
%
% Examples:
%{
    ---- To Use With Chunks Data
    1. ** create necessary file paths **
    2. ** parse Chunks as normal **

    [~, actual_nRows, actual_nCols] = size(chunks{i}.W.v)
    [X_mask_flat, X_mask_mat] = returnPixelIdx('X', 'nRows', actual_nRows, 'nCols', actual_nCols);
    
    ---- To Express Mean Signal of Each Frame for Different Color Channel
    myChunk = chunks{1};
    for tt=1:5400; 
    myFrame = squeeze(myChunk.W.v(tt,:,:)); 
    myVal(tt) = mean(myFrame(myIdx)); 

%}

arguments
    pixelClass (1,1) string {mustBeMember(pixelClass, {'R', 'G', 'B'})}
    options.nRows (1,1) {mustBePositive, mustBeInteger} = 480
    options.nCols (1,1) {mustBePositive, mustBeInteger} = 640
end

nRows = options.nRows;
nCols = options.nCols;

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
    % Mark positions for BLUE pixels 
    case 'B'
        for rr = 1:nRows
            for cc = 1:nCols
                if mod(rr - 1, 2) == 0 && mod(cc - 1, 2) == 0
                    temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
    otherwise
            error("Invalid pixelClass. Use 'R', 'G', or 'B'.");
end

% Flatten the 2D mask into a 1D column vector
idxMatrix_flat = temp_idxMatrix(:);
% Preserve 2D mask 
idxMatrix_mat = temp_idxMatrix;
