function idxMatrix = returnPixelIdx(chunks, targetChunk, pixelClass)
% Generates 1D pixel location mask for a specific and RBG value. 
%
% Syntax:l
%   idxMatrix = returnPixelIdx(chunks, targetChunk, pixelClass)
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
%   chunks             Cell array of structs (i.e., output from parse_chunks).
%                      Each struct is expected to contain a 'W' field, which
%                      in turn has a 'v' field with dimensions [num_frames x rows x cols].
%   targetChunk        Scalar integer. The 1-based index of the specific chunk
%                      within 'chunks' for which to determine mask dimensions.
%   pixelClass         String ('R', 'G', or 'B') specifying the color to identify.
%
% Output:
%   idxMatrix     Single 1D column vector (a mask).
%                 - '1' at indices where the specified color pixel is located, and '0' elsewhere.
%                 - Length is (nRows * nCols).
%
% Examples:
%{
    chunks = cell(3, 1);
    chunks{1}.W.v = uint8(rand(30, 480, 640) * 255);
    chunks{2}.W.v = uint8(rand(30, 480, 640) * 255);
    chunks{3}.W.v = uint8(rand(30, 480, 640) * 255);
    idxMatrix_RED_chunk1 = returnPixelIdx(chunks, 1, 'R');
    idxMatrix_GREEN_chunk2 = returnPixelIdx(chunks, 2, 'G');
    idxMatrix_BLUE_chunk3 = returnPixelIdx(chunks, 3, 'B');
%}


% If the input 'chunks' is empty, return an empty cell array. 
if isempty(chunks)
    idxMatrix = {};
    return;
end

% Navigate to the raw video frames for specified chunk.
targetData = chunks{targetChunk}.W.v;

% Extract the dimensions of a single frame from the 3D video data. Ignore
% number of frames.
[~, nRows, nCols] = size(targetData);

disp(['DEBUG: Inside returnPixelIdx - Derived nRows: ', num2str(nRows)]);
disp(['DEBUG: Inside returnPixelIdx - Derived nCols: ', num2str(nCols)])

% Initialize 2D Bayer mask.
temp_idxMatrix = zeros(nRows, nCols, 'uint8');

switch pixelClass
    % Mark positions for RED pixels. 
    case 'R'
        for rr = 1:nRows
            for cc = 1:nCols
                if mod(rr, 2) ~= 0 && mod(cc, 2) ~= 0
                   temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
    % Mark positions for GREEN pixels. 
    case 'G'
        for rr = 1:nRows
            for cc = 1:nCols
                if (mod(rr, 2) == 0 && mod(cc, 2) ~= 0) || ...
                   (mod(rr, 2) ~= 0 && mod(cc, 2) == 0)
                    temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
    % Mark positions for BLUE pixels. 
    case 'B'
        for rr = 1:nRows
            for cc = 1:nCols
                if mod(rr, 2) == 0 && mod(cc, 2) == 0
                    temp_idxMatrix(rr, cc) = 1;
                end
            end
        end
    otherwise
            error("Invalid pixelClass. Use 'R', 'G', or 'B'.");
end

disp(['DEBUG: Inside returnPixelIdx - Size of tempBayerMask2D before flattening: ', num2str(size(temp_idxMatrix))]);

% Flatten the 2D mask into a 1D column vector. 
idxMatrix = temp_idxMatrix;