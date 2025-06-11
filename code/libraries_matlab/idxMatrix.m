function idxMatrix = returnPixelIdx(chunks, targetChunk, pixelClass)
% Generates 1D pixel location mask for a specific and RBG value. 
%
% Syntax:
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
%                      (Case-insensitive).
%
% Output:
%   idxMatrix     Single 1D column vector (a mask).
%                 - '1' at indices where the specified color pixel is located, and '0' elsewhere.
%                 - Length is (nRows * nCols).
%
% Examples:
%{
    TBD
%}


% If the input 'chunks' is empty, return an empty cell array. 
if isempty(chunks)
    idxMatrix = {};
    return;
end

% Navigate to the raw video frames for specified chunk.
targetData = chunks{chunkTarget}.W.v;

% Extract the dimensions of a single frame from the 3D video data. Ignore
% number of frames.
[~, nRows, nCols] = size(targetData);

% Initialize 2D Bayer mask.
temp_idxMatrix = zeros(nRows, nCols, dtype=np.uint8)


switch pixelClass
    case 'R'
        for rr = 1:nRows; 
            for cc = 1:nCols;
                idxMatrix(rr,cc) = some calculation
        end
        end
    case 'G'
    case 'B'
end
end