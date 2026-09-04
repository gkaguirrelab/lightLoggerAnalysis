function [idxR, idxG, idxB] = returnBayerIndices(rawFrame, bayerPattern)

% Create a matrix of linear indices the same size as the raw frame
idx = reshape(1:numel(rawFrame), size(rawFrame));

switch upper(bayerPattern)

    case "BGGR"
        idxB  = idx(1:2:end, 1:2:end);
        idxG1 = idx(1:2:end, 2:2:end);
        idxG2 = idx(2:2:end, 1:2:end);
        idxR  = idx(2:2:end, 2:2:end);

    case "RGGB"
        idxR  = idx(1:2:end, 1:2:end);
        idxG1 = idx(1:2:end, 2:2:end);
        idxG2 = idx(2:2:end, 1:2:end);
        idxB  = idx(2:2:end, 2:2:end);

    case "GRBG"
        idxG1 = idx(1:2:end, 1:2:end);
        idxR  = idx(1:2:end, 2:2:end);
        idxB  = idx(2:2:end, 1:2:end);
        idxG2 = idx(2:2:end, 2:2:end);

    case "GBRG"
        idxG1 = idx(1:2:end, 1:2:end);
        idxB  = idx(1:2:end, 2:2:end);
        idxR  = idx(2:2:end, 1:2:end);
        idxG2 = idx(2:2:end, 2:2:end);

    otherwise
        error("Unknown Bayer pattern: %s", bayerPattern)

end

% Convert to linear index column vectors.
idxR = idxR(:);
idxG = [idxG1(:); idxG2(:)];
idxB = idxB(:);

end