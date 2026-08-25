function [R, G, B] = returnBayerIndices(rawFrame, bayerPattern)

% Create a matrix of linear indices the same size as the raw frame
idx = reshape(1:numel(rawFrame), size(rawFrame));

switch upper(bayerPattern)

    case "BGGR"
        B  = idx(1:2:end, 1:2:end);
        G1 = idx(1:2:end, 2:2:end);
        G2 = idx(2:2:end, 1:2:end);
        R  = idx(2:2:end, 2:2:end);

    case "RGGB"
        R  = idx(1:2:end, 1:2:end);
        G1 = idx(1:2:end, 2:2:end);
        G2 = idx(2:2:end, 1:2:end);
        B  = idx(2:2:end, 2:2:end);

    case "GRBG"
        G1 = idx(1:2:end, 1:2:end);
        R  = idx(1:2:end, 2:2:end);
        B  = idx(2:2:end, 1:2:end);
        G2 = idx(2:2:end, 2:2:end);

    case "GBRG"
        G1 = idx(1:2:end, 1:2:end);
        B  = idx(1:2:end, 2:2:end);
        R  = idx(2:2:end, 1:2:end);
        G2 = idx(2:2:end, 2:2:end);

    otherwise
        error("Unknown Bayer pattern: %s", bayerPattern)

end

% Convert to linear index column vectors.
R = R(:);
G = [G1(:); G2(:)];
B = B(:);

end