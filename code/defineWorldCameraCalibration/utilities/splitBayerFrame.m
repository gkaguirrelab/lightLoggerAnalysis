% Uses the specified bayerPattern to return the image partitioned by
% channel
function [R, G, B] = splitBayerFrame(rawFrame, bayerPattern)

switch upper(bayerPattern)

    case "BGGR"
        B  = rawFrame(1:2:end, 1:2:end);
        G1 = rawFrame(1:2:end, 2:2:end);
        G2 = rawFrame(2:2:end, 1:2:end);
        R  = rawFrame(2:2:end, 2:2:end);

    case "RGGB"
        R  = rawFrame(1:2:end, 1:2:end);
        G1 = rawFrame(1:2:end, 2:2:end);
        G2 = rawFrame(2:2:end, 1:2:end);
        B  = rawFrame(2:2:end, 2:2:end);

    case "GRBG"
        G1 = rawFrame(1:2:end, 1:2:end);
        R  = rawFrame(1:2:end, 2:2:end);
        B  = rawFrame(2:2:end, 1:2:end);
        G2 = rawFrame(2:2:end, 2:2:end);

    case "GBRG"
        G1 = rawFrame(1:2:end, 1:2:end);
        B  = rawFrame(1:2:end, 2:2:end);
        R  = rawFrame(2:2:end, 1:2:end);
        G2 = rawFrame(2:2:end, 2:2:end);

    otherwise
        error("Unknown Bayer pattern: %s", bayerPattern)

end

G = (G1 + G2) ./ 2;

end