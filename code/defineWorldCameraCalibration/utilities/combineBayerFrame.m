% Uses the specified bayerPattern to re-assemble a raw image from
% separate R, G, and B channels while preserving the input variable type.
function rawFrame = combineBayerFrame(R, G, B, bayerPattern)

    % Get dimensions of a single channel
    [rows, cols] = size(R);
    
    % Initialize the output frame using ones, which safely preserves 
    % data classes such as double, single, uint8, or uint16.
    rawFrame = ones(rows * 2, cols * 2, class(R));

    % Cast G to match the class of R just in case they differ
    G_typed = cast(G, class(R));

    switch upper(bayerPattern)

        case "BGGR"
            rawFrame(1:2:end, 1:2:end) = B;
            rawFrame(1:2:end, 2:2:end) = G_typed;
            rawFrame(2:2:end, 1:2:end) = G_typed;
            rawFrame(2:2:end, 2:2:end) = R;

        case "RGGB"
            rawFrame(1:2:end, 1:2:end) = R;
            rawFrame(1:2:end, 2:2:end) = G_typed;
            rawFrame(2:2:end, 1:2:end) = G_typed;
            rawFrame(2:2:end, 2:2:end) = B;

        case "GRBG"
            rawFrame(1:2:end, 1:2:end) = G_typed;
            rawFrame(1:2:end, 2:2:end) = R;
            rawFrame(2:2:end, 1:2:end) = B;
            rawFrame(2:2:end, 2:2:end) = G_typed;

        case "GBRG"
            rawFrame(1:2:end, 1:2:end) = G_typed;
            rawFrame(1:2:end, 2:2:end) = B;
            rawFrame(2:2:end, 1:2:end) = R;
            rawFrame(2:2:end, 2:2:end) = G_typed;

        otherwise
            error("Unknown Bayer pattern: %s", bayerPattern)

    end

end