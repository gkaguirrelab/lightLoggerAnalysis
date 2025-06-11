% SCRIPT TO EXPRESS MEAN SIGNAL OF EACH FRAME IN A CHUNK FOR DIFFERENT COLOR CHANNELS

% --- Ensure 'returnPixelIdx.m' file is in MATLAB path or current folder
% --- Ensure 'chunks' variable is loaded in your workspace
%     (from parse_chunks or saved mat file)


% Define a chunk to analyze. Adjust this index value as needed
myChunkIdx = 1;
myChunk = chunks{myChunkIdx};

% Get dimensions from video data
[numFrames, actual_nRows, actual_nCols] = size(myChunk.W.v);

mask_R = returnPixelIdx('R', 'nRows', actual_nRows, 'nCols', actual_nCols);
mask_G = returnPixelIdx('G', 'nRows', actual_nRows, 'nCols', actual_nCols);
mask_B = returnPixelIdx('B', 'nRows', actual_nRows, 'nCols', actual_nCols);

for tt = 1:numFrames
    % Extract current 2D frame
    myFrame = squeeze(myChunk.W.v(tt,:,:));
    
    % Flatten 2D frame into a 1D column vector
    myFlatFrame = myFrame(:);

    % Apply color mask
    pixels_R = myFlatFrame(mask_R == 1);
    pixels_G = myFlatFrame(mask_G == 1);
    pixels_B = myFlatFrame(mask_B == 1);

    % Calculate mean intensity of R/G/B color pixels
    avgOverTime_R(tt) = mean(pixels_R);
    avgOverTime_G(tt) = mean(pixels_G);
    avgOverTime_B(tt) = mean(pixels_B);

end

% Visualize the result
figure;
plot(1:numFrames, avgOverTime_R, 'R', 1:numFrames, avgOverTime_G, 'G', 1:numFrames, avgOverTime_B, 'B');
xlabel('Frame Number');
ylabel('Mean Intensity');
title('Mean Signal of Each Color Channel Over Time');
legend('Red Channel', 'Green Channel', 'Blue Channel');
