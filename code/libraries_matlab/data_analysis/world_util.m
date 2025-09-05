% Define constants
WORLD_CAM_FPS = 120;  % Capture at 120 FPS
WORLD_FRAME_SHAPE = [480, 640];  % (rows, cols)
WORLD_FRAME_DTYPE = 'uint8';
WORLD_METADATA_DTYPE = 'double';
WORLD_USE_AGC = true;
WORLD_SAVE_AGC_METADATA = true;
WORLD_AGC_METADATA_COLS = {'Again', 'Dgain', 'exposure'};
WORLD_AGC_CHANGE_INTERVAL = 0.250;
WORLD_AGC_SPEED_SETTING = 0.95;

% Placeholder for AGC settings (replace with actual implementation)
WORLD_AGC_SETTINGS_RANGES = retrieve_settings_ranges('W'); 
WORLD_AGC_DISCRETE_STATES = retrieve_discrete_states('W');

% Fixed NDF settings
WORLD_NDF_LEVEL_SETTINGS = containers.Map('KeyType', 'int32', 'ValueType', 'any');
WORLD_NDF_LEVEL_SETTINGS(0) = [1.0, 1.0, 747];
WORLD_NDF_LEVEL_SETTINGS(1) = [1.0, 1.0, 7085];
WORLD_NDF_LEVEL_SETTINGS(2) = [10.666, 1.14, 8333];
WORLD_NDF_LEVEL_SETTINGS(3) = [10.666, 4.297, 8333];
WORLD_NDF_LEVEL_SETTINGS(4) = [10.666, 10.0, 8333];
WORLD_NDF_LEVEL_SETTINGS(5) = [10.666, 10.0, 8333];
WORLD_NDF_LEVEL_SETTINGS(6) = [10.666, 10.0, 8333];

% RGB scalar multipliers
WORLD_RGB_SCALARS = [1, 0.8223775, 0.95367937];

% Create empty mask
WORLD_RGB_MASK = zeros(WORLD_FRAME_SHAPE, 'uint8');

% Pixel color assignment (MATLAB is 1-indexed)
[rows, cols] = ndgrid(1:WORLD_FRAME_SHAPE(1), 1:WORLD_FRAME_SHAPE(2));

R_mask = mod(rows, 2) == 1 & mod(cols, 2) == 1;
G_mask = (mod(rows, 2) == 0 & mod(cols, 2) == 1) | ...
         (mod(rows, 2) == 1 & mod(cols, 2) == 0);
B_mask = mod(rows, 2) == 0 & mod(cols, 2) == 0;

WORLD_RGB_MASK(R_mask) = 1;  % R=1
WORLD_RGB_MASK(G_mask) = 2;  % G=2
WORLD_RGB_MASK(B_mask) = 3;  % B=3

% Fielding function
WORLD_FIELDING_FUNCTION = ones(WORLD_FRAME_SHAPE);

%% Functions

function image_out = apply_color_weights(image_in, RGB_SCALARS, R_mask, G_mask, B_mask)
    image_double = double(image_in);
    image_double(R_mask) = image_double(R_mask) * RGB_SCALARS(1);
    image_double(G_mask) = image_double(G_mask) * RGB_SCALARS(2);
    image_double(B_mask) = image_double(B_mask) * RGB_SCALARS(3);
    image_out = uint8(min(max(round(image_double), 0), 255));
end

function [rgb_img, fig] = debayer_image(image, R_mask, G_mask, B_mask, WORLD_RGB_SCALARS, visualize_results)
    corrected = double(image);
    corrected(R_mask) = corrected(R_mask) * WORLD_RGB_SCALARS(1);
    corrected(G_mask) = corrected(G_mask) * WORLD_RGB_SCALARS(2);
    corrected(B_mask) = corrected(B_mask) * WORLD_RGB_SCALARS(3);
    corrected = uint8(min(max(round(corrected), 0), 255));

    % Debayering (assumes BayerRG pattern)
    debayered_no_corr = demosaic(image, 'rggb');
    debayered_with_corr = demosaic(corrected, 'rggb');

    fig = [];
    if visualize_results
        fig = figure;
        subplot(1, 3, 1); imshow(image); title('Before');
        subplot(1, 3, 2); imshow(debayered_no_corr); title('After [no color correction]');
        subplot(1, 3, 3); imshow(debayered_with_corr); title('After [with color correction]');
    end

    rgb_img = debayered_with_corr;
end

function scalars = calculate_color_weights(chunks, R_mask, G_mask, B_mask)
    R_means = [];
    G_means = [];
    B_means = [];

    for chunk = chunks
        frames = chunk.W.v;
        for i = 1:size(frames, 3)
            frame = frames(:, :, i);
            R_means(end+1) = mean(frame(R_mask), 'all');
            G_means(end+1) = mean(frame(G_mask), 'all');
            B_means(end+1) = mean(frame(B_mask), 'all');
        end
    end

    mean_R = mean(R_means);
    mean_G = mean(G_means);
    mean_B = mean(B_means);

    scalars = [1, mean_R / mean_G, mean_B / mean_R];
end

%% AGC Placeholder Functions
function settings = retrieve_settings_ranges(camera_id)
    % Placeholder - replace with actual implementation or MEX call
    settings = containers.Map();
    settings('Again') = [1.0, 10.0];
    settings('Dgain') = [1.0, 10.0];
    settings('exposure') = [100, 10000];
end

function states = retrieve_discrete_states(camera_id)
    % Placeholder - replace with actual implementation or MEX call
    states = containers.Map();
    states('Again') = struct('min', 1.0, 'max', 10.0);
    states('Dgain') = struct('min', 1.0, 'max', 10.0);
    states('exposure') = struct('min', 100, 'max', 10000);
end

%% Main
function main()
    % Placeholder for testing
    % Load an image and test debayer
    test_img = imread('bayer_image.png'); % Load a Bayer RG image
    [rgb_img, ~] = debayer_image(test_img, R_mask, G_mask, B_mask, WORLD_RGB_SCALARS, true);
end

% Uncomment the line below to run in script mode
% main();
