% Coordinates of targets
degPositions = [ ...
     0,  0;  -20, 20;   -20, -20;   20, 20;   20, -20; ...
     0, 20;   0, -20;   -20,   0;   20, 0; ...
    -15, 15;  15, 15;   -15, -15;   15, -15; ...
    -10,  10;   -10, -10;    10,  10;    10, -10; ...
     0,  10;    0, -10;   -10,   0;     10,   0; ...
     -5,   5;   5,   5;   -5,  -5;     5,  -5;     0,   0];

% Repeat each position 3 times
repeats = 3;
positionsRepeated = repmat(degPositions, repeats, 1);

% Create time vector: each point lasts 1 second
numPoints = size(positionsRepeated, 1);
time = 0:(numPoints - 1);

% Plot
figure;
plot(time, positionsRepeated(:,1), 'ro-', 'LineWidth', 1.5); hold on;
plot(time, positionsRepeated(:,2), 'bo-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (degrees)');
legend('Horizontal', 'Vertical');
title('Target Position Over Time');

gazeData = load('/Users/samanthamontoya/Aguirre-Brainard Lab Dropbox/Sam Montoya/FLIC_data/lightLogger/gaze_cal_pupil_features.mat').pupil_features
