% 1, drag the data into Neon Player. 
% 2, plugin manager > enable blinks and rawdataexporter.
% 3, ctrl + e to export the data to csv.
% load the data by dragging and dropping the blinks.csv nd gaze_positions
% csv
gazePos = table2array(gaze_positions);
blinks = table2array(blinks);
eyeStates = table2array(x3d_eye_states);

nTime = length(gazePos);

eyeInfo = nan(nTime, 7);

eyeInfo(:,1) = gazePos(:,4);%1 azi (deg)
eyeInfo(:,2) = gazePos(:,5);%2 ele (deg)
eyeInfo(:,3) = eyeStates(:,2);%3 pupil size left (mm)
eyeInfo(:,4) = eyeStates(:,3);%4 pupil size right (mm)
eyeInfo(:,5) = eyeStates(:,20);%5 eye openness arc length left (mm)
eyeInfo(:,6) = eyeStates(:,21);%6 eye openness arc length right (mm)
eyeInfo(:,7) = gazePos(:,1)./1e9;%7 timestamp (s)

% Convert blink nanoseconds to seconds)
b_start_sec = blinks(:, 2) ./ 1e9;
b_end_sec = blinks(:, 3) ./ 1e9;

% Loop through each detected blink
for i = 1:length(b_start_sec)
    % Find logical indices where the timestamp is within the blink window
    blinkMask = (eyeInfo(:,7) >= b_start_sec(i)) & (eyeInfo(:,7) <= b_end_sec(i));
    
    % Set columns 1 through 6 to NaN for these timepoints
    eyeInfo(blinkMask, 1:6) = NaN;
end

%% Save!!

%% Sanity check plots
% plot gaze positions
figure; hold on; plot(eyeInfo(:,7),eyeInfo(:,1)); plot(eyeInfo(:,7),eyeInfo(:,2))

%Get Y-axis limits to make the bars span the whole height
y_lims = ylim; 

% Loop through each blink and draw a transparent rectangle
for i = 1:length(start_times)
    % Define the 4 corners of the rectangle
    % X: [start, end, end, start]
    % Y: [bottom, bottom, top, top]
    x_patch = [b_start_sec(i), b_end_sec(i), b_end_sec(i), b_start_sec(i)];
    y_patch = [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];
    
    % Draw the patch
    patch(x_patch, y_patch, 'r', ...         % 'r' for red bars
          'FaceAlpha', 0.2, ...              % Transparency (0 to 1)
          'EdgeColor', 'none');              % Removes the black outline
end

hold off;
xlabel('Time (s)');
ylabel('Gaze Angle (deg)');
title('Gaze Pos with Blink Overlays');

figure; hold on; plot(eyeInfo(:,7),eyeInfo(:,3)); plot(eyeInfo(:,7),eyeInfo(:,4))
