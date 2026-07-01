function plotPupilCenters(frameSet, perimeter, whichFrames)
% Interactively display pupil perimeter points and centers for selected frames
%
% Syntax:
%   plotPupilCenters(frameSet, perimeter, whichFrames)
%
% Description:
%   For each index in whichFrames, retrieves the pupil perimeter data at
%   the corresponding frame from frameSet, plots the perimeter points and
%   their Euclidean center, and pauses for a keypress before advancing to
%   the next frame. Useful for visually inspecting pupil boundary quality
%   during gaze calibration.
%
% Inputs:
%   frameSet              - Numeric vector. Frame indices into the
%                           perimeter data cell array.
%   perimeter             - Struct. Perimeter data struct with a .data
%                           field containing a cell array of per-frame
%                           structs with fields Xp and Yp.
%   whichFrames           - Numeric vector. Indices into frameSet
%                           specifying which frames to display.
%
% Outputs:
%   none (creates an interactive figure)
%
% Examples:
%{
    load('perimeter.mat', 'perimeter');
    frameSet = [10349; 10666; 11507];
    plotPupilCenters(frameSet, perimeter, 1:3);
%}
figure;
all_centers = []; % Initialize array to store centers for plotting the trajectory later

for ii = whichFrames 
    % 1. Get Perimeter Points
    Xp = perimeter.data{frameSet(ii)}.Xp; 
    Yp = perimeter.data{frameSet(ii)}.Yp;
    
    % Check for empty data
    if ~isempty(Xp)
        
        % 2. Calculate Euclidean Center (Simple Average)
        center_x = mean(Xp);
        center_y = mean(Yp);
        
        % Store the center for later use
        all_centers = [all_centers; center_x, center_y];
        
        % 3. Plot the Perimeter and the Center
        plot(Xp, Yp, '.'); % Plot perimeter points
        hold on; 
        plot(center_x, center_y, 'o', 'MarkerSize', 8); % Plot center
        hold on; 
        title(['Frame: ', num2str(frameSet(ii))]); % Add frame number to title
        xlabel('X Position');
        ylabel('Y Position');
        
        % PAUSE HERE: Waits for any key press to advance to the next frame
        pause; 
        
    else
        fprintf('Frame %d has no perimeter points. Skipping.\n', frameSet(ii));
    end
end
end