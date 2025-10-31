function plotPupilCenters(frameSet, perimeter, whichFrames)
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