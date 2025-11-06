% Local function to find the min square figsize required to plot data 
function [rows, cols] = find_min_figsize(num_plots)
    % Iterate over the ints between 1 and num_plots 
    for ii = 1:num_plots
        rows = ii; 
        cols = ii; 

        % Determine if we have reached the target 
        if(rows * cols >= num_plots)
            return ; 
        end 

    end 

end 