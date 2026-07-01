% Local function to find the min square figsize required to plot data 
function [rows, cols] = find_min_figsize(num_plots)
% Compute a square-ish subplot grid large enough for a set of plots.
%
% Syntax:
%   [rows, cols] = find_min_figsize(num_plots)
%
% Description:
%   This helper increments a square grid size until the number of slots in
%   the grid is at least as large as the requested number of plots. It is
%   mainly used when downstream code wants a compact, approximately square
%   tiling without separately reasoning about rows and columns.
%
% Inputs:
%   num_plots                - Positive integer. Number of axes or panels
%                              that must fit in the layout.
%
% Outputs:
%   rows                     - Integer. Number of rows in the chosen grid.
%   cols                     - Integer. Number of columns in the chosen
%                              grid.
%
% Examples:
%{
    [rows, cols] = find_min_figsize(7);
%}

    for ii = 1:num_plots
        rows = ii; 
        cols = ii; 

        % Determine if we have reached the target 
        if(rows * cols >= num_plots)
            return ; 
        end 

    end 

end 
