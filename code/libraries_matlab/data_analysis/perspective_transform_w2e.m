function perspective_transform_w2e(frame, geometric_transform)
% VISUALIZE_PERSPECTIVE_3D Creates a 3D surface-like visualization of the 
% projective transformation.
%
% This function maps the 2D pixel grid of the image onto a 3D plot to show
% the warping effect defined by the projective transformation.

    arguments
        frame;                 % The original image frame.
        geometric_transform;   % The matlab.hgf.GeometricTransformation2D object.
    end

    % --- 1. Get the Image Size and Grid ---
    [rows, cols, ~] = size(frame);
    
    % Create a regular grid of coordinates for the original frame
    [X_in, Y_in] = meshgrid(1:cols, 1:rows);
    
    % Prepare the coordinates for transformation: [x, y] pairs
    input_points = [X_in(:), Y_in(:)];
    
    % --- 2. Apply the Forward Transformation ---
    % Transform the grid points using the geometric transformation
    output_points = transformPointsForward(geometric_transform, input_points);
    
    % Reshape the transformed points back into the grid structure
    X_out = reshape(output_points(:, 1), rows, cols);
    Y_out = reshape(output_points(:, 2), rows, cols);
    
    % --- 3. Create the 3D Plot (Mesh/Surface) ---
    figure('Name', '3D Visualization of Projective Transformation');
    
    % Create the Z-data for the plot. You can use the R, G, or B channel 
    % of the image as the "height" (Z-axis) to get a textured 3D surface.
    % Let's use the Red channel for Z-data.
    Z_data = double(frame(:, :, 1)); 

    % Use the Z_data to color the surface (or 'C' in surf/mesh).
    % Use the original X_in/Y_in as the base plane.
    % The X_out/Y_out are the 'unwarped' coordinates in the target space.
    
    % --- Option A: Plotting the Warped Surface in 3D Space ---
    % Plot the transformation: X_out and Y_out define the warped grid, 
    % and Z_data provides a "texture" or height from the original image.
    % Set Z to zero for a 'flat' projection visualization.
    Z_flat = zeros(rows, cols);

    h = surf(X_out, Y_out, Z_flat, Z_data, 'EdgeColor', 'none');
    
    % Optional: Plot as a mesh to see the grid lines clearly
    % h = mesh(X_out, Y_out, Z_flat, Z_data); 
    
    colormap('gray'); % Use a grayscale colormap based on the Red channel intensity
    axis tight;
    xlabel('Transformed X Coordinate');
    ylabel('Transformed Y Coordinate');
    zlabel('Z (Set to 0)');
    title('Projective Transform as a Warped 2D Surface (Mesh/Surf Plot)');
    
    % Improve the view
    view(3); % View from 3D perspective
    camlight; % Add lighting
    lighting phong; % Set lighting style
    
    % Overlay the original image as a texture if desired (more complex)
    % Set the CData (color data) of the surface to the entire RGB frame
    if size(frame, 3) == 3
        set(h, 'CData', frame);
        set(h, 'FaceColor', 'texturemap');
    end

end