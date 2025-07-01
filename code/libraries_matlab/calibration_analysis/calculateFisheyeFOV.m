function fov_degrees = calculateFisheyeFOV(fisheyeIntrinsics)
% calculateFisheyeFOV_Revised Calculates the field of view (FOV) of a fisheye camera.
%
%   fov_degrees = calculateFisheyeFOV_Revised(fisheyeIntrinsics) calculates the
%   diagonal field of view in degrees for a fisheye camera using the
%   Scaramuzza model. This version interprets the MappingCoefficients as
%   direct coefficients for the polynomial mapping from angle (theta) to
%   distorted pixel radial distance (r_pixels).
%
%   Inputs:
%       fisheyeIntrinsics - A camera.FisheyeIntrinsics object obtained
%                           from fisheye camera calibration in MATLAB.
%
%   Outputs:
%       fov_degrees       - The diagonal field of view of the camera in degrees.
%
% Usage:
%{
    data = load('~/FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/camera_intrinsics_calibration.mat');
    fisheyeIntrinsics = data.camera_intrinsics_calibration.camera_intrinsics.Intrinsics;
    fov_degrees = calculateFisheyeFOV(fisheyeIntrinsics);
%}

% Extract relevant properties from the fisheyeIntrinsics object
mappingCoefficients = fisheyeIntrinsics.MappingCoefficients; % [a0 a2 a3 a4]
distortionCenter = fisheyeIntrinsics.DistortionCenter;       % [cx cy]
imageSize = fisheyeIntrinsics.ImageSize;                     % [mrows ncols]

cx = distortionCenter(1);
cy = distortionCenter(2);
img_width = imageSize(2);
img_height = imageSize(1);

% 1. Calculate the maximum radial distance (r_max_pixels) in pixels.
%    This is the distance from the distortion center to the furthest image corner.
corners = [0, 0;
           img_width, 0;
           0, img_height;
           img_width, img_height];

radial_distances = zeros(4, 1);
for i = 1:4
    radial_distances(i) = sqrt((corners(i,1) - cx)^2 + (corners(i,2) - cy)^2);
end
r_max_pixels = max(radial_distances);

% Debugging: Print values
fprintf('Debug: r_max_pixels = %.4f\n', r_max_pixels);
fprintf('Debug: Mapping Coefficients [a0 a2 a3 a4] = [%e %e %e %e]\n', ...
        mappingCoefficients(1), mappingCoefficients(2), ...
        mappingCoefficients(3), mappingCoefficients(4));

% 2. Define the polynomial function and solve for theta_max using fzero.
%    Based on the Scaramuzza model with MATLAB's [a0 a2 a3 a4] coefficients,
%    the common interpretation is that r_pixels(theta) is:
%    r_pixels = a0*theta + a2*theta^3 + a3*theta^4 + a4*theta^5
%    We are solving for theta when r_pixels(theta) = r_max_pixels.
%    So, we need to find the root of:
%    f(theta) = a0*theta + a2*theta^3 + a3*theta^4 + a4*theta^5 - r_max_pixels = 0

a0 = mappingCoefficients(1); % This acts as the primary focal length / scale factor
a2 = mappingCoefficients(2); % Coefficient for theta^3 (not theta^2)
a3 = mappingCoefficients(3); % Coefficient for theta^4
a4 = mappingCoefficients(4); % Coefficient for theta^5

% Define the anonymous function for which we want to find the root
% The input 'theta_val' is the angle in radians
func_to_solve = @(theta_val) (a0 * theta_val + ...
                              a2 * theta_val.^3 + ...
                              a3 * theta_val.^4 + ...
                              a4 * theta_val.^5 - r_max_pixels);

% --- Add Plotting for Debugging ---
% The maximum possible angle from the optical axis for a fisheye lens
% is pi (180 degrees). Let's test a range slightly larger than that.
theta_test_vals = linspace(0, pi + 0.1, 200); % Test up to approx 185 degrees
y_vals = arrayfun(func_to_solve, theta_test_vals); % Evaluate the function

figure;
plot(rad2deg(theta_test_vals), y_vals, 'b-');
hold on;
plot(rad2deg(theta_test_vals), zeros(size(theta_test_vals)), 'r--'); % Zero line
xlabel('Angle (degrees)');
ylabel('r_{pixels}(theta) - r_{max\_pixels}');
title('Function for fzero Root Finding: r(theta) - r_{max}');
grid on;
legend('f(theta)', 'Zero Line', 'Location', 'best');

% Check if the function crosses zero in the expected range.
val_at_0 = func_to_solve(0);
val_at_pi = func_to_solve(pi);
fprintf('Debug: Function value at 0 rad: %.4f\n', val_at_0);
fprintf('Debug: Function value at pi rad: %.4f\n', val_at_pi);

% The root must be found between 0 and pi.
% Initial guess should be in the ballpark. A common FOV for fisheye is 150-180 degrees,
% meaning theta_max is 75-90 degrees (1.3 to 1.57 radians).
% We can use the sign check to refine the search interval or initial guess.

% If `func_to_solve(0)` is negative, it means `r_max_pixels` is greater than 0.
% If `func_to_solve(pi)` is positive, it means `r_pixels(pi)` is greater than `r_max_pixels`.
% This setup means we are looking for a root where the polynomial equals `r_max_pixels`.
% Typically, r(theta) increases with theta. So func_to_solve should be negative at 0
% and positive at pi, for a root to exist within [0, pi].

% Adjust initial guess or interval based on plot/values
theta_initial_guess = deg2rad(75); % A common central angle for fisheye

% The fzero interval [0, pi] should still be appropriate if the function behaves as expected.
% If func_to_solve(0) and func_to_solve(pi) have the same sign, fzero will fail.
% This could mean the root is outside [0, pi], or no real root exists.
% Given the coefficients, a0 (294) is dominant. The higher-order terms are tiny.
% At theta=0, func_to_solve(0) should be -r_max_pixels (a large negative number).
% So we need func_to_solve(pi) to be positive.

% Let's adjust the fzero call to provide a range based on sign change, if possible
try
    if func_to_solve(0) < 0 && func_to_solve(pi) > 0
        theta_max_rad = fzero(func_to_solve, [0, pi]);
    else
        % If no sign change in [0, pi], try a wider range or a different initial guess
        % or acknowledge that the calibration might be problematic.
        warning('calculateFisheyeFOV_Revised:NoSignChange', ...
                'Function does not change sign in [0, pi]. Root may not exist in this range or calibration is problematic.');
        % Try searching with an initial guess, allowing fzero to go beyond the interval
        % if it can find a root, but this is less robust for guaranteeing a physical angle.
        theta_max_rad = fzero(func_to_solve, theta_initial_guess);
    end
catch ME
    warning('calculateFisheyeFOV_Revised:FzeroFailed', ...
            'fzero failed to find a root. Check your calibration data and mapping coefficients. Error: %s', ME.message);
    theta_max_rad = NaN; % Return NaN if root not found
    fov_degrees = NaN;
    return;
end

% 3. Calculate the full FOV in degrees
%    The calculated theta_max_rad is the half-FOV (angle from optical axis to edge).
%    The full diagonal FOV is 2 * theta_max_rad.
fov_degrees = rad2deg(2 * theta_max_rad);

end


function f_norm = getFocalLengthFromFisheyeIntrinsics(fisheyeIntrinsics)
% getFocalLengthFromFisheyeIntrinsics Extracts an effective focal length
% from a fisheyeIntrinsics object by creating a virtual pinhole camera.
%
%   f_norm = getFocalLengthFromFisheyeIntrinsics(fisheyeIntrinsics)
%   uses the undistortFisheyeImage function to generate a virtual
%   cameraIntrinsics object, from which an effective focal length
%   (fx) is extracted. This focal length is suitable for normalizing
%   radial pixel distances in the Scaramuzza model.
%
%   Inputs:
%       fisheyeIntrinsics - A camera.FisheyeIntrinsics object.
%
%   Outputs:
%       f_norm            - The effective focal length in pixels (fx).

% Create a dummy image just to pass to undistortFisheyeImage
% The content of the image doesn't matter, only its size, which comes
% from the intrinsics.
dummyImage = zeros(fisheyeIntrinsics.ImageSize(1), fisheyeIntrinsics.ImageSize(2), 'uint8');

% Undistort the dummy image to get the virtual pinhole intrinsics
% The 'OutputView' can influence the resulting focal length slightly,
% 'same' is a good default for getting a comparable view.
[~, virtualPinholeIntrinsics] = undistortFisheyeImage(dummyImage, fisheyeIntrinsics, 'OutputView', 'same');

% Extract the focal length (fx) from the virtual pinhole camera
f_norm = virtualPinholeIntrinsics.FocalLength(1); % Use fx for normalization
end