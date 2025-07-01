function mapFisheyeFOV(v, fisheyeIntrinsics)
% projects the temporal mean of a fisheye video chunk onto a 3D visual field
% surface at 1 meter.
%
% mapFisheyeVideoChunkToVisualField(v, fisheyeIntrinsics)
%   
% takes a video chunk (frames x rows x cols) and the calibration
% fisheyeIntrinsics, and visualizes the temporally averaged frame
% remapped onto a hemisphere at 1 m.
%
% Inputs:
%   v                   - Video chunk (frames x rows x cols), double.
%   fisheyeIntrinsics   - camera.FisheyeIntrinsics object.
%
% Outputs:
%   Displays a 3D surface plot with angularly correct warping.
% 
% Example Usage:
%{
    data = load('~/FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/camera_intrinsics_calibration.mat');
    fisheyeIntrinsics = data.camera_intrinsics_calibration.camera_intrinsics.Intrinsics;
    mapFisheyeFOV(v, fisheyeIntrinsics)
%}

% Validate input
[frames, ~] = size(v);
if frames < 1
    error('Input video chunk has no frames.');
end

% Compute the temporal mean frame
meanFrame = squeeze(mean(v, 1)); % [rows x cols]

% Scale to [0, 1] for display
meanFrame = meanFrame - min(meanFrame(:));
meanFrame = meanFrame / max(meanFrame(:));

% Create a grayscale RGB image for visualization
image = repmat(meanFrame, [1, 1, 3]);

% Generate pixel grid
[imgHeight, imgWidth, ~] = size(image);
[xGrid, yGrid] = meshgrid(1:imgWidth, 1:imgHeight);

% Center relative to distortion center
cx = fisheyeIntrinsics.DistortionCenter(1);
cy = fisheyeIntrinsics.DistortionCenter(2);
xCentered = xGrid - cx;
yCentered = yGrid - cy;
rPixels = sqrt(xCentered.^2 + yCentered.^2);

% Invert the Scaramuzza polynomial to get theta
mappingCoefficients = fisheyeIntrinsics.MappingCoefficients;
a0 = mappingCoefficients(1);
a2 = mappingCoefficients(2);
a3 = mappingCoefficients(3);
a4 = mappingCoefficients(4);

theta = zeros(size(rPixels));
for i = 1:numel(rPixels)
    r_val = rPixels(i);
    func = @(theta_val) a0*theta_val + a2*theta_val^3 + ...
                         a3*theta_val^4 + a4*theta_val^5 - r_val;
    try
        theta(i) = fzero(func, [0, pi]);
    catch
        theta(i) = NaN; % In case the inversion fails, mask out
    end
end

% Calculate azimuth (phi)
phi = atan2(yCentered, xCentered);

% Project to 3D sphere at 1 meter
R = 1; % meters
X = R * sin(theta) .* cos(phi);
Y = R * sin(theta) .* sin(phi);
Z = R * cos(theta);

% Reshape for plotting
X = reshape(X, imgHeight, imgWidth);
Y = reshape(Y, imgHeight, imgWidth);
Z = reshape(Z, imgHeight, imgWidth);

% Mask NaN values
mask = ~isnan(X) & ~isnan(Y) & ~isnan(Z);
X(~mask) = NaN;
Y(~mask) = NaN;
Z(~mask) = NaN;

% Display the surface
figure;
surf(X, Y, Z, image, 'EdgeColor', 'none');
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Fisheye Video Chunk Temporal Mean Mapped onto 1m Visual Field');
view(3);
camlight;
lighting gouraud;
axis off;
end