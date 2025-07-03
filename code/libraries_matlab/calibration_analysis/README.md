

This README file provides the documentation for performing spatial calibration on videos recorded with the *GKA Lab High-Speed, Personal Light Logger*. The goal of this analysis is to align the spatial mapping of the device’s world camera with known calibration targets, enabling accurate extraction of spatial and temporal properties for experimental use.

Chiefly described is the methods and strategy implemented for spatial calibration of this device, primarily camera intrinsics estimation and camera FOV angle calculation. 

## Background
We implemented MATLAB's Single Camera Calibrator App to estimate camera intrinsic parameters, including focal length, principal point, and lens distortion coefficients, using a checkerboard calibration pattern. This was done to determine reliable spatial mappings for future applications involving spatiotemporal video analysis in experimental settings. This calibration corrects for lens distortion and enables precise pixel-to-degree conversions.

## Setup - Spatial Calibration
1. The custom detector pattern we chose for such calibration was a checkerboard, printed on Letter-sized (8.5 x 11 in) paper and attatched to the front of a hard cover textbook. The paper included a white border (~1 in) around the pattern, with labeled x- and y-axes for clarity.
2. The Light Logger was placed on a tripod at a height of ~100 cm, level with the ground. Using a tape measure and sticky notes, we measured a 1 m (100 cm) arc around the camera, placing sticky notes at regular intervals to guide consistent distance during calibration. This arc extended from one side of the camera to the other.
3. After powering on the camera and allowing AGC to stabilize, we began navigating the visual field parameters along the desired axis. Beginning with the far left side of the arc, one person held the detector pattern at range of heights, holding the pattern angled upward when positioned low, level when at camera height, and angled downward when positioned high. At about every six inches of vertical movement, the pattern was tilted ~30° left and right. This process continued across the arc. All movements were made slowly to ensure sharp feature visibility and optimal edge detection for calibration.

## Setup - FOV Measurement
1. The Light Logger’s height and position on the tripod remained constant during FOV calculation, this time centered and directed at a blank wall for measurement.
2. Using a tape measure, we estimated the camera’s center point on the wall and marked it with a sticky note. On either side of this center point, we placed sticky notes at 10 cm intervals for several meters.
3. After powering on the camera and allowing AGC to stabilize, we slowly rotated the device along its x and y axes to minimize error, ensuring we captured a frame where the camera faced the wall head-on.
4. In the captured image, we determined the horizontal field by counting the number of sticky notes visible on each side of the central marker (17 on each side). Calculating 17 x 10 = 170 cm per side plus 7.6 cm to account for an additional note, we measured a total of 355.2 cm at a 100 cm distance. We used this value to trigonometrically derive the horizontal FOV angle as ~120°
5. Given the horizontal FOV angle and the camera’s known image aspect ratio (480 x 640), we geometrically and trigonometrically derived the diagonal FOV angle as ~150°. We also captured an indirect calculation with the device’s fisheye projection model using the intrinsics determined by the calibration. This came out to ~167°.