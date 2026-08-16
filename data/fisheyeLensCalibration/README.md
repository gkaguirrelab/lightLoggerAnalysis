ArduCam B0392 IMX219 Wide Angle M12 fisheye lens calibration

This README file documents the procedure for intrinsics calibration of the *GKA Lab High-Speed, Personal Light Logger*. The goal of this process is to determine the optical properties of the device’s world camera, including focal length, principal point, and lens distortion coefficients, to enable accurate lens distortion correction and the reliable mapping of image pixels to real-world angles for spatiotemporal analysis.

The Light Logger’s wide-FOV camera exhibits significant lens distortion, requiring careful intrinsics measurement to correct image warping and enable precise data extraction in experimental workflows. We used MATLAB’s Single Camera Calibrator App to estimate the camera’s intrinsics. The custom detector pattern we chose for such calibration was a checkerboard, printed on Letter-sized (8.5 x 11 in) paper with a ~1 in white border and labeled x and y axes for clarity, then attached to a rigid surface (hard-cover textbook) for data collection.

The Light Logger was mounted on a tripod at a height of approximately 100 cm above the ground, leveled with the floor. Using a tape measure and sticky notes, we marked a 100 cm radius arc around the camera’s position to ensure consistent distance during calibration, extending from one side of the camera’s field of view to the other. After powering on the camera and allowing AGC to stabilize, we slowly moved the checkerboard pattern across the arc, beginning at the far left and progressing to the far right. At each horizontal position, the pattern was moved vertically in increments, held angled upward when positioned low, level at camera height, and angled downward when positioned high. At each position, the pattern was also tilted ~30° to the left and right for robust orientations. All movements were performed slowly to ensure sharp feature capture for the calibration software.

Using the collected calibration data, we computed and exported the camera’s intrinsics, enabling lens distortion correction and accurate pixel-to-degree conversions for subsequent analyses.

The images and calibration calculation can be displayed in matlab using the command `cameraCalibrator('intrinsics_calibration_session.mat')`.

The resulting camera intrinsics file is saved in the derived folder as `arducamB0392_cameraInstrinsics.mat`.