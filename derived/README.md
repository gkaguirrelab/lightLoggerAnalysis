Stored here are results of various calibration operations. These parameters are subsequently used in the processing pipeline of recordings using the light logger.

## `cameraAGCLag.mat`

The output of
[`deriveAGCLag.py`](../code/defineWorldCameraCalibration/deriveAGCLag.py). The
MATLAB struct `cameraAGCLag` contains the shared camera-response lag in seconds
and included / skipped recording counts. A fresh raw-recording run also
populates the candidate lag grid and mean correlation at each candidate; those
diagnostic arrays are empty in the migrated legacy result until it is rerun.
Positive lag means that the camera response follows the minispect signal. The
empirical AGC/illuminance data-prep script reads this file before exporting its
calibration point cloud. Regenerate it by passing the raw `GKA` recording paths
directly to `deriveAGCLag.py`.

## `cameraAGCToIlluminanceFit.mat`

The output of
[`fitCameraAGCToIlluminance.m`](../code/defineWorldCameraCalibration/utilities/fitCameraAGCToIlluminance.m).
The `cameraAGCToIlluminanceFit` struct contains the four coefficients of the
continuous two-slope piecewise model in log10 space, named coefficient fields,
the calibration sample count and ranges, the shared temporal lag, and fitting
provenance. Its fields are:

- `parameterVector`: coefficients in the order slope below the breakpoint,
  intercept below the breakpoint, slope above the breakpoint, and the log10
  camera-score breakpoint;
- `slopeBelowBreakpoint`, `interceptBelowBreakpoint`,
  `slopeAboveBreakpoint`, and `log10CameraScoreBreakpoint`: named copies of
  those coefficients;
- `sampleCount`, `cameraScoreRange`, and `illuminanceRange`: the fitting-data
  extent;
- `sharedLagSeconds`: the temporal alignment inherited from the empirical
  point cloud; and
- `sourceDataFile`, `model`, and `objective`: fitting provenance.

The diagnostic script
[`fitEmpircalAGCtoIlluminance.m`](../code/defineWorldCameraCalibration/utilities/fitEmpircalAGCtoIlluminance.m)
calls the same fitter and therefore refreshes this file while producing its
plot and console report. `defineAllScript.m` also regenerates it as part of the
complete world-camera definition sequence.

Python's `video_io.camera_scores_to_illuminance` loads this file and evaluates
the saved model without starting MATLAB or repeating the fit. Its result feeds
`video_io.video_to_illuminance`, which scales processed world-camera frames
into illuminance values. Minispect count calibration remains a separate MATLAB
operation.
