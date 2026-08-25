# World Camera Calibration Utilities

This directory contains the MATLAB fitting and diagnostic routines used at the
end of the empirical world-camera AGC calibration workflow.

## `fitCameraAGCToIlluminance.m`

This is the production fitting function. It:

1. Loads `data/empircalAGCAndIlluminance.mat` and reads the `empiralAGC`
   structure created by
   `dataPrep/deriveEmpircalAGCAndIlluminance.py`.
2. Keeps finite, positive camera-score and minispect-illuminance pairs.
3. Transforms both quantities into log10 space.
4. Fits a continuous two-slope piecewise-linear model by minimizing the sum of
   absolute log10 illuminance residuals.
5. Writes the fitted calibration to
   `derived/cameraAGCToIlluminanceFit.mat` every time the function runs.

The saved MAT file contains a structure named
`cameraAGCToIlluminanceFit`. It records the original four-parameter vector,
named coefficient fields, the calibration ranges and sample count, the shared
AGC lag, the source-data path, and descriptions of the model and objective.

The function may also be given camera scores to evaluate. It returns predicted
illuminance in lux together with the fitted parameter vector and the cleaned
calibration vectors. Nonfinite and nonpositive input scores produce `NaN`.

`defineAllScript.m` calls this function after the other world-camera
calibrations, so running the complete definition sequence regenerates the
derived fit artifact.

## `fitEmpircalAGCtoIlluminance.m`

This is the human-facing diagnostic script. It calls
`fitCameraAGCToIlluminance`, so it also regenerates
`derived/cameraAGCToIlluminanceFit.mat`. It then:

- prints the fitted piecewise equation;
- plots the empirical calibration points and fitted curve;
- marks the fitted breakpoint; and
- overlays the independent integrating-sphere measurements.

It does not create a second calibration file. Its plot and console output are
diagnostics for inspecting the production fit saved by the shared function.
The historical `Empircal` spelling is retained in the filename to avoid
breaking existing calls.

## Inputs and Downstream Use

```text
data/empircalAGCAndIlluminance.mat
    -> fitCameraAGCToIlluminance.m
    -> derived/cameraAGCToIlluminanceFit.mat
    -> video_io.camera_scores_to_illuminance
    -> video_io.video_to_illuminance
```

Python loads the saved coefficients and evaluates the same piecewise equation
with NumPy. It does not start MATLAB or repeat the camera fit during video
processing. MATLAB is still used separately by the video pipeline when
minispect counts must be converted to calibrated illuminance. The derived MAT
file is the authoritative source for these coefficients; `world_util.py` no
longer carries a duplicate hard-coded parameter vector.
