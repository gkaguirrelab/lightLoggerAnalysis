# Python MATLAB-I/O Libraries

This directory contains Python utilities that read light-logger data and
interoperate with MATLAB calibration routines and MAT files.

## Camera-score to mean-luminance conversion

`video_io.camera_scores_to_mean_luminance` converts the world-camera AGC
product

```text
analog gain × digital gain × exposure
```

to average scene luminance. It loads the `cameraScore` and
`avgSceneLuminance` lookup vectors from
`derived/cameraScoreToAverageLuminance.mat` and interpolates in log10 space, as
specified by `defineAGCToMeanLuminance.m`. Nonfinite and nonpositive camera
scores remain `NaN`, as do scores outside the calibrated range.

`video_io.camera_scores_to_illuminance` applies the corresponding Lambertian
conversion (illuminance = π × mean luminance) for the existing video-to-lux
workflow.

The video workflow also converts minispect counts through
`ms_util.ms_counts_to_illuminance`. That is a separate calibration path and
delegates to MATLAB.
