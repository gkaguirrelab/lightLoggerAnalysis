# Fit AGC to Illuminance

This directory contains the analysis used to estimate environmental
illuminance from the world camera automatic-gain-control settings. The goal
is to convert the camera AGC product,

```text
camera_score = analog_gain * digital_gain * exposure
```

into an implied steady-state illuminance value. The minispectrometer provides
the external illuminance reference, while the camera AGC settings provide the
quantity that will be available during ordinary video processing.

## Directory Contents

- `fit_agc_to_illuminance_util.py`: Python workflow for loading raw recordings,
  converting minispect AS counts to illuminance, applying the temporal AGC
  response model, measuring world-video saturation, selecting usable points,
  and writing cached processing data.
- `fitEmpircalAGCtoIlluminance.m`: MATLAB routine that loads
  `camera_agc_illuminance_linear_scale.mat`, fits a robust piecewise
  log-log relationship, plots the result, and reports a conversion equation.
- `fitCameraAGCToIlluminance.m`: Shared MATLAB function containing the exact
  robust piecewise fitting procedure and conversion. It is called by both the
  MATLAB plotting script and Python's `video_to_illuminance` through the MATLAB
  Engine API.
- `camera_agc_illuminance_linear_scale.mat`: Linear-scale x/y point cloud used
  by the MATLAB fitting routine. This file is intentionally stored at the top
  level of this analysis directory, not in the cache subfolder.
- `fit_agc_to_illuminance_processing.ipynb`: Formal Python processing notebook
  that selects recordings, runs the empirical AGC-kernel processing fit, writes
  cached data, renders the calibration-selection dashboard, and generates
  `camera_agc_illuminance_linear_scale.mat`.
- `cached_processing_data/`: Generated data that allow downstream plotting and
  diagnostics without reprocessing all raw recordings.

## Method Summary

The method follows the group discussion from July 29 through August 11, 2026
with Geoff Aguirre, Zachary Kelly, Sam Montoya, and Sophia Mirabal. Sophia's
initial work characterized the temporal response of the AGC algorithm using
`adaptiveControlDemo.m`, comparing a first-order low-pass model with an
empirical kernel derived from simulated dark-to-bright and bright-to-dark step
responses. The current workflow uses the empirical AGC kernel stored in
`misc/agc_simulation/agc_empirical_kernels.mat`.

For each raw `GKA` recording, `fit_agc_to_illuminance_util.py`:

1. Loads world-camera metadata and extracts analog gain, digital gain, and
   exposure.
2. Computes the AGC product as the camera score.
3. Locates the matching processed `W.avi` video and measures each frame's
   spatial saturation.
4. Loads minispect AS chunks, flattens the samples, and converts AS counts to
   calibrated illuminance using `msCounts2Illuminance.m`.
5. Applies the same empirical AGC kernel to every recording and chooses one
   shared lag that maximizes mean recording-level correlation, unless a fixed
   lag is supplied.
6. Aligns filtered illuminance and camera score on a common timebase.
7. Saves plot-ready cached point clouds so saturation thresholds and display
   choices can be revisited without rerunning the full video/sensor workflow.

The formal notebook `fit_agc_to_illuminance_processing.ipynb` is the intended
entry point for running this Python stage of the analysis. It records the
recording-selection scope, invokes `fit_agc_to_illuminance`, and then reloads
the cached point cloud to produce the final linear-scale `.mat` file for the
MATLAB model fit.

The final linear-scale `.mat` export is produced by
`plot_frame_saturation_from_processed_data`. The current export keeps points
that are finite, positive, below or equal to 40% frame spatial saturation, and
not among the first 100 cached samples of a recording. These exported points
are then loaded by `fitEmpircalAGCtoIlluminance.m`, transformed into log10
space, and fit with a continuous two-slope piecewise linear model using an L1
objective.

The Python diagnostic dashboard also shows a correlation-qualified contextual
fit restricted to recordings with shared-model correlation at least 0.9. This
contextual fit helps inspect the most temporally coherent recordings, but it
does not change which points are written to `camera_agc_illuminance_linear_scale.mat`.

## Applying the Fit to Video

`video_to_illuminance` in
`code/library/matlabIO/python_libraries/video_io.py` applies the final fit to a
processed `W.avi` and its matching raw `GKA` directory. For each captured
frame, it:

1. Multiplies analog gain, digital gain, and exposure from the raw world
   metadata.
2. Converts that camera score to the illuminance represented by the AGC target
   using the fitted continuous two-slope log-log equation.
3. Computes the mean nonsaturated sensor value in the processed, linearized
   RGB frame.
4. Divides that mean by the linearized camera target and scales the fitted
   frame illuminance by the resulting relative-intensity factor.

Each valid RGB channel is converted to illuminance before any spatial
reduction. By default, the result is one spatially averaged illuminance estimate
in lux per video frame. Passing `result_as_mean=False` instead returns the full
`(frames, rows, columns, 3)` array with every RGB channel expressed in lux.
Inserted dummy frames, completely saturated frames, and saturated spatial
pixels are represented by `NaN`. The function requires the standard processed
world video with camera response linearization enabled; it is not valid for a
raw or gamma-encoded video. The camera-score conversion is performed by
`fitCameraAGCToIlluminance.m`, so video processing and the MATLAB diagnostic
script always use the same fitting procedure and calibration point cloud. The
linearized AGC target is calculated at runtime by passing the raw target of 127
through `world_util.linearize_camera_responsivity`.

## Current Included Recordings

The current `camera_agc_illuminance_linear_scale.mat` file includes
point-filtered data from 111 processed recordings, spanning 17 subjects and
10 activities:

```text
Subjects: FLIC_20, FLIC_21, FLIC_28, FLIC_39, FLIC_42, FLIC_51,
          FLIC_1029, FLIC_1034, FLIC_1038, FLIC_1044, FLIC_1047,
          FLIC_2001, FLIC_2002, FLIC_2003, FLIC_2004, FLIC_2005,
          FLIC_2006
Activities: chat, gazeCalibration, lunch, phone, read, sitBiopond,
            walkBiopond, walkIndoor, walkOutdoor, work
```

The correlation-qualified diagnostic subset contains 20 recordings from 15
subjects and 6 activities:

```text
FLIC_20: read, walkBiopond
FLIC_21: read
FLIC_28: read
FLIC_42: read
FLIC_51: read
FLIC_1029: read
FLIC_1034: read, sitBiopond
FLIC_1038: read, sitBiopond
FLIC_1044: read
FLIC_1047: read
FLIC_2002: gazeCalibration
FLIC_2003: work
FLIC_2004: work
FLIC_2005: gazeCalibration, walkOutdoor, work
FLIC_2006: work
```

Dark recordings are processed far enough to expose conversion failures, but
they are excluded from pooled fits when minispect counts cannot be converted
to illuminance. Lower-correlation recordings remain represented in the cached
processing archive so thresholds can be re-evaluated later.

## Notes From Discussion

The analysis originally used the reciprocal of the AGC product as the camera
score, but the group switched to the natural product of the AGC values on
August 7, 2026. This makes the relationship with illuminance have a negative
slope, but preserves the direct interpretation of the AGC settings.

The group also investigated whether different apparent slopes could arise from
camera saturation, AGC transients, or different exposure ranges in recordings
collected at 120 FPS versus 180 FPS. The current workflow therefore includes
frame saturation diagnostics, derivative-colored AGC product plots, and an
option to exclude the initial cached samples from each recording.
