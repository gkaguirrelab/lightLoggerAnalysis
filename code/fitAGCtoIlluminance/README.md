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

The active MATLAB fitting and diagnostic files now live in
`../defineWorldCameraCalibration/utilities/`. This directory retains the formal
workflow notebook and its analysis documentation. The production fit is stored
in `derived/cameraAGCToIlluminanceFit.mat`, rather than being recalculated by
Python for each processed video.

## Directory Contents

- `../defineWorldCameraCalibration/deriveAGCLag.py`: Runnable temporal-alignment
  stage that estimates one shared AGC lag and writes
  `derived/cameraAGCLag.mat`.
- `../defineWorldCameraCalibration/dataPrep/deriveEmpircalAGCAndIlluminance.py`:
  Runnable entry point that processes raw recordings in memory and writes the
  final MATLAB calibration data.
- `../defineWorldCameraCalibration/utilities/fitEmpircalAGCtoIlluminance.m`:
  MATLAB routine that loads
  `data/empircalAGCAndIlluminance.mat`, fits a robust piecewise
  log-log relationship, plots the result, reports a conversion equation, and
  writes the derived fit artifact.
- `../defineWorldCameraCalibration/utilities/fitCameraAGCToIlluminance.m`:
  Shared MATLAB fitter containing the robust piecewise fitting procedure. It
  writes `derived/cameraAGCToIlluminanceFit.mat` for Python consumers.
- `../../data/empircalAGCAndIlluminance.mat`: Linear-scale point cloud used by the MATLAB
  fitting routine. It contains the `empiralAGC` struct with the fields
  `cameraScoreLinear`, `msIlluminance`, and `sharedLagSeconds`.
- `../../derived/cameraAGCToIlluminanceFit.mat`: Saved piecewise fit
  coefficients and fitting provenance loaded by Python video processing.
- `fit_agc_to_illuminance_processing.ipynb`: Formal Python processing notebook
  that selects recordings and runs the two Python derivation stages before the
  MATLAB model fit.

## Method Summary

The method follows the group discussion from July 29 through August 11, 2026
with Geoff Aguirre, Zachary Kelly, Sam Montoya, and Sophia Mirabal. Sophia's
initial work characterized the temporal response of the AGC algorithm using
`adaptiveControlDemo.m`, comparing a first-order low-pass model with an
empirical kernel derived from simulated dark-to-bright and bright-to-dark step
responses. The current workflow uses the empirical AGC kernel stored in
`data/agc_empirical_kernels.mat`.

For each raw `GKA` recording, `deriveAGCLag.py`:

1. Loads world-camera metadata and extracts analog gain, digital gain, and
   exposure.
2. Computes the AGC product as the camera score.
3. Loads minispect AS chunks, flattens the samples, and converts AS counts to
   calibrated illuminance using `msCounts2Illuminance.m`.
4. Applies the same empirical AGC kernel to every recording and chooses one
   shared lag that maximizes mean recording-level correlation.
5. Writes the lag and its supporting search data to
   `derived/cameraAGCLag.mat`.

The data-prep exporter then reads that lag, aligns filtered illuminance and
camera score, measures processed-video saturation, applies the point-selection
filters, and writes the final MATLAB point cloud.

The formal notebook `fit_agc_to_illuminance_processing.ipynb` records the
recording-selection scope used to develop the analysis. Run
`../defineWorldCameraCalibration/deriveAGCLag.py` first to derive the lag, then
run `../defineWorldCameraCalibration/dataPrep/deriveEmpircalAGCAndIlluminance.py`
with the raw `GKA` recording paths to produce the linear-scale MATLAB file.

The final linear-scale `.mat` export is produced by
`deriveEmpircalAGCAndIlluminance.py`. The current export keeps points
that are finite, positive, below or equal to 40% frame spatial saturation, and
not among the first 100 samples of a recording. These exported points are then
loaded by `../defineWorldCameraCalibration/utilities/fitEmpircalAGCtoIlluminance.m`,
transformed into log10 space, and fit with a continuous two-slope piecewise
linear model using an L1 objective. The shared fitter saves the coefficients
and fit provenance to `derived/cameraAGCToIlluminanceFit.mat`.

The exported MATLAB struct records the shared lag alongside the matched
camera-score and illuminance vectors so the temporal-alignment provenance is
carried into the model-fitting stage.

## Applying the Fit to Video

`video_to_illuminance` in
`code/library/matlabIO/python_libraries/video_io.py` applies the final fit to a
processed `W.avi` and its matching raw `GKA` directory. For each captured
frame, it:

1. Multiplies analog gain, digital gain, and exposure from the raw world
   metadata.
2. Converts that camera score to the illuminance represented by the AGC target
   using the fitted continuous two-slope log-log equation.
3. Divides each nonsaturated value in the processed, linearized RGB frame by
   the linearized camera target and scales the fitted frame illuminance by the
   resulting relative-intensity factor.
4. Optionally averages the resulting per-channel illuminance values to produce
   one estimate per video frame.

Each valid RGB channel is converted to illuminance before any spatial
reduction. By default, the result is one spatially averaged illuminance estimate
in lux per video frame. Passing `result_as_mean=False` instead returns the full
`(frames, rows, columns, 3)` array with every RGB channel expressed in lux.
Inserted dummy frames, completely saturated frames, and saturated spatial
pixels are represented by `NaN`. The function requires the standard processed
world video with camera response linearization enabled; it is not valid for a
raw or gamma-encoded video. The camera-score conversion loads
`derived/cameraAGCToIlluminanceFit.mat` and evaluates the saved coefficients in
NumPy; it no longer starts MATLAB or repeats the fit during video processing.
The linearized AGC target is calculated at runtime by passing the raw target of
127 through `world_util.linearize_camera_responsivity`.

The function also loads minispect chunks from the matching raw GKA recording,
selects the production spectral channels `AS_0-AS_7`, and converts their counts
with `ms_util.ms_counts_to_illuminance`. This shared helper is also used during
fitting and delegates calibration to MATLAB's `msCounts2Illuminance.m`. The
function returns `(camera_illuminance, camera_t, ms_illuminance, ms_t)`, with
each sensor retaining its native timestamps and sample rate. For a segmented
video, MS samples are restricted to the selected camera interval. Passing
`visualize_results=True` plots finite camera and MS samples against these real
timestamps through `video_io.plot_video_illuminance`.

When the processed video is a `tag` or `task` segment, `video_to_illuminance`
loads the sibling `tag_task_start_end.mat` file and selects the matching field
from its `tag_task_start_end` struct. These frame ranges are stored as
one-based, inclusive MATLAB indices and are converted to a zero-based Python
slice before the raw metadata are aligned with the segmented video.

## Current Included Recordings

The current `data/empircalAGCAndIlluminance.mat` file includes
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
to illuminance.

## Notes From Discussion

The analysis originally used the reciprocal of the AGC product as the camera
score, but the group switched to the natural product of the AGC values on
August 7, 2026. This makes the relationship with illuminance have a negative
slope, but preserves the direct interpretation of the AGC settings.

The group also investigated whether different apparent slopes could arise from
camera saturation, AGC transients, or different exposure ranges in recordings
collected at 120 FPS versus 180 FPS. The current workflow therefore includes
frame saturation diagnostics, derivative-colored AGC product plots, and an
option to exclude the initial samples from each recording.
