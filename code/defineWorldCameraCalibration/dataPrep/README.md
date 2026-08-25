Description of routines used to prepare raw files and populate the `data`
directory of the lightLoggerAnalysis repository.

## `deriveEmpircalAGCAndIlluminance.py`

This runnable Python derivation creates the empirical data consumed by the
AGC-to-illuminance fitter. It reads the shared temporal alignment from
`derived/cameraAGCLag.mat`, processes the supplied raw `GKA` recordings in
memory, aligns the camera score with filtered minispect illuminance, and
selects the calibration points.

It writes `data/empircalAGCAndIlluminance.mat`. The output contains the MATLAB
struct `empiralAGC` with these fields:

- `cameraScoreLinear`: analog gain × digital gain × exposure;
- `msIlluminance`: matched, temporally filtered illuminance in lux; and
- `sharedLagSeconds`: the lag loaded from `derived/cameraAGCLag.mat`.

The script does not fit the final camera conversion model. Its output feeds
[`fitCameraAGCToIlluminance.m`](../utilities/fitCameraAGCToIlluminance.m), which
writes the production coefficients and fitting metadata to
`derived/cameraAGCToIlluminanceFit.mat`. That derived fit is subsequently
loaded by `video_io.camera_scores_to_illuminance` during video processing.

```text
derived/cameraAGCLag.mat
    -> deriveEmpircalAGCAndIlluminance.py
    -> data/empircalAGCAndIlluminance.mat
    -> ../utilities/fitCameraAGCToIlluminance.m
    -> derived/cameraAGCToIlluminanceFit.mat
```
