Description of routines used to prepare raw files and populate the `data`
directory of the lightLoggerAnalysis repository.

## `deriveEmpircalAGCAndIlluminance.py`

This runnable Python derivation creates the empirical data consumed by the
AGC-to-illuminance fitter. It reads the shared temporal alignment from
`derived/MSIlluminanceToAGCLag.mat`, processes the supplied raw `GKA` recordings in
memory, aligns the camera score with filtered minispect illuminance, and
selects the calibration points.

It writes `data/empircalAGCAndIlluminance.mat`. The output contains the MATLAB
struct `empiralAGC` with these fields:

- `cameraScoreLinear`: analog gain × digital gain × exposure;
- `msIlluminance`: matched, temporally filtered illuminance in lux; and
- `sharedLagSeconds`: the lag loaded from `derived/MSIlluminanceToAGCLag.mat`.

The script does not define the final camera conversion. Its output is plotted
by [`defineAGCToMeanLuminance.m`](../defineAGCToMeanLuminance.m) alongside the
integrating-sphere calibration points. That MATLAB definition saves the
authoritative camera-score/mean-luminance lookup vectors in
`derived/cameraScoreToAverageLuminance.mat`.

```text
data/agc_empirical_kernels.mat
    -> ../defineMSIlluminanceToAGCKernel.py
    -> derived/MSIlluminanceToAGCKernel.mat
    -> ../defineMSIlluminanceToAGCLag.py
    -> derived/MSIlluminanceToAGCLag.mat
    -> deriveEmpircalAGCAndIlluminance.py
    -> data/empircalAGCAndIlluminance.mat
    -> ../defineAGCToMeanLuminance.m
    -> derived/cameraScoreToAverageLuminance.mat
```
