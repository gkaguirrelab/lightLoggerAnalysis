Stored here are results of various calibration operations. These parameters are subsequently used in the processing pipeline of recordings using the light logger. Every MAT artifact includes a top-level `README` variable describing its contents and provenance. MATLAB generators must write artifacts through `saveDerivedFile`, Python generators must use `save_derived_mat`, and complete derivation workflows should call `validateDerivedFilesHaveREADME` before finishing.

## `deltaSteradians.mat`

The output of
[`defineDeltaSteradians.m`](../code/defineWorldCameraCalibration/defineDeltaSteradians.m).
The 480-by-640 double array `deltaSteradians` is aligned with raw world-camera
`[row, column]` coordinates and gives the calibrated solid angle represented by
each pixel, in steradians. The definition uses the fisheye intrinsics to map
pixel centers to unit viewing directions, then calculates the local
spherical-area Jacobian with the same second-order finite-difference method as
`world_util.world_frame_visual_angle_to_steradians`. The MAT file also contains
a `README` variable recording the summed pixel solid angle and the independently
integrated field of view used for validation.

## `MSIlluminanceToAGCKernel.mat`

The output of
[`defineMSIlluminanceToAGCKernel.py`](../code/defineWorldCameraCalibration/defineMSIlluminanceToAGCKernel.py).
The MATLAB struct `MSIlluminanceToAGCKernel` is the standalone temporal-filter
contract loaded by `defineMSIlluminanceToAGCLag.py`. Its fields are:

- `sourceFile` and `definition`: kernel provenance;
- `timeSeconds` and `normalizedWeights`: the exact timebase and unit-sum
  weights passed to convolution;
- `rawMeanKernelSum`, `sampleIntervalSeconds`, `durationSeconds`, and
  `sampleCount`: normalization and sampling properties;
- `stepResponseTauFitsSeconds`: descriptive exponential time constants for
  the two source step responses (the empirical weights, not these fitted
  exponentials, are used for convolution); and
- `convolutionPadding` and `convolutionMode`: the boundary and output rules.

Regenerate this file before the lag artifact whenever the historical source
kernel changes.

## `MSIlluminanceToAGCLag.mat`

The output of
[`defineMSIlluminanceToAGCLag.py`](../code/defineWorldCameraCalibration/defineMSIlluminanceToAGCLag.py). The
MATLAB struct `MSIlluminanceToAGCLag` contains the shared camera-response lag in seconds
and included / skipped recording counts. A fresh raw-recording run also
populates the candidate lag grid and mean correlation at each candidate; those
diagnostic arrays are empty in the migrated legacy result until it is rerun.
The `model` field records that the transformation is empirical FIR smoothing
followed by a pure lag shift. The nested `kernel` struct contains the exact
normalized weights, timebase, sampling properties, fitted source-response time
constants, provenance, and convolution rules loaded from the artifact named by
`kernelArtifact`. Positive lag means that the camera response follows the minispect signal. The
empirical AGC/illuminance data-prep script reads this file before exporting its
calibration point cloud. Regenerate it by passing the raw `GKA` recording paths
directly to `defineMSIlluminanceToAGCLag.py`.

## `cameraScoreToAverageLuminance.mat`

The output of
[`defineAGCToMeanLuminance.m`](../code/defineWorldCameraCalibration/defineAGCToMeanLuminance.m).
The top-level vectors `cameraScore` and `avgSceneLuminance` contain the fixed
integrating-sphere calibration points. Interpolating between them in log10
space maps analog-gain × digital-gain × exposure to average scene luminance in
cd/m². The file also contains a top-level `README` describing this contract.
