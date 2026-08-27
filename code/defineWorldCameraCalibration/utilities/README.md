# World Camera Calibration Utilities

This directory contains shared helpers used by the world-camera calibration
definitions.

`saveDerivedFile.m` is the required writer for MATLAB artifacts placed in the
project's `derived/` directory. It requires a substantive top-level `README`,
saves the requested variables through a temporary file, and atomically replaces
the destination only after verifying the payload.

`validateDerivedFilesHaveREADME.m` checks every MAT artifact in `derived/` and
fails if any file lacks valid README metadata. `defineAllScript.m` runs this
validation after generating the complete calibration set.

The other utilities split and combine Bayer frames, return Bayer-channel
indices, and apply the IMX219 sensor-count linearization used by the individual
definition routines.
