# RGB Radiometric-Correction Measurements

This directory contains the paired cloudy-sky camera and PR670 measurements used by `defineRadiometricWeights.m` to derive the IMX219 RGB radiometric weights.

The raw world-camera chunks are reconstructed as a gap-filled temporary video with digital gain, response linearization, and color weighting disabled. Ten consecutive raw-Bayer grayscale frames are selected from the lossless temporary video:

| Output TIFF | Zero-based temporary-video frame index |
| --- | ---: |
| `rawFrames/0.tiff` | 8000 |
| `rawFrames/1.tiff` | 8001 |
| `rawFrames/2.tiff` | 8002 |
| `rawFrames/3.tiff` | 8003 |
| `rawFrames/4.tiff` | 8004 |
| `rawFrames/5.tiff` | 8005 |
| `rawFrames/6.tiff` | 8006 |
| `rawFrames/7.tiff` | 8007 |
| `rawFrames/8.tiff` | 8008 |
| `rawFrames/9.tiff` | 8009 |

Because the video reconstruction fills timestamp gaps, these are temporary-video indices, not necessarily indices counting only physical frames in the raw chunks.

`CloudySkySPD_37degSolarElevation.mat` is the primary PR670 spectral measurement and cannot be recreated from the camera chunks. `cropExample.tiff` is an illustrative crop retained with the measurement inputs.
