# Dark Noise Measurements

This directory contains raw dark-frame measurements for the world camera at
fixed automatic gain control (AGC) operating states. The measurements are used
to estimate the camera dark signal and dark noise separately for each AGC
state, because the effective camera noise depends on analog gain, digital gain,
and exposure.

Each `AGCstate*` folder contains a stack of TIFF dark frames acquired with the
camera covered with both the lens cap and a black shroud, and the entire room in darkness. The frames are linearly spaced through the recording. These frames should be processed within each AGC state, not
pooled across AGC states.

## Dropbox Location of Raw Data 

Dropbox location:

```text
FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/darkNoiseCalibrations
```

## Recorded Camera Settings

The settings below were read from every row of the world-camera metadata in
the raw dark-noise recordings. The six metadata columns, in order, are
`timestamp`, `cameraAgain`, `AGCDgain`, `cameraExposure`, `AGCAgain`, and
`AGCExposure`. The table reports the three settings that were actually applied
by the camera: `cameraAgain`, `AGCDgain`, and `cameraExposure`.

| Folder | AGC state | `cameraAgain` | `AGCDgain` | `cameraExposure` |  |
| --- | ---: | ---: | ---: | ---: | :---: |
| `AGCstate1` | 1 | 1 | 1 | 39 |  |
| `AGCstate2` | 2 | 1 | 1 | 4180 |  |
| `AGCstate3` | 3 | 1 | 1 | 8290 |  |
| `AGCstate4` | 4 | 5.5652174949645996 | 1 | 8290 |  |
| `AGCstate5` | 5 | 10.666666984558105 | 1 | 8290 |

The applied camera settings were fixed within each AGC-state recording.
