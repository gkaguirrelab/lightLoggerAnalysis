# Dark Noise Measurements

This directory contains raw dark-frame measurements for the world camera at
fixed automatic gain control (AGC) operating states. The measurements are used
to estimate the camera dark signal and dark noise separately for each AGC
state, because the effective camera noise depends on analog gain, digital gain,
and exposure.

Each `AGCstate*` folder contains a stack of TIFF dark frames acquired with the
camera covered. These frames should be processed within each AGC state, not
pooled across AGC states.

## Dropbox Location of Raw Data 

Dropbox location:

```text
FLIC_admin/Equipment/ArduCam B0392 IMX219 Wide Angle M12/darkNoiseCalibrations/AGCstate1
```

## AGC States

The AGC states used for these measurements are:

| Folder | AGC state | Analog gain | Digital gain | Exposure |
| --- | ---: | ---: | ---: | ---: |
| `AGCstate0` | 0 | 1.00000000 | 1.00000000 | 1466 |
| `AGCstate1` | 1 | 1.80281687 | 1.00000000 | 8333 |
| `AGCstate2` | 2 | 10.66600000 | 1.58156674 | 8333 |
| `AGCstate3` | 3 | 10.66600000 | 5.96331605 | 8333 |
| `AGCstate4` | 4 | 10.66600000 | 7.97561560 | 8333 |
| `AGCstate5` | 5 | 10.66600000 | 10.00000000 | 8333 |

