# Camera Chromatic Verification Measurement

## Overview

This measurement was performed in the **Flicker Room** to verify the chromatic response of the Light Logger camera against a Macbeth ColorChecker, using a PR-670 spectroradiometer as the reference instrument.

## Equipment

| Item | Details |
| --- | --- |
| Spectroradiometer | Photo Research PR-670, unit #2 |
| Color target | Macbeth ColorChecker |
| Camera | Light Logger N02, camera set 2 |

## Experimental Setup

The Macbeth ColorChecker was positioned upright against a box in the Flicker Room. The PR-670 was placed in front of the chart and aligned with the selected color patches. The Light Logger camera was positioned nearby to capture the same target for comparison.

![Experimental setup showing the Light Logger camera, Macbeth ColorChecker, and PR-670 spectroradiometer](setup.jpeg)

*Figure 1. Experimental setup in the Flicker Room.*

## Selected ColorChecker Patches

Five patches were selected for measurement. Patch locations use **one-indexed row and column coordinates**:

| Patch | Row | Column |
| :---: | :---: | :---: |
| R1C3 | 1 | 3 |
| R2C3 | 2 | 3 |
| R1C6 | 1 | 6 |
| R2C6 | 2 | 6 |
| R4C4 | 4 | 4 |

![Macbeth ColorChecker with the five selected patches marked in red](selectedColors.jpeg)

*Figure 2. Selected ColorChecker patches, marked in red.*

## Procedure

1. Position the Macbeth ColorChecker upright in the Flicker Room, supported by a box.
2. Place the PR-670 in front of the ColorChecker and align it with the target patch.
3. Run `measureUncontrolledSourceSpectrum.m` to acquire **five measurements** from each selected patch.
4. Repeat the measurement sequence for all five patches listed above.
5. Record from various distances with the world camera and MS on the light logger 
6. Extract suitable world camera frames

## Objective

The objective of this procedure was to verify the chromatic properties of the Light Logger camera by comparing its response to the known ColorChecker patches and the corresponding PR-670 spectral measurements.
