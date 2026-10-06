Data collected using the IMX219 camera are passed through a pre-processing pipeline. The routines in this directory were used to characterize the properties of the camera needed to create this pipeline.

All programmatically generated MAT files in `derived/` use the shared
`utilities/saveDerivedFile.m` contract, which requires a substantive top-level
`README` variable and verifies the completed artifact. Python-derived MAT files
use the equivalent `derived_io.save_derived_mat` helper. `defineAllScript.m`
validates the entire derived directory before completing.

We have configured the IMX219 to save RAW, 8-bit images (640x480). Camera sensitivity is under the control of a custom automatic gain control that operates upon the mean sensor value across each image frame, adjusting analog and digital gain, and exposure time. The camera chip is behind a fisheye lens, which produces both spatial variation in image intensity and image distortion. The chip has R, G, and B sensors, which differ in their radiometric sensitivity.

The processing pipeline adjusts the data to account for these properties of the measurement. In order of application, we characterize:

defineFullWellCapacityEffect -- The camera chip reports sensor values that are linearly related to photon capture within each photodiode well. There is a roll-off in this linear function, however, at high sensor values which we attribute to a "full well capacity" effect. To characterize this effect, we measured camera sensor values across an ~1.4 log unit range. This routine analyzes these measurements, and demonstrates that the sensor values are well described by a saturating exponential function. This function is then used to implement a a linearization function that transforms raw sensor values to a linearized form. In practice, this adjusts the original value range of 16-254 to the linearized range of 0-254 (removing the dark value).

defineFlatFieldingFunction -- We assume that the camera chip itself has uniform spatial sensitivity to light. The lens which directs light to the chip, however, introduces spatial variation in light intensity across the chip surface. A "flat fielding function" is used to correct for this effect. We characterized the fielding function for our camera by directing the camera towards the zenith of a uniformly illuminated, 60' diameter hemisphere (the Fels Planetarium at the Franklin Institute). The camera was rotated about its optical axis while we recorded images. The images were linearized, and then averaged across rotations; this step distributes any spatial non-uniformity in the light source across the resulting data. This routine loads these images and then fits a modified Gaussian to the 2D surface of sensor values separately for the R, G, and B channels. Wavelength variation in transverse chromatic aberration may be appreciated in these fitting functions. The resulting functions are used to correct (flatten) the images in the pre-processing pipeline.

defineRadiometricCorrection -- The R, G, and B channels on the camera chip are created by placing each pixel behind a transmittance filter that limits the spectral sensitivity of that class of photodiode. The sensor values reported by the chip across the channel classes do not necessarily reflect the relative radiometric intensity of light falling upon the sensor for each class. To correct for this, we measured the SPD of a cloudy sky using the PR670, and obtained images of this same sky using the IMX219 world camera. These data to produce the derived file radiometricCorrectionRGB.mat.

defineAGCToMeanLuminance -- The camera chip reports 8 bit sensor values. These values reflect variation in radiance across the scene around the set-point of the automatic gain control. Before linearization, a sensor value of 127 represents the mid-point of the sensor range and is the set-point target of the AGC. The linearization step maps 127 --> 57. We wish to express the sensor values in absolute radiometric units, instead of relative sensor values. In this routine we characterize how AGC settings are related to the mean luminance of the field to which the camera is exposed. With knowledge of this relationship, we can use the AGC settings to identify the set-point of the camera in units of luminance, and then interpret each pixel as absolute luminance using: AGC mean luminance * (pixel value / 57). The resulting lookup vectors are saved in `derived/cameraScoreToAverageLuminance.mat`.

defineMSIlluminanceToAGCKernel -- Converts the historical AGC simulation output
in `data/agc_empirical_kernels.mat` into the standalone derived contract
`derived/MSIlluminanceToAGCKernel.mat`.

defineMSIlluminanceToAGCLag -- The runnable Python temporal-alignment stage. It
loads `derived/MSIlluminanceToAGCKernel.mat`, applies that kernel to raw minispect signals, selects one shared lag across
recordings, and writes the lag together with the complete kernel definition to
`derived/MSIlluminanceToAGCLag.mat`. The empirical AGC and
illuminance data-prep script reads this derived lag when building the MATLAB
calibration point cloud.

defineFisheyeCameraIntrinsics -- The world camera is a ArduCam B0392 IMX219 Wide Angle M12. This is a wide-angle, fisheye lens system. This routine works upon a set of images taken with the camera to derive the file: arducamB0392cameraIntrinsics.mat

defineCameraToVisualAngles -- Uses the calibrated fisheye intrinsics to map every
640-by-480 world-camera pixel center to a unit viewing direction. It reproduces
the finite-difference solid-angle calculation in
`world_util.world_frame_visual_angle_to_steradians`, validates the summed pixel
areas against an independent integration of the calibrated camera field of
view, and saves the 480-by-640 `deltaSteradians` array to
`derived/cameraToVisualAngles`. Each element gives the solid angle represented
by the corresponding raw world-camera pixel, in steradians.


## Python calibration constants

The calibration block in `code/library/sensor_utility/world_util.py` loads the
same derived files as MATLAB `pipeline/reconstructionPipeline.m`. Constants are
loaded once when the module is imported. Restart Python, or reload the module,
after regenerating the calibration files.

| Constants | Source and purpose |
| --- | --- |
| `WORLD_DARK_SIGNAL` | `derived/darkSignal.mat`: dark offset subtracted from the quantization-corrected 8-bit counts. |
| `WORLD_FULL_WELL_CLIPPING_EXPONENT` | `derived/nonLinearClippingExponent.mat`: fitted exponent used to invert the sensor's full-well response. |
| `WORLD_MAX_ALLOWED_LINEARIZATION_DERIVATIVE` | Matches MATLAB's value of 4.0. Samples at or above the derived raw-count threshold become `Inf` and are imputed. |
| `WORLD_FIELDING_FUNCTIONS` | `derived/flatFieldingFunction.mat`: spatial sensitivity correction maps, indexed by `(rows, cols)`. |
| `WORLD_RADIOMETRIC_CALIBRATION` | Contents of `derived/radiometricCorrectionRGB.mat`, loaded once for the channel weights and full map. |
| `WORLD_RGB_SCALARS` | Red, green, and blue calibration weights from `radiometricCorrectionRGB`. |
| `WORLD_RADIOMETRIC_CORRECTION_MAP` | MATLAB's full Bayer correction map, used to calculate the harmonic mean for normalization. |
| `WORLD_BAYER_CORRECTION_MATRICES` | Read-only Python correction maps built with `generate_RGB_mask` and the RGB weights. The BGGR pattern has blue at even/even indices and red at odd/odd indices. |
| `WORLD_INTEGRATED_RADIANCE_CALIBRATION` | Contents of `derived/cameraScoreToIntegratedRadiance.mat`. |
| `WORLD_AGC_TO_RADIANCE_P` | Polynomial coefficients mapping log10 camera score to log10 mean integrated radiance. |
| `WORLD_MEAN_SPATIAL_CORRECTIONS` | Product of the harmonic means of the fielding and RGB maps. This scales the linearized AGC set point. |
| `WORLD_TRIANGULATION_SHEAR` | Small affine coordinate shear used by SciPy interpolation to resolve ties on the regular Bayer grid. It is an implementation setting, not a measured calibration. |

The current reconstruction adds 0.375 counts before dark subtraction, preserves
negative linearized counts, and imputes non-positive or infinite samples before
fielding and RGB correction. Digital gain affects the AGC set point rather than
multiplying the image. Camera score is exposure times analog gain divided by the
spatially corrected, linearized set point. The fit then supplies integrated
radiance in W/m²/sr. See `pipeline/reconstructionPipeline.m` for the reference
calculation and `definePipelineParameters/defineAGCToIntegratedRadianceViaMacbeth.m`
for the current radiance calibration.

The affine shear preserves interpolation weights within a given triangle, but
MATLAB and SciPy can still choose different triangles around missing samples
and make different boundary choices. See
`validation/validatePythonReconstructionPipeline.ipynb` for measured agreement
and its limits; the shear does not guarantee exact parity for every image.


### Saving and displaying Python reconstruction stages

Pass an empty list to collect independent snapshots, including the untouched
raw counts followed by all six processing stages:

```python
stages: list[np.ndarray] = []
radiance = world_util.world_transformation_pipeline(
    raw_frame, agc_settings, stagewise_output=stages
)
figure, axes = world_util.plot_reconstruction_stages(stages)
```

Omitting `stagewise_output` avoids the snapshot copies. Buffer input produces
one full buffer per saved stage; select the same frame from each stage before
plotting. `visualize_pipeline(recording_path, chunk_range=(0, 1))` selects the
middle frame of a recording chunk, runs this same pipeline, and plots the saved
stages. Its two-row layout follows MATLAB's `plotReconstructionStages`: images
above RGB histograms, ceiling samples in red, floor samples in blue, and original
values in the image-coordinate readout.

The in-place reconstruction helpers (`linearize_camera_responsivity`,
`apply_fielding_function`, `apply_color_correction`, and
`world_counts_to_radiance`) return `None`. They still show before/after figures
when `visualize_results=True`. Imputation creates a separate output and therefore
returns that new array.
