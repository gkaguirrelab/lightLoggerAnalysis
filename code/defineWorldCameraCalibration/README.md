Data collected using the IMX219 camera are passed through a pre-processing pipeline. The routines in this directory were used to characterize the properties of the camera needed to create this pipeline.

We have configured the IMX219 to save RAW, 8-bit images (640x480). Camera sensitivity is under the control of a custom automatic gain control that operates upon the mean sensor value across each image frame, adjusting analog and digital gain, and exposure time. The camera chip is behind a fisheye lens, which produces both spatial variation in image intensity and image distortion. The chip has R, G, and B sensors, which differ in their radiometric sensitivity.

The processing pipeline adjusts the data to account for these properties of the measurement. In order of application, we characterize:

defineFullWellCapacityEffect -- The camera chip reports sensor values that are linearly related to photon capture within each photodiode well. There is a roll-off in this linear function, however, at high sensor values which we attribute to a "full well capacity" effect. To characterize this effect, we measured camera sensor values across an ~1.4 log unit range. This routine analyzes these measurements, and demonstrates that the sensor values are well described by a saturating exponential function. This function is then used to implement a a linearization function that transforms raw sensor values to a linearized form. In practice, this adjusts the original value range of 16-254 to the linearized range of 0-254 (removing the dark value).

defineFlatFieldingFunction -- We assume that the camera chip itself has uniform spatial sensitivity to light. The lens which directs light to the chip, however, introduces spatial variation in light intensity across the chip surface. A "flat fielding function" is used to correct for this effect. We characterized the fielding function for our camera by directing the camera towards the zenith of a uniformly illuminated, 60' diameter hemisphere (the Fels Planetarium at the Franklin Institute). The camera was rotated about its optical axis while we recorded images. The images were linearized, and then averaged across rotations; this step distributes any spatial non-uniformity in the light source across the resulting data. This routine loads these images and then fits a modified Gaussian to the 2D surface of sensor values separately for the R, G, and B channels. Wavelength variation in transverse chromatic aberration may be appreciated in these fitting functions. The resulting functions are used to correct (flatten) the images in the pre-processing pipeline.

defineRadiometricCorrection -- The R, G, and B channels on the camera chip are created by placing each pixel behind a transmittance filter that limits the spectral sensitivity of that class of photodiode. The sensor values reported by the chip across the channel classes do not necessarily reflect the relative radiometric intensity of light falling upon the sensor for each class. To correct for this, we measured the SPD of a cloudy sky using the PR670, and obtained images of this same sky using the IMX219 world camera. These data to produce the derived file radiometricCorrectionRGB.mat.

defineAGCToMeanLuminance -- The camera chip reports 8 bit sensor values. These values reflect variation in radiance across the scene around the set-point of the automatic gain control. Before linearization, a sensor value of 127 represents the mid-point of the sensor range and is the set-point target of the AGC. The linearization step maps 127 --> 57. We wish to express the sensor values in absolute radiometric units, instead of relative sensor values. In this routine we characterize how AGC settings are related to the mean luminance of the field to which the camera is exposed. With knowledge of this relationship, we can use the AGC settings to identify the set-point of the camera in units of luminance, and then interpret each pixel as absolute luminance using: AGC mean luminance * (pixel value / 57).

deriveAGCLag -- The runnable Python temporal-alignment stage. It applies the
empirical AGC kernel to raw minispect signals, selects one shared lag across
recordings, and writes `derived/cameraAGCLag.mat`. The empirical AGC and
illuminance data-prep script reads this derived lag when building the MATLAB
calibration point cloud.

utilities/fitCameraAGCToIlluminance -- Fits the selected empirical camera-score
and illuminance point cloud from `data/empircalAGCAndIlluminance.mat` with a
continuous two-slope model in log10 space and an L1 objective. It writes the
coefficients, fitting ranges, sample count, shared lag, and provenance to the
`cameraAGCToIlluminanceFit` struct in
`derived/cameraAGCToIlluminanceFit.mat`. `defineAllScript.m` invokes this
function after the other world-camera definitions. Python's
`video_io.camera_scores_to_illuminance` loads this artifact directly, so the
camera fit is not repeated during video processing.

utilities/fitEmpircalAGCtoIlluminance -- Runs the shared fit and produces the
diagnostic plot comparing the fitted curve, empirical recordings, breakpoint,
and integrating-sphere measurements. Because it calls the shared fitter, it
also updates `derived/cameraAGCToIlluminanceFit.mat`; it does not produce a
separate diagnostic MAT file. See `utilities/README.md` for the complete input,
output, and downstream-consumer contract.

defineFisheyeCameraIntrinsics -- The world camera is a ArduCam B0392 IMX219 Wide Angle M12. This is a wide-angle, fisheye lens system. This routine works upon a set of images taken with the camera to derive the file: arducamB0392cameraIntrinsics.mat

defineDeltaSteradians -- Uses the calibrated fisheye intrinsics to map every
640-by-480 world-camera pixel center to a unit viewing direction. It reproduces
the finite-difference solid-angle calculation in
`world_util.world_frame_visual_angle_to_steradians`, validates the summed pixel
areas against an independent integration of the calibrated camera field of
view, and saves the 480-by-640 `deltaSteradians` array to
`derived/deltaSteradians.mat`. Each element gives the solid angle represented
by the corresponding raw world-camera pixel, in steradians.
