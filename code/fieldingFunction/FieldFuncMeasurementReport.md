**Flat Fielding Function Report and Data Collection/Calibration Notes** 

*Planetarium Measurement Procedure (Friday, June 26, 2025)* 

The camera and tripod were positioned near the approximate center of the Fels Planetarium at the Franklin Institute. 
During data collection, the planetarium dome displayed a projected image that appeared spatially uniform in spectral 
content.

The camera was oriented toward the zenith (the center of the dome) and remained fixed in tilt throughout the recording. 
To sample the projected illumination across different portions of the dome, the camera was slowly rotated clockwise 
about its optical axis. Four complete 360° rotations were performed over the course of the recording.

The first approximately four minutes of recording were collected with no camera movement. This stationary period was 
included to allow the camera's automatic gain control (AGC) to stabilize and adapt to the ambient luminance conditions 
within the planetarium.

Following this stabilization period, the camera was rotated through a full revolution approximately once per minute. 
Rotations were performed in discrete increments of roughly 45°, with each step occurring every 6–7 seconds. One exception 
was the first rotation, which required approximately 1.5 minutes to complete and was performed using less regular angular 
increments. At the start of each rotation, the north side of the planetarium corresponded approximately to the bottom of 
the camera field of view.

The camera's pointing direction remained centered on the zenith throughout the recording; only rotational orientation 
about the optical axis was changed between positions.


Additional Calibration Notes (Geoff):

I obtained three measurements of the luminance of the dome. At the zenith (center of the camera image) the luminance was 
25.63 cd/m2. Two measurements (corresponding to the initial center bottom of the camera image and to the center left side 
of the initial image) obtained at ~20° elevation from the horizon (so, 70° eccentricity in the camera image) were 22.56 and 22.64. 
This is an average change of 3 cd/m2 over 70°. I unfortunately did not obtain a measurement in between, so we will have to assume 
a linear fall-off of intensity with eccentricity of about 0.0428 cd/m2/deg as we move away from the zenith.


*Generation of the Average Flat-Field Image* 

Selected video frames were loaded into MATLAB and converted to floating-point representation. To gather an ideal represention of
rotated images, snippets of (approximately one minute-long) video footage were assessed seperately, as well as a collection of
frames from throughout the video. The mean image was then computed on a pixel-by-pixel basis across all selected frames for each
section, all of which produced near-identical results. We observed a horizontal displacement of the hot spot (toward the left side
of the frame) but this was consistent across all of the observed individual images across the rotation of the camera. As such, this 
spatial structure was attributed to the lens of the camera and not to the planetarium dome. One minute-long snippet was selected as
a reliable representation of the data collected and kept for subsequent analysis.

The resulting average image represented the spatial sensitivity pattern of the camera-lens system. Surface plots of the 
averaged image were generated to visualize large-scale intensity variations across the field of view.


*Separation of Bayer Color Channels*

The camera sensor uses a Bayer color filter array. To determine whether the spatial sensitivity pattern differed across color channels, 
the averaged raw data were separated into red, green, and blue pixel classes according to the sensor Bayer pattern (BGGR). Individual 
red, green, and blue channel images were constructed by extracting the corresponding pixel locations from each frame and averaging 
them separately across all selected frames.


*Hotspot Model Fitting*

Inspection of the average image revealed a broad hotspot centered near the optical axis. To obtain a smooth estimate of the spatial 
sensitivity pattern, a two-dimensional parametric model was fit separately to the red, green, and blue channel images.

The model consisted of a Gaussian-like spatial profile with an additional hyperbolic tangent compression term:

$I(x,y)=B+Atanh(kexp[−(2σx2​(x−x0​)2​+2σy2​(y−y0​)2​)])$

where B is the baseline intensity, A is the amplitude, (x0,y0) is the hotspot center, σx and σy describe the hotspot width, and k 
controls the degree of flattening near the peak. Parameters were estimated using nonlinear least-squares optimization (lsqcurvefit). 
Separate fits were performed for the red, green, and blue channel images.

To evaluate the quality of the fits, horizontal, vertical, and diagonal cross-sections through the hotspot center were compared 
between the measured data and the fitted model. The profiles were normalized using the minimum and maximum values of each fitted 
surface to facilitate comparison across color channels. The resulting profiles showed substantial overlap among the red, green, 
and blue channels, indicating minimal chromatic dependence of the hotspot structure.


*Construction of the Fielding Correction Map*

The fitted hotspot surfaces were used to derive multiplicative flat-field correction factors. For each color channel, a correction 
factor was calculated at every pixel location as:

$C(x,y)=I(x,y)Imax​​$

where Imax is the peak value of the fitted hotspot surface and I(x,y) is the fitted sensitivity at the corresponding pixel location. 
This procedure produces correction factors equal to one at the hotspot peak and greater than one in less-sensitive regions of the 
image.

Separate correction maps were generated for the red, green, and blue Bayer channels and then recombined according to the sensor 
Bayer pattern to create a full-resolution correction map. Verification was performed by multiplying each fitted hotspot surface by 
its corresponding correction map, which produced an approximately uniform image as expected.