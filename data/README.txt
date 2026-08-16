Directory contents:

- ASM7341_spectralSensitivity.mat: spectral sensitivity functions for the first 10 channels of the "minispect" chip. These values were taken from the spreadsheet "AS7341_Filter_Templates.xlsx" which was supplied by the manufacturer. The values in this table correspond to Figures 18 and 19 of the white paper AMS Datasheet DS000504 "AS7341; 11-Channel Multi-Spectral Digital Sensor; v3-00 • 2020-Jun-25"

- TSL2591_spectralSensitivity.mat: spectral sensitivity functions for the first 2 channels of the TSL2591 chip. These values were generated via web digitizing figure 12 from the TSL2591 datasheet, then feeding these csv files into generate_spectral_sensitivity_plot.py. This returns a dataframe, which was then converted to a MATLAB table.

- IMX219_spectralSensitivity.mat: spectral sensitivity functions for the three channels (RGB) of the IMX camera chip. These values (provided by the first author) correspond to Figure 18 of the paper Pagnutti 2017 J. Electron. Imaging 26(1), 013014. The SPDs were then scaled by the radiometric calibration coefficients that we measured in the lab, with the scalar weights being (RGB): [1.0 0.83 1.06].

- D65_SPD.mat: spectral power distribution for the CIE D65 illuminant standard as a table variable with wavelengths and relative power. Downloaded from:  http://files.cie.co.at/204.xls

- CIEDaylightComponents_T.mat: basis set for reconstruction of daylight illuminant. https://cie.co.at/datatable/components-relative-spectral-distribution-daylight

- camera_linearity_ND0_ND0p4_rgb_means.mat: full-well capacity effect calibration data. We measured the CombiLED with the IR filter in the optical pathway, made measurements at 0 and 0.4 NDF with AGC set for 0.4 NDF, 0.5 contrast, and 127 value, then measured both NDFs across a variety of CombiLED setting levels. We fit a correction with fitFullWellCapacityEffect.m and derived an exponent from this. This is the data used to do that fitting.

- fieldingFunctions/: fielding correction functions derived from the planetarium recording post linearization.

- camera_agc_illuminance_linear_scale.mat: The output of an analysis of empirical light logger recordings. Recordings from multiple tasks / subjects were analyzed to find the low-pass temporal filter that best related environmental illuminance (as recorded by the minispect) to camera sensitivity (the product of the AGC settings). These data are plotted by the routine defineAGCToMeanLuminance.m to be compared to measurements made within the light integrating sphere. We would like to replace this data file with a call to a routine that processes the files.

- CloudySkySPD_37degSolarElevation.mat: Geoff used the PR670 to make a measurement of the cloudy sky outside of Goddard Labs when the sun was at 37° elevation. Separately, the IMX219 camera was used to record a video of the sky (Zach is looking for these raw images). These measurements are used in the calculation of the radiometric correction of the RGB channels in the function defineRadiometricWeights.

- fisheyeLensCalibration: Directory that contains the images and "session file" used by the matlab camera calibration app to derive the intrinsics for the world camera lens. Additional details available in the readme file within this directory.