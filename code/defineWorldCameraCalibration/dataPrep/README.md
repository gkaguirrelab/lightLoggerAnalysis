Description of routines used to prep raw data files and populate the "data" directory of the lightLoggerAnalysis repo

deriveEmpircalAGCAndIlluminance.py -- Runnable Python derivation for the empirical AGC calibration input. With no arguments, it selects the calibrated camera-score and minispect-illuminance point cloud from the checked-in processing cache and writes `data/empircalAGC.mat`. The output contains the MATLAB struct `empiralAGC` with the fields `cameraScoreLinear` and `msIlluminance`. Optional raw `GKA` recording paths may be supplied to rebuild the processing cache before export.
