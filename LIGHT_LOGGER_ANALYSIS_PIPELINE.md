# Light Logger Analysis Pipeline

This document summarizes the recording-to-analysis workflow implemented in `code/preprocessRecordingData/preprocessing_pipeline.py` and the MATLAB/Python routines it calls. The pipeline is interacted with through the associated Jupyter notebook, [`code/preprocessRecordingData/data_preprocessing_pipeline.ipynb`](code/preprocessRecordingData/data_preprocessing_pipeline.ipynb). The Python module is a collection of callable pipeline functions and its `main()` function is empty, so the notebook—not direct execution of `preprocessing_pipeline.py`—is the operational entry point.

## How the pipeline is run

`data_preprocessing_pipeline.ipynb` presents the stages in workflow order and provides the editable calls used to run them. The notebook imports `preprocessing_pipeline`, commonly calls `importlib.reload(preprocessing_pipeline)` before a stage so local Python changes are picked up without restarting the kernel, and then invokes the corresponding function with explicit subject/activity filters, overwrite behavior, verbosity, color mode, and other stage-specific options. This makes the notebook the interactive control surface and run log, while `preprocessing_pipeline.py` contains the reusable implementation. Long-running cells can be executed one stage at a time, inspected before continuing, and selectively rerun without automatically repeating the complete pipeline.

<small>**Primary interaction files:** [`code/preprocessRecordingData/data_preprocessing_pipeline.ipynb`](code/preprocessRecordingData/data_preprocessing_pipeline.ipynb), [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py)</small>

## Required Python environments

All notebook stages **except the Egocentric Video Mapper stage** must be run in the lab's modified version of the `pylids` environment. This is the notebook's normal kernel and supplies the scientific Python, MATLAB-engine, video-I/O, and project-specific dependencies used by data transfer, validation, world-video generation, virtual foveation, SPD analysis, aggregation, and actigraphy.

The Egocentric Video Mapper stage must instead be run using the separate `egocentric_video_mapper` environment/kernel. That environment contains the mapper's dedicated Pupil Labs package, feature-matching and optical-flow dependencies, model support, and GPU-oriented stack. After the mapper outputs have been generated, switch the notebook kernel back to the modified `pylids` environment before continuing with tag/task segmentation and all later stages. Variables and imports should be recreated by rerunning the notebook's setup/import cell after each kernel change.

| Notebook work | Required environment/kernel |
|---|---|
| Stages 1–6 and 8–15 | Lab-modified `pylids` environment |
| Stage 7: Egocentric Video Mapper | `egocentric_video_mapper` environment |

<small>**Relevant files:** [`code/preprocessRecordingData/data_preprocessing_pipeline.ipynb`](code/preprocessRecordingData/data_preprocessing_pipeline.ipynb), [`code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/README.rst`](code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/README.rst), [`code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/pyproject.toml`](code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/pyproject.toml)</small>

## Storage locations and naming conventions

The current defaults divide the data into three tiers:

| Tier | Default location | Purpose |
|---|---|---|
| Acquisition media | `/Volumes/T7 Shield` and Pupil Cloud | Temporary light-logger storage and remote Neon recordings |
| Study metadata | `/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_subject/NEWscriptedIndoorOutdoorVideos2026` | Subject-summary and activity-definition spreadsheets |
| Raw NAS | `/Volumes/FLIC_raw/NEWscriptedIndoorOutdoorVideos2026` | Immutable-ish source chunks and exported Neon time series |
| Processing NAS | `/Volumes/FLIC_processing/NEWscriptedIndoorOutdoorVideos2026` | Playable world videos, gaze mappings, segmentation files, and derived videos |
| Analysis output | `.../FLIC_analysis/lightLogger/NEWscriptedIndoorOutdoorVideos2026` | Per-recording SPD results, aggregate `.mat` files, and figures |

Recording names are expected to follow `FLIC_<subject number>_<activity>_<recording number>`. Every recording number is retained as a separate attempt beneath the corresponding sensor directory. Subject and activity include/exclude arguments use `_is_desired_item`: a nonempty `*_to_process`/`*_to_transfer` collection is an allowlist; otherwise, everything except the `*_to_skip` collection is included. The subject-summary and activity spreadsheets that provide the expected cohort and task names during integrity checking are stored in `/Users/zacharykelly/Aguirre-Brainard Lab Dropbox/Zachary Kelly/FLIC_subject/NEWscriptedIndoorOutdoorVideos2026`. `get_subject_ids()` selects the spreadsheet whose filename contains `subjectsummarymigraine`, while `get_activity_names()` selects the spreadsheet whose filename contains `activitiesmigraine`; both filename searches are case-insensitive.

The normal raw and processed hierarchy is:

```text
FLIC_<subject>/
└── <activity>/
    ├── GKA/
    │   └── <recording number>/      # raw chunks, or W.avi in processing
    ├── Neon/
    │   └── <recording number>/      # raw Neon export or processed mapping output
    ├── tag_task_start_end.mat       # processing tier
    ├── ..._task_justProjection.avi  # processing tier
    └── ..._task_virtuallyFoveated.avi
```

## Numerically ordered pipeline

In normal use, each numbered stage below corresponds to a section or executable cell group in `data_preprocessing_pipeline.ipynb`. The notebook calls into `preprocessing_pipeline.py`, which in turn dispatches to the lower-level Python and MATLAB files listed for that stage.

## 1. Acquire and name the two synchronized recording streams

Each session begins with two independently stored streams. The GKA light logger writes a directory of chunked world-camera and multispectral-sensor arrays to the acquisition SSD. The Neon glasses recording is uploaded to Pupil Cloud under the matching `FLIC_<subject>_<activity>_<recording number>` name. The repository does not implement either device's recording process or the upload to Pupil Cloud; `preprocessing_pipeline.py` begins once the GKA directory is visible on the SSD and the Neon recording is already present in the cloud. `analyze_rerecording_frequency()` can audit the cloud naming/re-recording history before selection.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py)</small>

## 2. Transfer the latest GKA recording from the SSD to the raw NAS

`transfer_light_logger_recordings()` parses every SSD directory name into subject, activity, and recording number; retains every recording number for each subject/activity; applies the requested allow/skip filters; and copies each source directory to `FLIC_raw/.../FLIC_<subject>/<activity>/GKA/<recording_number>`. Existing numbered destinations are skipped unless `overwrite_existing=True`. This establishes the raw-NAS side of the canonical hierarchy without transforming the recorded chunk data. `rename_world_recordings()` is not part of this normal path: it is a legacy helper for older source trees whose light-logger directory was not already named `GKA`.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py)</small>

## 3. Select and download the latest Neon recording from Pupil Cloud

`download_pupil_cloud_recordings()` authenticates with the local Pupil Cloud API key, lists recordings in the selected workspace, parses the same subject/activity/recording-number convention, and retains every numbered attempt. It requests each raw-data export as `Timeseries Data + Scene Video.zip`, streams it in 1 MiB chunks, and retries blocked NAS opens or writes after a delay. Each archive is extracted independently into `<activity>/Neon/<recording_number>/`, then removed after successful extraction. A typical numbered result contains `enrichment_info.txt`, `sections.csv`, one recording subdirectory, the Neon scene video, and CSV tables for gaze, fixations, blinks, saccades, events, IMU, eye states, and timestamps. The companion `get_pupil_cloud_apikey()` reads the secret from the expected untracked file `code/apiTokens/pupil_cloud.txt`.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py)</small>

## 4. Validate raw GKA and Neon data before processing

`verify_data_integrity()` compares the NAS contents against the expected subjects and activities, checks that both `GKA` and `Neon` exist and are nonempty, and validates the light-logger `config.pkl`. By default it expects world and multispectral sensors, a 120 Hz world camera, a 640-by-480 frame, and custom AGC with target 127. It also reads chunk metadata, supports both legacy (`Again`, `Dgain`, `exposure`) and modern (`cameraAgain`, `AGCDgain`, `cameraExposure`) schemas, and can plot analog gain, digital gain, and exposure by activity. Neon validation checks required export files and a scene `.mp4`; `display_integrity_issues()` formats the nested problem report. The lighter `verify_neon_integrity()` provides warning-only Neon checks.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/library/sensor_utility/world_util.py`](code/library/sensor_utility/world_util.py), [`code/library/matlabIO/python_libraries/video_io.py`](code/library/matlabIO/python_libraries/video_io.py)</small>

## 5. Convert raw world-camera chunks into a corrected playable video

`generate_world_videos()` calls `video_io.world_chunks_to_video()` for every numerical raw `GKA/<recording_number>` directory and writes the corresponding `FLIC_processing/.../<activity>/GKA/<recording_number>/W.avi`. Existing videos are skipped independently, so one completed attempt does not prevent other attempts from being generated. The converter reads naturally ordered world metadata/frame pairs through a two-buffer shared-memory loader and obtains frame geometry and FPS from `config.pkl`. The transformations are applied in this order: (1) convert the raw Bayer frames to `float64` and mark raw saturated samples with an infinity sentinel while retaining raw zero-valued floor samples; (2) linearize camera responsivity using the dark-noise offset and nonlinear full-well correction, then reapply the saved saturation mask; (3) apply the measured spatial flat-field correction; (4) apply the per-Bayer-site RGB radiometric weights; (5) multiply each frame by its metadata-derived digital gain; (6) de-Bayer the corrected image to RGB while preserving a copy of the corrected raw Bayer contributors and carrying saturation sentinels through interpolation; (7) use the de-Bayer provenance map to force any RGB output pixel with a raw floor or ceiling contributor to the output floor or ceiling; (8) resolve remaining sentinels, round and clip values to `[0, 255]`, convert to `uint8`, and convert RGB to BGR for OpenCV; and (9) interleave black placeholder frames wherever the timestamp-derived frame distribution indicates a dropped frame, then write the real and placeholder frames to the FFV1 AVI at the configured camera rate. Timestamp values are converted from nanoseconds to seconds when the gap distribution is calculated.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/library/matlabIO/python_libraries/video_io.py`](code/library/matlabIO/python_libraries/video_io.py), [`code/library/sensor_utility/world_util.py`](code/library/sensor_utility/world_util.py), [`derived/flatFieldingFunction.mat`](derived/flatFieldingFunction.mat), [`derived/nonLinearClippingExponent.mat`](derived/nonLinearClippingExponent.mat), [`derived/radiometricCorrectionRGB.mat`](derived/radiometricCorrectionRGB.mat)</small>

## 6. Review world/Neon pairing when synchronization is uncertain (QC)

`verify_world_neon_pairing()` is a manual visual check that displays the first GKA frame beside the first non-gray Neon scene frame. It confirms that the two streams placed under a subject/activity are the intended pair before gaze-dependent processing continues. This is a QC branch rather than a data-producing dependency.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/library/matlabIO/python_libraries/video_io.py`](code/library/matlabIO/python_libraries/video_io.py)</small>

## 7. Map Neon gaze and fixation data onto the light-logger world video

`generate_egocentric_mapper_results()` pairs the raw Neon time-series directory with the corrected `W.avi`. Its wrapper converts the AVI to a temporary MP4, then runs the bundled Pupil Labs Egocentric Video Mapper with the requested gaze/fixation mode, optical-flow method, feature matcher, and refresh threshold. The mapper aligns the Neon scene video to the alternative GKA camera and produces gaze/fixation coordinates and alternative-camera timestamps in Neon time. Results are copied from a local temporary directory to `FLIC_processing/.../<activity>/Neon/egocentric_mapper_results/`; downstream stages principally consume `alternative_camera_gaze.csv` and `alternative_camera_timestamps.csv`.

**Environment:** Switch the notebook to the `egocentric_video_mapper` environment/kernel for this stage only. Once mapping is complete, switch back to the lab-modified `pylids` environment and rerun the notebook setup/import cell before proceeding.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/applyVirtualFoveation/pythonCode/virtual_foveation.py`](code/applyVirtualFoveation/pythonCode/virtual_foveation.py), [`code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/__main__.py`](code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/__main__.py), [`code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/gaze_mapper.py`](code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/gaze_mapper.py), [`code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/sync_videos.py`](code/applyVirtualFoveation/pythonCode/egocentric_video_mapper/src/pupil_labs/egocentric_video_mapper/sync_videos.py)</small>

## 8. Divide each recording into AprilTag and task ranges

`generate_tag_task_start_ends()` inspects the corrected video's FPS and frame count and writes `tag_task_start_end.mat`. For recordings longer than the default four-minute task duration, the last four minutes are labeled `task` and all earlier frames are labeled `tag`; for shorter recordings, the entire video is used as the task and `[1, 1]` is used as a placeholder tag range. Indices are stored one-based for direct MATLAB use. This stage supplies the exact frame ranges used by virtual foveation.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/library/matlabIO/python_libraries/video_io.py`](code/library/matlabIO/python_libraries/video_io.py)</small>

## 9. Generate retinally projected and virtually foveated videos

`generate_virtually_foveated_videos()` starts MATLAB and calls `generateVirtuallyFoveatedVideos()` for each selected tag/task range. The MATLAB driver loads `W.avi`, alternative-camera timestamps and gaze coordinates, raw Neon blink events, the range from `tag_task_start_end.mat`, and the calibrated fisheye intrinsics. Gaze pixels are converted to angular coordinates and paired to world frames by nearest timestamp. `virtuallyFoveateVideo()` processes frames in cached chunks and parallelizes the per-frame work: invalid gaze, angles beyond 45 degrees, blinks, or empty frames become blank output frames; valid frames are passed to `virtuallyFoveateFrame()`, which rotates/orients the image, maps fisheye pixels into degrees, and interpolates a 480-by-480 retinal grid spanning ±60 degrees around gaze. The `justProjection` control holds gaze at its median, while `virtuallyFoveated` follows gaze. The normal products are a foveated tag video plus foveated and projection-only task videos.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/applyVirtualFoveation/matlabCode/generateVirtuallyFoveatedVideos.m`](code/applyVirtualFoveation/matlabCode/generateVirtuallyFoveatedVideos.m), [`code/applyVirtualFoveation/matlabCode/virtuallyFoveateVideo.m`](code/applyVirtualFoveation/matlabCode/virtuallyFoveateVideo.m), [`code/applyVirtualFoveation/matlabCode/virtuallyFoveateFrame.m`](code/applyVirtualFoveation/matlabCode/virtuallyFoveateFrame.m), [`code/applyVirtualFoveation/matlabCode/anglesFromIntrinsics.m`](code/applyVirtualFoveation/matlabCode/anglesFromIntrinsics.m), [`code/library/matlabIO/@videoIOWrapper/videoIOWrapper.m`](code/library/matlabIO/@videoIOWrapper/videoIOWrapper.m), [`derived/arducamB0392cameraIntrinsics.mat`](derived/arducamB0392cameraIntrinsics.mat)</small>

## 10. Validate the derived-video set (QC)

`verify_virtually_foveated_video_integrity()` checks that each activity has the expected derived AVI set (normally at least three files: tag foveated, task foveated, and task projection-only; gaze calibration expects fewer) and checks the task products against the target duration. This gate catches missing or prematurely terminated MATLAB renders before spectral analysis begins.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/library/matlabIO/python_libraries/video_io.py`](code/library/matlabIO/python_libraries/video_io.py)</small>

## 11. Compute per-recording temporal spectral power distributions

`generate_spds()` locally copies each processed activity to avoid long NAS reads and calls MATLAB once for every requested projection type and color mode. `mapSPDs()` reads the task AVI in overlapping one-second windows stepped every 0.5 seconds. Within each temporal window it scans overlapping 24-by-24-pixel spatial patches at a 12-pixel step, converts RGB into the requested signal (`L+M+S`, `L-M`, `GRAY`, `a`, `c_lm`, or `c_s`), treats black/ceiling frames as missing, and calculates a temporal PSD for each patch. A log-log straight-line fit supplies the spectral exponent; the integrated fitted spectrum supplies a variance estimate. Results are averaged over time and saved as `..._<projection type>_SPDResults.mat`, containing the exponent map, variance map, SPD by region, frequencies, median image, and frame-drop vector; companion PDFs show the maps and regional spectra.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/analyzeSPD/matlabCode/processSPDs.m`](code/analyzeSPD/matlabCode/processSPDs.m), [`code/analyzeSPD/matlabCode/mapSPDs.m`](code/analyzeSPD/matlabCode/mapSPDs.m), [`code/analyzeSPD/matlabCode/calcTemporalSPD.m`](code/analyzeSPD/matlabCode/calcTemporalSPD.m), [`code/analyzeSPD/matlabCode/plotSPDs.m`](code/analyzeSPD/matlabCode/plotSPDs.m), [`code/library/matlabIO/@videoIOWrapper/readFrame.m`](code/library/matlabIO/@videoIOWrapper/readFrame.m)</small>

## 12. Standardize SPD axes and compare projection modes (optional plotting)

The `_find_spd_axes_*()` helpers load the per-recording `.mat` results and calculate common display limits across all data, per subject, or per activity. Map bounds are restricted to the 480-by-480 elliptical analyzed field, zeros are excluded, the lower map bound uses the 15th percentile, exponent signs are reversed for display, and frequency/regional-SPD bounds are forced to valid positive log-axis values. `adjust_spd_axes()` then regenerates figures with these common limits and can place `justProjection` beside/overlaid with `virtuallyFoveated`; `adjust_spd_axes_copy()` is a parallel legacy plotting path.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/analyzeSPD/matlabCode/plotSPDs.m`](code/analyzeSPD/matlabCode/plotSPDs.m), [`code/analyzeSPD/matlabCode/plotSPDs_copy.m`](code/analyzeSPD/matlabCode/plotSPDs_copy.m)</small>

## 13. Aggregate SPD results across subjects, activity groups, or all activities

The first aggregation level, `generate_spds_across_subject()`, calls `processSPDsAcrossSubjects()` to stack per-recording maps and regional spectra for a given activity and projection type, average them across subjects, preserve participant-level region averages, and write `..._SPDResultsAcrossSubjects.mat` plus figures. `generate_spds_across_groups()` then calls `processSPDAcrossActivities()` for configured sets such as indoor (`work`, `chat`, `walkIndoor`) and outdoor (`walkOutdoor`, `walkBiopond`, `sitBiopond`). `generate_spds_across_all()` performs the same operation over every retained activity. These are alternative summaries of the same per-recording results, rather than serial transforms that must all be run.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/analyzeSPD/matlabCode/processSPDsAcrossSubjects.m`](code/analyzeSPD/matlabCode/processSPDsAcrossSubjects.m), [`code/analyzeSPD/matlabCode/processSPDAcrossActivities.m`](code/analyzeSPD/matlabCode/processSPDAcrossActivities.m), [`code/analyzeSPD/matlabCode/plotSPDs.m`](code/analyzeSPD/matlabCode/plotSPDs.m)</small>

## 14. Assemble grouped panels and statistical mean products (optional reporting)

`group_spds_per_subject()` and `group_spds_across_subjects()` organize indoor/outdoor activities into tiled exponent, variance, and regional-SPD figures using `groupSPDs()`, either for one participant or from the across-subject result files. The newer statistical path uses `spd_util.load_spds()` plus `generate_mean_spds()`/`meanSPDs()` to calculate element-wise means and standard deviations across subjects or activities, including center/periphery best-fit curves. `combine_spds()` overlays multiple color modes and projection types for each subject/activity, and `plot_mean_spds()` plots the across-subject mean products. The `activity_then_subject` mean path and two `plot_mean_spds()` dimensions remain explicitly unimplemented.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/analyzeSPD/pythonCode/spd_util.py`](code/analyzeSPD/pythonCode/spd_util.py), [`code/analyzeSPD/matlabCode/groupSPDs.m`](code/analyzeSPD/matlabCode/groupSPDs.m), [`code/analyzeSPD/matlabCode/meanSPDs.m`](code/analyzeSPD/matlabCode/meanSPDs.m), [`code/analyzeSPD/matlabCode/combineSPDs.m`](code/analyzeSPD/matlabCode/combineSPDs.m), [`code/analyzeSPD/matlabCode/plotMeanSPDs.m`](code/analyzeSPD/matlabCode/plotMeanSPDs.m)</small>

## 15. Generate actigraphy, eye-state, and light-state summaries (parallel analysis branch)

`generate_actigraphy_graphs()` branches after Neon mapping because it needs raw Neon `imu.csv`, `3d_eye_states.csv`, `blinks.csv`, and `gaze.csv`, plus the alternative-camera timestamps created by the mapper. `plotParticipantState()` aligns those signals in Neon time, computes ENMO activity and head rotational velocity, masks blinks, summarizes gaze and eye aperture/pupil diameter, classifies active/inactive periods, and exports the participant-state and correlation figures. Its optional luminance path also uses raw multispectral/world chunks, the tag/task ranges, the known world-to-AS offset, and sensor calibration to classify indoor/outdoor periods.

<small>**Relevant files:** [`code/preprocessRecordingData/preprocessing_pipeline.py`](code/preprocessRecordingData/preprocessing_pipeline.py), [`code/analyzeActigraphy/plotParticipantState.m`](code/analyzeActigraphy/plotParticipantState.m), [`code/analyzeActigraphy/classifyActiveInactivePeriods.m`](code/analyzeActigraphy/classifyActiveInactivePeriods.m), [`code/analyzeMS/matlabCode/classifyIndoorOutdoorPeriods.m`](code/analyzeMS/matlabCode/classifyIndoorOutdoorPeriods.m), [`code/analyzeMS/matlabCode/msCounts2Illuminance.m`](code/analyzeMS/matlabCode/msCounts2Illuminance.m), [`code/library/sensor_utility/ms_util.py`](code/library/sensor_utility/ms_util.py)</small>

## Final products

The core path produces four main artifact families:

1. A canonical raw recording on the NAS containing untouched GKA chunks and an unpacked Neon export.
2. A corrected, playable `W.avi` plus gaze/fixation mappings expressed in the GKA camera coordinate system and Neon time.
3. Gaze-centered and projection-only task videos with matching temporal ranges.
4. Per-recording and optional aggregate SPD `.mat` files and publication-ready PDFs, with actigraphy/eye-state summaries as a parallel branch.
