# Eye Cleaning Explorer

Use this page to look through FLIC eye recordings, examine default performance, and try different cleaning
settings. It runs the existing `cleanEyeTracking.py` code and shows its results
in the browser. Changing a setting does not overwrite the recording's data,
PDF, or metadata. The screenshots and demo below use synthetic data, not participant recordings.

![Recording selection, cleaning settings, and the six comparison plots](media/overview.png)

## Run the viewer

1. Open a terminal in the `lightLoggerAnalysis` repository folder and activate the Python environment you use for the analysis pipeline. Copy and paste this command, then press **Return**:

   ```sh
   python code/analyzeEyeTracking/webCode/server.py
   ```

2. Leave the terminal running. When it prints `FLIC eye viewer: http://127.0.0.1:8765`, open your browser and paste this into the address bar:

   ```text
   http://127.0.0.1:8765
   ```

3. Click **Choose folder…** and select your data folder as described below.

To stop the viewer, press **Control+C** in the terminal. To restart it, stop it first, run the same command again, and refresh the browser page.

## Select the data folder

Click **Choose folder…** and select the experiment folder directly containing
the subject folders. For example:

```text
scriptedIndoorOutdoor_2026_photophobia/
└── FLIC_0021/
    └── chat/
        └── 1/
            └── eyeData/
```

Select `scriptedIndoorOutdoor_2026_photophobia` in this example, not `FLIC_0021`
or `eyeData`.

The viewer needs these four files in each `eyeData` directory:

| File | Used for |
| --- | --- |
| `3d_eye_states.csv` | Pupil diameter and eyelid measurements. |
| `blinks.csv` | Blink start and end times. |
| `alternative_camera_gaze.csv` | Gaze mapped into world-camera coordinates. |
| `eye_cleaning_metadata.pkl` | The settings saved when the recording was cleaned. |

For Dropbox folders, make sure the files are downloaded to this computer.
The viewer does not need a NAS connection once the input files are local.

## 1. Select Recording

Choose a **Subject**, **Activity**, and **Recording number** in the left panel.
The plot loads automatically. Selecting another recording also loads that
recording's saved cleaning settings; adjustments do not carry over.

## 2. Examine or Set Settings

The initial values come from `eye_cleaning_metadata.pkl`. These are the settings
used to generate the recording's cleaned data.

Edit a value and click **Apply & replot**. Until you apply it, the plot still
shows the previous result. Click **Restore saved settings** to discard your
changes and recalculate with the original values. Scroll within the left panel
to reach the remaining settings and buttons.

| Setting | What it controls |
| --- | --- |
| Pupil diameter | Measurements below **Low** or above **High** are rejected. Values exactly at either bound are allowed. |
| Sudden eye changes | How much pupil diameter or mean eyelid aperture can change between consecutive samples before being flagged as a spike. A change must be greater than the threshold; equality is allowed. |
| Gaze position bounds | The largest allowed X and Y coordinates in pixels. Coordinates from zero through the maximum are allowed, including both endpoints. |
| Sudden gaze movement | The distance gaze moves between consecutive samples. Distances greater than this threshold are flagged. The cleaner combines this with blink detection when rejecting gaze. |
| Blink padding | How far to extend each blink before its start and after its end. Units are nanoseconds: 100,000,000 ns is 100 ms. Both extended endpoints are included. |
| Fill limits | How long a missing stretch can be before the cleaner leaves it unfilled. Limits use sample-index distance, not seconds. Interior gaps include the bounding valid samples; edge gaps measure from the edge to the valid sample. A distance equal to the limit is allowed. |
| Mask width | The window used to extend rejected regions to neighboring samples. It must be a positive whole number. |
| Smoothing width | The number of samples in the smoothing window. It must be positive and odd, such as 1, 3, or 5. |

### Example: change a pupil bound

This step-by-step demo changes **Pupil diameter → High** from 8.5 to 5 mm,
applies the change, then restores the saved settings. The example values are
for demonstrating the controls, not recommended cleaning settings.

![Demo: edit the pupil upper bound, apply it, and restore saved settings](media/adjust-settings.gif)

The tall pupil excursion in the example is rejected after applying the lower
upper bound. The screenshot below shows the recalculated result. The vertical
axis rescales automatically, so check its numbers when comparing views.

![Recalculated plots after changing the pupil upper bound to 5 mm](media/settings-applied.png)

## 3. Visualize Plots

The six panels show pupil diameter, mean eyelid aperture, left-minus-right
aperture difference, gaze X, gaze Y, and the original blink events.

Black points are the cleaner's raw reference and red lines are the cleaned
signals. The raw reference has already been aligned and passed through the
cleaner's basic bounds checks; it is not an untouched plot of every CSV row.
The last panel shows the original blink intervals without the added padding.
Time is in seconds from the first gaze sample overlapping the eye-state data.

Use the toolbar below the figure:

- **Home:** restore the original plot limits. This does not restore parameter values.
- **Pan:** click the crossed-arrows button, then drag within a panel.
- **Zoom:** click the magnifying glass, then drag a rectangle around a region.
- **Back / Forward:** move through the toolbar's saved views.
- **Download:** choose a format and click the save icon to download the current figure.

You can also scroll over a plot to zoom its time axis. The time axes are linked,
so all six panels stay aligned. Scroll outside the plot to move through the
page. Browser zoom can make the compact interface larger if needed.
