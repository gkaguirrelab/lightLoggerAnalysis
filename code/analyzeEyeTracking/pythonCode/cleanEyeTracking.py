"""Clean one Neon recording, or run directly to select its folders."""

import pathlib
import warnings
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter1d


def clean_eye_data(
    raw_neon_path: str,
    processing_neon_path: str,
    verbose: bool = False,
    pupil_diameter_thresholds: tuple[float, float] = (1.5, 8.5),
    pupil_change_threshold: float = 0.25,
    aperture_change_threshold: float = 0.3,
    plot_zoom_seconds: float = 2.0,
) -> np.ndarray:
    """Clean bilateral eye measurements and align mapped gaze to eye samples.

    Apply pupil validity bounds, blink/spike masking, linear gap filling,
    and smoothing. Combine available left/right eye measurements and
    align gaze with the eye-state timeline using linear interpolation.

    Args:
        raw_neon_path (str): Folder containing 3d_eye_states.csv and blinks.csv.
        processing_neon_path (str): Folder containing
            alternative_camera_gaze.csv for the same recording.
        verbose (bool): Display each transformation as a row of before (black)
            and after (red) plots, ending with raw/final comparisons.
            Defaults to False.
        pupil_diameter_thresholds (tuple[float, float]): Inclusive minimum and
            maximum pupil diameters in millimeters. Defaults to (1.5, 8.5).
            Each left/right pupil measurement outside these limits is replaced
            with NaN before averaging the eyes. Eyelid aperture is unaffected
            by these limits.
        pupil_change_threshold (float): Absolute change in combined pupil
            diameter between consecutive samples, in millimeters, above which
            a sample is flagged as a spike. Defaults to 0.25.
        aperture_change_threshold (float): Absolute change in combined eyelid
            aperture between consecutive samples, in millimeters, above which
            a sample is flagged as a spike. Defaults to 0.3. Either threshold
            being exceeded triggers spike masking for all signals.
        plot_zoom_seconds (float): Width of each verbose plot's detail window
            in seconds. Defaults to 2.0; must be positive and finite.

    Returns:
        np.ndarray: Float array of shape (N, 5), with one row per eye-state
            sample. Columns are elapsed time [s], combined pupil diameter
            [mm], combined eyelid aperture [mm], gaze x [px], and gaze y [px].
            Elapsed time starts at the first eye-state timestamp. Unfilled
            missing values remain NaN. This function does not write files.

    Raises:
        ValueError: Either timeline has fewer than two samples, timestamps
            are not strictly increasing, or a blink ends before it starts.

    Notes:
        Processing follows Geoff's October 20, 2026 final Gemini script:
        https://upennvisionlabs.slack.com/archives/D076YU3T22U/p1792522033305849
        Cleaning pads blinks by 100 ms before and 150 ms after, then dilates
        by five samples on each side. Verbose plots show each cleaning stage.
        All input timestamps must use the same clock and nanosecond units.
        Eye and gaze samples need not have identical timestamps or coverage;
        different sampling rates are supported but are not required. Camera
        frame timestamps are a separate timeline and are not used here.
    """
    if verbose and (not np.isfinite(plot_zoom_seconds) or plot_zoom_seconds <= 0):
        raise ValueError("plot_zoom_seconds must be positive and finite.")

    # 1. Load the data files.    
    eye_columns: list[str] = [
            "pupil diameter left [mm]", "pupil diameter right [mm]",
            "eyelid aperture left [mm]", "eyelid aperture right [mm]",
    ]
    eye_data: pd.DataFrame = pd.read_csv(
        pathlib.Path(raw_neon_path) / "3d_eye_states.csv",
        usecols=["timestamp [ns]", *eye_columns],
    )

    gaze_columns: list[str] = ["gaze x [px]", "gaze y [px]"]
    gaze_data: pd.DataFrame = pd.read_csv(
        pathlib.Path(processing_neon_path) / "alternative_camera_gaze.csv",
        usecols=["timestamp [ns]", *gaze_columns],
    )
    blink_data: pd.DataFrame = pd.read_csv(
        pathlib.Path(raw_neon_path) / "blinks.csv",
        usecols=["start timestamp [ns]", "end timestamp [ns]"],
    )

    # 2. Extract timestamps on the shared recording clock.
    # Eye states and mapped gaze can share sample times but have different lengths.
    # In the checked FLIC_20/chat export, both are ~200 Hz; mapped gaze ends earlier.
    # alternative_camera_timestamps.csv describes video frames (~120 Hz), not gaze.
    timestamps: np.ndarray = eye_data["timestamp [ns]"].to_numpy(dtype=float)
    gaze_timestamps: np.ndarray = gaze_data["timestamp [ns]"].to_numpy(dtype=float)

    # Validate the timestamps. Must be more than 1 sample to do diff and they must be monotonic 
    if len(timestamps) < 2 or len(gaze_timestamps) < 2:
        raise ValueError("Eye and gaze exports must each contain at least two samples.")
    if np.any(np.diff(timestamps) <= 0) or np.any(np.diff(gaze_timestamps) <= 0):
        raise ValueError("Eye and gaze timestamps must be strictly increasing and unique.")

    # Match MATLAB readtable's default double timestamps, including epoch-time
    # rounding in interpolation and blink-boundary comparisons.
    # Gaze may start before the first eye sample; retain that negative offset.
    time: np.ndarray = (timestamps - timestamps[0]) / 1e9
    gaze_time: np.ndarray = (gaze_timestamps - timestamps[0]) / 1e9

    # eyes is the working four-column array: left/right pupil, then left/right aperture.
    eyes: np.ndarray = eye_data[eye_columns].to_numpy(dtype=float)

    # Keep a separate copy because pupil validity filtering below modifies eyes.
    # This copy supplies the original black dots in the comparison plots.
    raw_eyes: np.ndarray = eyes.copy()

    # raw_gaze already contains camera-mapped x/y, but still at gaze timestamps.
    # It stays unchanged; aligned_gaze below has one row per eye-state sample.
    raw_gaze: np.ndarray = gaze_data[gaze_columns].to_numpy(dtype=float)



    # 1. Read in the data 
    # 2. put gaze angles onto a common time space with the other things (just make sure this is true, the sampling rate will be different regardless)
    #       - do we want to interpolate gaze angles outside of 120hz? We're going to go back to 120hz anyway
    # 3. NaN out blinks from all of these things because now they are on a common timescale
    #      - maybe padding? add a little bit to the start and end of the blinks 
    # 4. Filter pupil diameter from the low and high thresholds
    # 5. combine the eye measurements (after from both eyes) for pupil diameter and aperature 
    #       - if we have 1 eye open 1 eye closed, what does this count as? Blink? 
    #       - if we have 1 nan 1 not, what do we do with this? (think we use the remaining one, makes sense for pupil diameter)


    # 3. Put gaze on the eye timeline so all signals can use the same blink masks.
    # Matching timestamps return the original values; uncovered times remain NaN.
    aligned_gaze: np.ndarray = np.empty((len(timestamps), 2))
    for column in range(2):
        aligned_gaze[:, column] = np.interp(
            timestamps, gaze_timestamps, raw_gaze[:, column], left=np.nan, right=np.nan
        )

    # Do not invent a continuous gaze trace across a long source dropout.
    max_source_gap_ns: float = 15 * float(np.median(np.diff(timestamps)))
    gap_indices: np.ndarray = np.flatnonzero(np.diff(gaze_timestamps) > max_source_gap_ns)
    gap_starts: np.ndarray = np.searchsorted(
        timestamps, gaze_timestamps[gap_indices], side="right"
    )
    gap_ends: np.ndarray = np.searchsorted(
        timestamps, gaze_timestamps[gap_indices + 1], side="left"
    )
    for start, end in zip(gap_starts, gap_ends):
        aligned_gaze[start:end] = np.nan

    # 4. Use the recorded blink intervals for cleaning.
    blink_starts: np.ndarray = blink_data["start timestamp [ns]"].to_numpy(dtype=float)
    blink_ends: np.ndarray = blink_data["end timestamp [ns]"].to_numpy(dtype=float)

    if np.any(blink_ends < blink_starts):
        raise ValueError("Blink end timestamps must not precede their starts.")

    # Geoff's cleaning mask includes asymmetric padding.
    padded_starts: np.ndarray = np.searchsorted(
        timestamps, blink_starts - 100_000_000, side="left"
    )
    padded_ends: np.ndarray = np.searchsorted(
        timestamps, blink_ends + 150_000_000, side="right"
    )
    padded_changes: np.ndarray = np.zeros(len(timestamps) + 1, dtype=np.int64)
    np.add.at(padded_changes, padded_starts, 1)
    np.add.at(padded_changes, padded_ends, -1)
    padded_blinks: np.ndarray = np.cumsum(padded_changes[:-1]) > 0

    # 5. Apply pupil bounds independently, then average the available eyes.
    # pupils is a view into eyes, so these assignments also update eyes.
    pupils: np.ndarray = eyes[:, :2]
    minimum_pupil_diameter: float = pupil_diameter_thresholds[0]
    maximum_pupil_diameter: float = pupil_diameter_thresholds[1]
    pupils[(pupils < minimum_pupil_diameter) | (pupils > maximum_pupil_diameter)] = np.nan

    def average_available_eyes(eye_measurements: np.ndarray) -> np.ndarray:
        """Average left/right values per sample; keep NaN if both are missing."""
        available_eye_count: np.ndarray = np.sum(~np.isnan(eye_measurements), axis=1)
        measurement_sum: np.ndarray = np.nansum(eye_measurements, axis=1)
        has_measurement: np.ndarray = available_eye_count > 0

        # Start with missing values, then fill samples with at least one eye.
        average: np.ndarray = np.full(len(eye_measurements), np.nan)
        average[has_measurement] = (
            measurement_sum[has_measurement] / available_eye_count[has_measurement]
        )
        return average

    # Combine cleaned measurements and prepare the raw measurements for plotting.
    pupil: np.ndarray = average_available_eyes(eyes[:, :2])
    aperture: np.ndarray = average_available_eyes(eyes[:, 2:])
    raw_pupil: np.ndarray = average_available_eyes(raw_eyes[:, :2])
    raw_aperture: np.ndarray = average_available_eyes(raw_eyes[:, 2:])

    # signals collects pupil, aperture, gaze x, and gaze y on the same timeline.
    signals: np.ndarray = np.column_stack((pupil, aperture, aligned_gaze))

    # 6. Detect spikes, expand spike rejection, and censor all signals.
    # Compute changes before masking so lid transitions can still trigger rejection.
    pupil_changes: np.ndarray = np.r_[0.0, np.abs(np.diff(pupil))]
    aperture_changes: np.ndarray = np.r_[0.0, np.abs(np.diff(aperture))]
    is_spike: np.ndarray = (
        (pupil_changes > pupil_change_threshold)
        | (aperture_changes > aperture_change_threshold)
    )

    # Expand blink and spike masks independently by five samples on either side.
    blink_dilated: np.ndarray = maximum_filter1d(padded_blinks, 11, mode="constant", cval=0)
    spike_dilated: np.ndarray = maximum_filter1d(is_spike, 11, mode="constant", cval=0)

    if verbose:
        before_masking: np.ndarray = signals.copy()
    signals[blink_dilated | spike_dilated] = np.nan

    # 7. Fill general gaps up to 50 samples and blink gaps up to 400 samples.
    # These limits also apply at recording edges: short uncovered tails may be
    # extrapolated, but a long tail (such as the example's ~24 s) remains missing.
    # The smoothing step can extend values a few samples into a missing region.
    # imputed holds gap-filled data; cleaned below holds the final smoothed data.
    imputed: np.ndarray = signals.copy()
    for column in range(4):
        values: np.ndarray = signals[:, column]
        valid: np.ndarray = np.flatnonzero(~np.isnan(values))
        # MATLAB linear fillmissing needs two valid samples, even at edges.
        if valid.size < 2:
            continue

        # Find each consecutive NaN run, including runs touching recording boundaries.
        missing: np.ndarray = np.isnan(values)
        transitions: np.ndarray = np.diff(np.r_[False, missing, False].astype(np.int8))
        starts = np.flatnonzero(transitions == 1)
        ends = np.flatnonzero(transitions == -1)
        lengths: np.ndarray = ends - starts

        # MATLAB MaxGap counts the distance between valid bounding samples.
        gap_sizes: np.ndarray = lengths + ((starts > 0) & (ends < len(values)))
        missing_indices: np.ndarray = np.flatnonzero(missing)

        # Fill short gaps everywhere, then allow longer fills only within blink masks.
        for max_gap in (50, 400):
            eligible: np.ndarray = np.repeat(gap_sizes <= max_gap, lengths)
            indices: np.ndarray = missing_indices[eligible]
            if max_gap == 400:
                indices = indices[blink_dilated[indices]]

            imputed[indices, column] = np.interp(indices, valid, values[valid])

            # Match MATLAB's linear extrapolation for eligible edge gaps.
            if valid.size >= 2:
                for first, second, edge_indices in (
                    (valid[0], valid[1], indices[indices < valid[0]]),
                    (valid[-1], valid[-2], indices[indices > valid[-1]]),
                ):
                    slope: float = float((values[second] - values[first]) / (second - first))
                    imputed[edge_indices, column] = values[first] + (edge_indices - first) * slope

    # 8. Smooth states with a median and gaze with quadratic Savitzky–Golay.
    cleaned: np.ndarray = np.full_like(imputed, np.nan)
    offsets: np.ndarray = np.arange(-2, 3, dtype=float)
    position_bits: np.ndarray = 1 << np.arange(5)

    for column in range(4):
        values = imputed[:, column]
        windows: np.ndarray = np.lib.stride_tricks.sliding_window_view(
            np.pad(values, (2, 2), constant_values=np.nan), 5
        )

        # NaN padding makes the median ignore nonexistent samples at the edges.
        if column < 2:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
                cleaned[:, column] = np.nanmedian(windows, axis=1)
            continue

        # Encode which of the five samples are valid, then reuse one fit per pattern.
        # This avoids solving a separate least-squares problem at every gaze sample.
        patterns: np.ndarray = (~np.isnan(windows)).astype(np.uint8) @ position_bits
        for pattern in np.unique(patterns):
            if pattern == 0:
                continue

            # Fewer than three valid points require a lower-degree fit.
            present: np.ndarray = (pattern & position_bits) != 0
            degree: int = min(2, int(np.count_nonzero(present)) - 1)
            design: np.ndarray = np.vander(offsets[present], degree + 1, increasing=True)
            rows: np.ndarray = np.flatnonzero(patterns == pattern)
            # The constant coefficient is the fitted value at the window center.
            cleaned[rows, column] = windows[rows][:, present] @ np.linalg.pinv(design)[0]

        # MATLAB uses full one-sided windows for the first/last two gaze samples.
        width: int = min(5, len(values))
        edge_indices: list[int] = sorted(
            set(range(min(2, len(values))))
            | set(range(max(0, len(values) - 2), len(values)))
        )

        for index in edge_indices:
            start = 0 if index < 2 else len(values) - width
            edge: np.ndarray = values[start:start + width]
            present = ~np.isnan(edge)
            if not present.any():
                continue

            degree = min(2, int(np.count_nonzero(present)) - 1)
            edge_offsets: np.ndarray = np.arange(start, start + width, dtype=float) - index
            design = np.vander(edge_offsets[present], degree + 1, increasing=True)
            cleaned[index, column] = np.linalg.pinv(design)[0] @ edge[present]

    if verbose:
        names = ("Pupil (combined)", "Eyelid aperture (combined)", "Gaze x", "Gaze y")
        units = ("Diameter [mm]", "Aperture [mm]", "Position [px]", "Position [px]")
        raw = (raw_pupil, raw_aperture, raw_gaze[:, 0], raw_gaze[:, 1])

        # Each panel describes its own input timeline, so original gaze samples
        # remain at their actual timestamps. Two-column inputs show both eyes.
        stages = [
            ("Gaze alignment to eye timestamps (with source-gap masking)", [
                (names[c + 2], units[c + 2], gaze_time, raw_gaze[:, c], aligned_gaze[:, c])
                for c in range(2)
            ]),
            ("Pupil diameter bounds", [
                (f"Pupil ({side})", units[0], time, raw_eyes[:, c], eyes[:, c])
                for c, side in enumerate(("left", "right"))
            ]),
            ("Average available left/right eyes", [
                (names[c], units[c], time, eyes[:, 2*c:2*c+2], before_masking[:, c])
                for c in range(2)
            ]),
            ("Blink and spike masking (including padding and dilation)", [
                (names[c], units[c], time, before_masking[:, c], signals[:, c])
                for c in range(4)
            ]),
            ("Linear gap filling (general and blink gaps)", [
                (names[c], units[c], time, signals[:, c], imputed[:, c])
                for c in range(4)
            ]),
            ("Smoothing (median for eye states; Savitzky–Golay for gaze)", [
                (names[c], units[c], time, imputed[:, c], cleaned[:, c])
                for c in range(4)
            ]),
            ("Final comparison to raw data (available-eye averages for eye states)", [
                (names[c], units[c], time if c < 2 else gaze_time, raw[c], cleaned[:, c])
                for c in range(4)
            ]),
        ]

        figure = plt.figure(figsize=(20, 3 * len(stages)), constrained_layout=True)
        rows = figure.subfigures(nrows=len(stages), ncols=1)
        shared_axis = None
        for row_index, (row, (title, panels)) in enumerate(zip(rows, stages)):
            row.suptitle(title, fontweight="bold")
            axes = row.subplots(1, len(panels), squeeze=False)[0]
            final_row = row_index == len(stages) - 1
            for axis, (name, unit, input_time, before, after) in zip(axes, panels):
                if shared_axis is None:
                    shared_axis = axis
                else:
                    axis.sharex(shared_axis)
                before_lines = axis.plot(input_time, before, "k.", markersize=2)
                # Both eyes are previous-stage inputs with the same styling;
                # use one legend handle even when two black traces are drawn.
                before_lines[0].set_label("Raw" if final_row else "Previous")
                axis.plot(time, after, "r.-", markersize=1.5, linewidth=0.5,
                          label="Final" if final_row else "Result")
                axis.set(title=name, ylabel=unit, xlabel="Time [s]")
                axis.legend(fontsize="small", loc="upper left",
                            bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
                axis.grid(alpha=0.25)

                # Score changes on the output timeline while leaving the plotted
                # input on its original timestamps. For bilateral input, compare
                # each eye to the result and keep the larger change per sample.
                input_columns = before[:, None] if before.ndim == 1 else before
                aligned_before = np.column_stack([
                    np.interp(time, input_time, values, left=np.nan, right=np.nan)
                    for values in input_columns.T
                ])
                finite_before = np.isfinite(aligned_before)
                finite_after = np.isfinite(after)[:, None]
                comparable = finite_before & finite_after
                changes = np.zeros_like(aligned_before)
                np.subtract(aligned_before, after[:, None], out=changes, where=comparable)
                changes = np.abs(changes)
                # Removing or filling a sample is a change too. Give these events
                # priority so masking/imputation rows zoom into a changed gap.
                missing_changed = finite_before != finite_after
                changes[missing_changed] = max(1.0, float(changes.max()) * 2)
                scores = changes.max(axis=1)
                if np.any(scores > 0):
                    focus_index = int(np.argmax(scores))
                else:
                    # An unchanged signal has no transformation event to rank.
                    # Show its middle available sample (or the recording middle).
                    available = np.flatnonzero(finite_after[:, 0] | finite_before.any(axis=1))
                    focus_index = int(available[len(available) // 2]) if available.size else len(time) // 2

                # Keep a small, independently scaled detail view inside each
                # overview. Its time limits never join the shared overview axes.
                duration = time[-1] - time[0]
                zoom_width = min(plot_zoom_seconds, duration)
                zoom_start = np.clip(time[focus_index] - zoom_width / 2,
                                     time[0], time[-1] - zoom_width)
                zoom_stop = zoom_start + zoom_width
                detail = axis.inset_axes([0.58, 0.58, 0.38, 0.34])
                detail.set_facecolor("white")
                detail.patch.set_alpha(0.95)
                # Include one neighboring sample on either side to retain lines
                # that cross the detail window without a sample inside its edge.
                input_start = max(0, np.searchsorted(input_time, zoom_start) - 1)
                input_stop = min(len(input_time), np.searchsorted(input_time, zoom_stop, side="right") + 1)
                output_start = max(0, np.searchsorted(time, zoom_start) - 1)
                output_stop = min(len(time), np.searchsorted(time, zoom_stop, side="right") + 1)
                detail.plot(input_time[input_start:input_stop], before[input_start:input_stop],
                            "k.", markersize=2)
                detail.plot(time[output_start:output_stop], after[output_start:output_stop],
                            "r.-", markersize=2, linewidth=0.7)
                detail.set_xlim(zoom_start, zoom_stop)
                detail.set_title(f"Detail: {zoom_start:.2f}–{zoom_stop:.2f} s", fontsize=7)
                detail.tick_params(labelsize=6, pad=1)
                detail.locator_params(axis="both", nbins=3)
                detail.grid(alpha=0.2)
                axis.axvspan(zoom_start, zoom_stop, color="tab:blue", alpha=0.10)

        plt.show(block=False)

    return np.column_stack((time, cleaned))



def main() -> None:
    """Select input folders, optionally save a CSV, and display cleaned data."""
    root = tk.Tk()
    root.title("Eye Tracking · Data Cleaning")
    root.geometry("780x480")
    root.minsize(680, 440)

    style = ttk.Style(root)
    style.configure("Title.TLabel", font=("Helvetica", 22, "bold"))
    style.configure("Subtitle.TLabel", font=("Helvetica", 11))
    style.configure("Section.TLabel", font=("Helvetica", 11, "bold"))
    style.configure("Folder.TButton", padding=(12, 18))
    style.configure("Run.TButton", padding=(18, 8))

    raw_path = tk.StringVar(root, value="")
    mapped_path = tk.StringVar(root, value="")
    output_path = tk.StringVar(root, value="")
    status = tk.StringVar(root, value="Choose both input folders to continue.")
    run_requested: bool = False

    def select_folder(path: tk.StringVar, title: str) -> None:
        selected: str = filedialog.askdirectory(
            parent=root, title=title, mustexist=True,
            initialdir=path.get() or str(pathlib.Path.home()),
        )
        if selected:
            path.set(selected)
        ready: bool = bool(raw_path.get() and mapped_path.get())
        run_button.configure(state="normal" if ready else "disabled")
        status.set("Ready to clean and plot." if ready else "Choose both input folders to continue.")

    def select_output() -> None:
        current = pathlib.Path(output_path.get()) if output_path.get() else None
        selected: str = filedialog.asksaveasfilename(
            parent=root, title="Save cleaned eye data",
            initialdir=str(current.parent) if current else mapped_path.get() or str(pathlib.Path.home()),
            initialfile=current.name if current else "cleaned_eye_data.csv",
            defaultextension=".csv", filetypes=[("CSV files", "*.csv")],
        )
        if selected:
            output_path.set(selected)

    def run() -> None:
        nonlocal run_requested
        run_requested = True
        root.quit()

    content = ttk.Frame(root, padding=24)
    content.pack(fill="both", expand=True)
    content.columnconfigure((0, 1), weight=1, uniform="folders")
    content.rowconfigure(3, weight=1)

    ttk.Label(content, text="Eye tracking", style="Title.TLabel").grid(
        row=0, column=0, columnspan=2, sticky="w")
    ttk.Label(content, text="Select your recordings to clean and compare the signals.",
              style="Subtitle.TLabel").grid(row=1, column=0, columnspan=2, sticky="w", pady=(6, 20))

    for column, (label, path, description) in enumerate((
        ("RAW", raw_path, "Eye states and blink exports"),
        ("MAPPED", mapped_path, "Alternative camera gaze export"),
    )):
        frame = ttk.Frame(content, padding=14, relief="groove", borderwidth=1)
        frame.grid(row=2, column=column, rowspan=2, sticky="nsew",
                   padx=(0, 8) if column == 0 else (8, 0))
        ttk.Label(frame, text=label, style="Section.TLabel").pack(anchor="w")
        ttk.Label(frame, text=description).pack(anchor="w", pady=(4, 12))
        box = ttk.Button(
            frame, text="Choose folder…", style="Folder.TButton",
            command=lambda value=path, name=label: select_folder(value, f"Select {name} folder"),
        )
        box.pack(fill="x")
        ttk.Label(frame, textvariable=path, wraplength=290).pack(fill="x", pady=(10, 0))

    output_frame = ttk.Frame(content)
    output_frame.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(20, 0))
    output_frame.columnconfigure(0, weight=1)
    ttk.Label(output_frame, text="CSV OUTPUT · Optional", style="Section.TLabel").grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
    output_entry = ttk.Entry(output_frame, textvariable=output_path, state="readonly")
    output_entry.grid(row=1, column=0, sticky="ew")
    output_entry.bind("<Button-1>", lambda event: select_output())
    ttk.Button(output_frame, text="Save as…", command=select_output).grid(row=1, column=1, padx=(8, 0))
    ttk.Button(output_frame, text="Clear", command=lambda: output_path.set("")).grid(row=1, column=2, padx=(4, 0))
    ttk.Label(output_frame, text="Leave blank to plot without saving a file.").grid(
        row=2, column=0, columnspan=3, sticky="w", pady=(6, 0))

    ttk.Separator(content).grid(row=5, column=0, columnspan=2, sticky="ew", pady=(18, 12))
    ttk.Label(content, textvariable=status).grid(row=6, column=0, sticky="w")
    run_button = ttk.Button(content, text="Clean and plot", style="Run.TButton", command=run, state="disabled")
    run_button.grid(row=6, column=1, sticky="e")
    root.protocol("WM_DELETE_WINDOW", root.quit)

    try:
        root.mainloop()
        if not run_requested:
            return
        root.withdraw()
        cleaned_data: np.ndarray = clean_eye_data(raw_path.get(), mapped_path.get(), verbose=True)
        if output_path.get():
            cleaned_df: pd.DataFrame = pd.DataFrame(cleaned_data, columns=[
                "time [s]", "pupil diameter [mm]", "eyelid aperture [mm]",
                "gaze x [px]", "gaze y [px]",
            ])
            cleaned_df.to_csv(output_path.get(), index=False)
    except Exception as error:
        messagebox.showerror("Eye tracking", str(error), parent=root)
        return
    finally:
        root.destroy()

    plt.show(block=True)


if __name__ == "__main__":
    main()
