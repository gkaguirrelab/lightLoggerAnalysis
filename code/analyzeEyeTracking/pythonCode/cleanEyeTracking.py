"""Geoff's eye-cleaning pipeline, with MATLAB reference plots and folder UI."""

import pathlib
import warnings
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import maximum_filter1d

# These defaults reproduce the executable MATLAB source, including values whose
# older comments still refer to a lower sampling rate. No downsampling is done.
BLINK_PRE_NS = 100_000_000
BLINK_POST_NS = 150_000_000
GAZE_FRAME_WIDTH = 640
GAZE_FRAME_HEIGHT = 480
GAZE_SPIKE_THRESHOLD = 40.0
DILATION_SAMPLES = 5
GENERAL_MAX_GAP = 25
BLINK_MAX_GAP = 200
SMOOTHING_SAMPLES = 3
EYE_COLUMNS = (
    "pupil diameter left [mm]", "pupil diameter right [mm]",
    "eyelid aperture left [mm]", "eyelid aperture right [mm]",
)
GAZE_COLUMNS = ("gaze x [px]", "gaze y [px]")
OUTPUT_COLUMNS = (
    "time [s]", "pupil diameter [mm]", "eyelid aperture [mm]",
    "eyelid aperture difference [mm]", "gaze x [px]", "gaze y [px]",
)


def _mean_available(values: np.ndarray) -> np.ndarray:
    """MATLAB mean(..., 2, 'omitnan'), retaining NaN for an all-missing row."""
    counts = np.sum(~np.isnan(values), axis=1)
    return np.divide(np.nansum(values, axis=1), counts,
                     out=np.full(len(values), np.nan), where=counts > 0)


def _fill_missing(values: np.ndarray, method: str, max_gap: int) -> np.ndarray:
    """Match fillmissing with default sample indices and default end filling."""
    result = values.copy()
    valid = np.flatnonzero(~np.isnan(values))
    if valid.size == 0 or (method == "pchip" and valid.size < 2):
        return result
    missing = np.isnan(values)
    transitions = np.diff(np.r_[False, missing, False].astype(np.int8))
    starts = np.flatnonzero(transitions == 1)
    stops = np.flatnonzero(transitions == -1)
    lengths = stops - starts
    # Interior gap size is the distance between bounding valid samples;
    # an edge gap is measured from the vector endpoint to its valid neighbor.
    gap_sizes = lengths + ((starts > 0) & (stops < len(values)))
    indices = np.flatnonzero(missing)[np.repeat(gap_sizes <= max_gap, lengths)]
    if not indices.size:
        return result
    if method == "pchip":
        result[indices] = PchipInterpolator(valid, values[valid], extrapolate=True)(indices)
    else:
        # A tie goes to the later sample, as in MATLAB nearest interpolation.
        right = np.minimum(np.searchsorted(valid, indices), len(valid) - 1)
        left = np.maximum(right - 1, 0)
        nearest = np.where(indices - valid[left] < valid[right] - indices, left, right)
        result[indices] = values[valid[nearest]]
    return result


def _smooth_signals(values: np.ndarray) -> np.ndarray:
    """Three-sample median for physiology and quadratic SG for gaze."""
    result = np.full_like(values, np.nan)
    half = SMOOTHING_SAMPLES // 2
    offsets = np.arange(-half, half + 1, dtype=float)
    bits = 1 << np.arange(SMOOTHING_SAMPLES)
    for column in range(values.shape[1]):
        samples = values[:, column]
        windows = np.lib.stride_tricks.sliding_window_view(
            np.pad(samples, (half, half), constant_values=np.nan), SMOOTHING_SAMPLES)
        if column < 3:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="All-NaN slice encountered", category=RuntimeWarning)
                result[:, column] = np.nanmedian(windows, axis=1)
            continue
        patterns = (~np.isnan(windows)).astype(np.uint8) @ bits
        for pattern in np.unique(patterns):
            if pattern == 0:
                continue
            present = (pattern & bits) != 0
            degree = min(2, int(present.sum()) - 1)
            design = np.vander(offsets[present], degree + 1, increasing=True)
            rows = np.flatnonzero(patterns == pattern)
            result[rows, column] = windows[rows][:, present] @ np.linalg.pinv(design)[0]
        # MATLAB uses a full one-sided window at the recording endpoints.
        width = min(SMOOTHING_SAMPLES, len(samples))
        for index in sorted(set(range(half)) | set(range(max(0, len(samples) - half), len(samples)))):
            start = 0 if index < half else len(samples) - width
            edge = samples[start:start + width]
            present = ~np.isnan(edge)
            if present.any():
                degree = min(2, int(present.sum()) - 1)
                design = np.vander((np.arange(start, start + width) - index)[present], degree + 1, increasing=True)
                result[index, column] = np.linalg.pinv(design)[0] @ edge[present]
    return result


def clean_eye_data(
    raw_neon_path: str,
    processing_neon_path: str,
    verbose: bool = True,
    pupil_diameter_thresholds: tuple[float, float] = (1.5, 8.5),
    pupil_change_threshold: float = 0.25,
    aperture_change_threshold: float = 0.3,
    *,
    plot_label: str | None = None,
    plot_id: str | None = None,
    plot_output_path: str | pathlib.Path | None = None,
    diagnostics: dict | None = None,
) -> np.ndarray:
    """Reproduce processEyeData.m on the overlapping native gaze timeline.

    Read 3d_eye_states.csv and blinks.csv from raw_neon_path, and
    alternative_camera_gaze.csv from processing_neon_path. Return (N, 6):
    elapsed seconds, pupil diameter, mean eyelid aperture, left-minus-right
    eyelid aperture, gaze x, gaze y. Millimeters and pixels are retained.
    Time zero is the first retained gaze sample, NOT the first eye sample.

    The defaults follow the updated processEyeData (1).m supplied on
    October 26, 2026, saved unchanged as matlabCode/processEyeData.m.
    Gaze is cropped to sensor overlap; eye states are interpolated onto it.
    Physiology uses blink/physiology-spike masking, leading/trailing-value carry,
    PCHIP fills and median smoothing. Gaze uses sensor bounds, blink-AND-
    gaze-spike masking, nearest fills and Savitzky–Golay smoothing.
    Apertures remain separate until the final mean and difference.

    Threshold arguments optionally override MATLAB's pupil bounds and
    physiological spike thresholds. verbose shows the MATLAB reference layout:
    five raw/cleaned signal panels followed by unpadded blink events. Its raw
    vectors are captured after alignment and bounds filtering, before masking,
    exactly as in processEyeData.m. plot_label and plot_id identify the recording.
    When verbose is true (the default), save eye_cleaned_vs_raw.pdf in
    processing_neon_path, or at plot_output_path when supplied, then display it.
    The PDF preserves vector graphics. Missing output values remain NaN.
    If supplied, diagnostics is populated with gaze arrays and masks from each
    cleaning step for auditing; it does not change the returned data.

    Raises ValueError for invalid/unsorted timelines, no overlapping samples,
    or reversed blink intervals.
    """
    eye_data = pd.read_csv(pathlib.Path(raw_neon_path) / "3d_eye_states.csv",
                           usecols=["timestamp [ns]", *EYE_COLUMNS])
    gaze_data = pd.read_csv(pathlib.Path(processing_neon_path) / "alternative_camera_gaze.csv",
                            usecols=["timestamp [ns]", *GAZE_COLUMNS])
    blink_data = pd.read_csv(pathlib.Path(raw_neon_path) / "blinks.csv",
                             usecols=["start timestamp [ns]", "end timestamp [ns]"])
    # MATLAB readtable uses double timestamps. Convert before subtraction,
    # interpolation and blink comparisons to reproduce epoch-time rounding.
    eye_timestamps = eye_data["timestamp [ns]"].to_numpy(dtype=float)
    gaze_timestamps = gaze_data["timestamp [ns]"].to_numpy(dtype=float)
    if len(eye_timestamps) < 2 or not len(gaze_timestamps):
        raise ValueError("Need at least two eye samples and one gaze sample.")
    for timestamps in (eye_timestamps, gaze_timestamps):
        if not np.all(np.isfinite(timestamps)) or np.any(np.diff(timestamps) <= 0):
            raise ValueError("Eye and gaze timestamps must be finite and strictly increasing.")
    overlap = ((gaze_timestamps >= max(gaze_timestamps[0], eye_timestamps[0]))
               & (gaze_timestamps <= min(gaze_timestamps[-1], eye_timestamps[-1])))
    timestamps = gaze_timestamps[overlap]
    if not len(timestamps):
        raise ValueError("Eye and gaze timelines have no overlapping gaze samples.")
    time = (timestamps - timestamps[0]) / 1e9
    raw_eyes = eye_data[list(EYE_COLUMNS)].to_numpy(dtype=float)
    eyes = np.column_stack([np.interp(timestamps, eye_timestamps, column) for column in raw_eyes.T])
    raw_gaze = gaze_data[list(GAZE_COLUMNS)].to_numpy(dtype=float)[overlap]
    gaze = raw_gaze.copy()
    out_of_bounds = ((gaze[:, 0] < 0) | (gaze[:, 0] > GAZE_FRAME_WIDTH)
                     | (gaze[:, 1] < 0) | (gaze[:, 1] > GAZE_FRAME_HEIGHT))
    gaze[out_of_bounds] = np.nan

    blink_starts = blink_data["start timestamp [ns]"].to_numpy(dtype=float)
    blink_ends = blink_data["end timestamp [ns]"].to_numpy(dtype=float)
    if np.any(blink_ends < blink_starts):
        raise ValueError("Blink end timestamps must not precede their starts.")
    padded_changes = np.zeros(len(timestamps) + 1, dtype=np.int64)
    np.add.at(padded_changes, np.searchsorted(timestamps, blink_starts - BLINK_PRE_NS, side="left"), 1)
    np.add.at(padded_changes, np.searchsorted(timestamps, blink_ends + BLINK_POST_NS, side="right"), -1)
    padded_blinks = np.cumsum(padded_changes[:-1]) > 0

    pupils = eyes[:, :2]
    pupils[(pupils < pupil_diameter_thresholds[0]) | (pupils > pupil_diameter_thresholds[1])] = np.nan
    pupil = _mean_available(pupils)
    aperture_temp = _mean_available(eyes[:, 2:])
    # The working columns retain each aperture separately through smoothing.
    signals = np.column_stack((pupil, eyes[:, 2:], gaze))
    before_masking = signals.copy() if verbose else None
    physio_spike = ((np.r_[0., np.abs(np.diff(pupil))] > pupil_change_threshold)
                    | (np.r_[0., np.abs(np.diff(aperture_temp))] > aperture_change_threshold))
    gaze_changes = np.vstack((np.zeros(2), np.abs(np.diff(gaze, axis=0))))
    gaze_spike = np.sqrt(np.sum(gaze_changes ** 2, axis=1)) > GAZE_SPIKE_THRESHOLD
    blink_dilated = maximum_filter1d(padded_blinks, DILATION_SAMPLES, mode="constant", cval=0)
    physio_spike_dilated = maximum_filter1d(physio_spike, DILATION_SAMPLES, mode="constant", cval=0)
    gaze_spike_dilated = maximum_filter1d(gaze_spike, DILATION_SAMPLES, mode="constant", cval=0)
    physio_mask = blink_dilated | physio_spike_dilated
    gaze_mask = blink_dilated & gaze_spike_dilated
    signals[physio_mask, :3] = np.nan
    signals[gaze_mask, 3:] = np.nan
    # Match Geoff's endpoint protection before PCHIP: carry the last valid
    # physiology value forward, then the first valid value backward. Both
    # carries are uncapped; all-missing columns remain missing. This prevents
    # extrapolation at either recording edge without changing interior gaps.
    for column in range(3):
        valid = np.flatnonzero(~np.isnan(signals[:, column]))
        if valid.size:
            signals[valid[-1] + 1:, column] = signals[valid[-1], column]
    for column in range(3):
        valid = np.flatnonzero(~np.isnan(signals[:, column]))
        if valid.size:
            signals[:valid[0], column] = signals[valid[0], column]
    imputed = signals.copy()
    for column in range(5):
        method = "pchip" if column < 3 else "nearest"
        longer_fill = _fill_missing(signals[:, column], method, BLINK_MAX_GAP)
        imputed[:, column] = _fill_missing(signals[:, column], method, GENERAL_MAX_GAP)
        overwrite = blink_dilated if column < 3 else gaze_mask
        imputed[overwrite, column] = longer_fill[overwrite]
    smoothed = _smooth_signals(imputed)
    aperture_mean = _mean_available(smoothed[:, 1:3])
    aperture_difference = smoothed[:, 1] - smoothed[:, 2]
    cleaned = np.column_stack((smoothed[:, 0], aperture_mean, aperture_difference, smoothed[:, 3:]))

    if diagnostics is not None:
        diagnostics.update(
            time=time.copy(),
            mapped_gaze=raw_gaze.copy(),
            bounded_gaze=gaze.copy(),
            masked_gaze=signals[:, 3:].copy(),
            filled_gaze=imputed[:, 3:].copy(),
            cleaned_gaze=cleaned[:, 3:].copy(),
            out_of_bounds=out_of_bounds.copy(),
            blink_mask=blink_dilated.copy(),
            gaze_spike_mask=gaze_spike_dilated.copy(),
            gaze_mask=gaze_mask.copy(),
        )

    if verbose:
        # Match MATLAB's pupil_raw/aperture_*_raw/gaze_*_raw snapshots. These
        # are the aligned, bounded signals before blink/spike censoring.
        plot_raw = np.column_stack((
            before_masking[:, 0], _mean_available(before_masking[:, 1:3]),
            before_masking[:, 1] - before_masking[:, 2], before_masking[:, 3:],
        ))
        # MATLAB's last panel shows original blink intervals, without the
        # asymmetric padding or dilation used by the cleaning algorithm.
        raw_blink_changes = np.zeros(len(timestamps) + 1, dtype=np.int64)
        np.add.at(raw_blink_changes, np.searchsorted(timestamps, blink_starts, side="left"), 1)
        np.add.at(raw_blink_changes, np.searchsorted(timestamps, blink_ends, side="right"), -1)
        raw_blinks = np.cumsum(raw_blink_changes[:-1]) > 0
        _plot_matlab_summary(time, plot_raw, cleaned, raw_blinks,
                             plot_label or pathlib.Path(raw_neon_path).name, plot_id,
                             pathlib.Path(plot_output_path) if plot_output_path is not None
                             else pathlib.Path(processing_neon_path) / "eye_cleaned_vs_raw.pdf")
    return np.column_stack((time, cleaned))



def _plot_matlab_summary(time, raw, cleaned, raw_blinks, recording_label, recording_id,
                         output_path):
    """Reproduce the six-panel raw/cleaned figure in the supplied MATLAB/PDF.

    Use the same signal definitions, panel order, titles, units and line styles
    as MATLAB section 9. The wide layout follows the supplied PDF rather than
    the portrait window dimensions in the script. No values are clipped to
    expected physiological bounds: extrapolation artifacts remain visible.
    """
    titles = (
        "Pupil Diameter", "Mean Eyelid Aperture (L/R Combined)",
        "Eyelid Aperture Difference (Left - Right)", "Gaze X", "Gaze Y",
        "Raw Blink Events",
    )
    units = ("Diameter [mm]", "Aperture [mm]", "Difference [mm]",
             "Position [px]", "Position [px]", "Blink Active")
    figure, axes = plt.subplots(6, 1, figsize=(18, 11), sharex=True,
                                constrained_layout=True)
    figure_id = f"{recording_id}-F1" if recording_id else "F1"
    figure.suptitle(f"{figure_id} | {recording_label}", fontsize=12, fontweight="bold")
    for column, axis in enumerate(axes[:5]):
        axis.plot(time, raw[:, column], "k.", markersize=3, label="Raw")
        axis.plot(time, cleaned[:, column], "r-", linewidth=1, label="Cleaned")
    axes[0].legend(loc="upper right", framealpha=1, fancybox=False)
    axes[5].plot(time, raw_blinks.astype(float), "k-", linewidth=1.5)
    axes[5].set_ylim(-0.1, 1.1)
    axes[5].set_xlabel("Time [s]")
    for axis, title, unit in zip(axes, titles, units):
        axis.set_title(title, fontweight="bold")
        axis.set_ylabel(unit)
        axis.grid(alpha=0.3)
        axis.tick_params(axis="x", labelbottom=True, direction="in")
        axis.tick_params(axis="y", direction="in")
        axis.spines[["top", "right"]].set_visible(False)
        axis.margins(x=0)
    # Keep the notebook's existing figure-caption interface.
    figure.eye_cleaning_title = f"{figure_id} | {recording_label}"
    figure.eye_cleaning_figure_id = figure_id
    figure.eye_cleaning_review = []
    figure.set_label("Eye Data Processing: Raw vs Cleaned")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, format="pdf")
    plt.show(block=False)
    return figure


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
        cleaned_data: np.ndarray = clean_eye_data(
            raw_path.get(), mapped_path.get(), verbose=True,
            plot_output_path=(pathlib.Path(output_path.get()).with_name("eye_cleaned_vs_raw.pdf")
                              if output_path.get() else None),
        )
        if output_path.get():
            cleaned_df: pd.DataFrame = pd.DataFrame(cleaned_data, columns=OUTPUT_COLUMNS)
            cleaned_df.to_csv(output_path.get(), index=False)
    except Exception as error:
        messagebox.showerror("Eye tracking", str(error), parent=root)
        return
    finally:
        root.destroy()

    plt.show(block=True)


if __name__ == "__main__":
    main()
