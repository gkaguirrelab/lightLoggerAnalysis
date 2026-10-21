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
) -> np.ndarray:
    """Return [time_s, pupil_mm, aperture_mm, gaze_x_px, gaze_y_px].

    Read eye states/blinks from raw_neon_path and alternative_camera_gaze.csv
    from processing_neon_path. Match the MATLAB cleaning routine, including
    NaN and edge behavior. Set verbose=True for the six-panel comparison.
    """
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

    #2. Extract the timestamps 
    timestamps: np.ndarray = eye_data["timestamp [ns]"].to_numpy(dtype=np.int64)
    gaze_timestamps: np.ndarray = gaze_data["timestamp [ns]"].to_numpy(dtype=np.int64)

    # Validate the timestamps. Must be more than 1 sample to do diff and they must be monotonic 
    if len(timestamps) < 2 or len(gaze_timestamps) < 2:
        raise ValueError("Eye and gaze exports must each contain at least two samples.")
    if np.any(np.diff(timestamps) <= 0) or np.any(np.diff(gaze_timestamps) <= 0):
        raise ValueError("Eye and gaze timestamps must be strictly increasing and unique.")

    # Subtract integer timestamps before conversion to preserve precision.
    # time is the eye-state/output timeline; gaze_time keeps the gaze sample times.
    # Both use the FIRST EYE timestamp as zero so their recording offset is retained.
    time: np.ndarray = (timestamps - timestamps[0]).astype(float) / 1e9
    gaze_time: np.ndarray = (gaze_timestamps - timestamps[0]).astype(float) / 1e9
    
    # eyes is the working four-column array: left/right pupil, then left/right aperture.
    eyes: np.ndarray = eye_data[eye_columns].to_numpy(dtype=float)

    # Keep a separate copy because pupil validity filtering below modifies eyes.
    # This copy supplies the original black dots in the comparison plots.
    raw_eyes: np.ndarray = eyes.copy()

    # raw_gaze already contains camera-mapped x/y, but still at gaze timestamps.
    # It stays unchanged; mapped_gaze below is the working, time-aligned version.
    raw_gaze: np.ndarray = gaze_data[gaze_columns].to_numpy(dtype=float)

    # 3. Synchronize gaze, excluding large gaps in the source timestamps.
    # Resample gaze onto eye times so one mask can address all four signals.
    # Samples outside gaze coverage remain NaN rather than being extrapolated.
    mapped_gaze: np.ndarray = np.column_stack([
        np.interp(time, gaze_time, raw_gaze[:, column], left=np.nan, right=np.nan)
        for column in range(2)
    ])

    # Interpolation alone would draw a line across long periods with no gaze data.
    max_source_gap_ns: float = 15 * float(np.median(np.diff(timestamps)))
    gap_indices: np.ndarray = np.flatnonzero(np.diff(gaze_timestamps) > max_source_gap_ns)

    # Binary searches find the eye samples inside each gap, retaining its endpoints.
    gap_starts: np.ndarray = np.searchsorted(
        timestamps, gaze_timestamps[gap_indices], side="right"
    )
    gap_ends: np.ndarray = np.searchsorted(
        timestamps, gaze_timestamps[gap_indices + 1], side="left"
    )

    for start, end in zip(gap_starts, gap_ends):
        mapped_gaze[start:end] = np.nan

    # 4. Build asymmetric blink masks (100 ms before, 150 ms after).
    blink_starts: np.ndarray = blink_data["start timestamp [ns]"].to_numpy(dtype=np.int64)
    blink_ends: np.ndarray = blink_data["end timestamp [ns]"].to_numpy(dtype=np.int64)

    if np.any(blink_ends < blink_starts):
        raise ValueError("Blink end timestamps must not precede their starts.")

    # The unpadded mask is for display; padding catches partial eyelid closure.
    blink_masks: list[np.ndarray] = []
    for pre_padding, post_padding in ((0, 0), (100_000_000, 150_000_000)):
        starts: np.ndarray = np.searchsorted(timestamps, blink_starts - pre_padding)
        ends: np.ndarray = np.searchsorted(timestamps, blink_ends + post_padding, side="right")

        # +1 at each start and -1 after each end count how many blinks cover a sample.
        # The cumulative sum merges overlapping intervals without full-timeline scans.
        changes: np.ndarray = np.zeros(len(timestamps) + 1, dtype=np.int64)
        np.add.at(changes, starts, 1)
        np.add.at(changes, ends, -1)
        blink_masks.append(np.cumsum(changes[:-1]) > 0)

    unpadded_blinks, is_blink = blink_masks

    # 5. Apply pupil bounds independently, then average the available eyes.
    # pupils is a view into eyes, so these assignments also update eyes.
    pupils: np.ndarray = eyes[:, :2]
    pupils[(pupils < 1.5) | (pupils > 8.5)] = np.nan

    # Average available eyes; keep NaN if both are missing. Also prepare raw plot data.
    combined: list[np.ndarray] = []
    for bilateral in (eyes[:, :2], eyes[:, 2:], raw_eyes[:, :2], raw_eyes[:, 2:]):
        valid_counts: np.ndarray = np.sum(~np.isnan(bilateral), axis=1)
        combined.append(np.divide(
            np.nansum(bilateral, axis=1), valid_counts,
            out=np.full(len(time), np.nan), where=valid_counts > 0,
        ))

    pupil, aperture, raw_pupil, raw_aperture = combined

    # signals collects pupil, aperture, gaze x, and gaze y on the same timeline.
    signals: np.ndarray = np.column_stack((pupil, aperture, mapped_gaze))

    # 6. Detect spikes, dilate masks independently, and censor all signals.
    # Compute changes before masking so lid transitions can still trigger rejection.
    pupil_changes: np.ndarray = np.r_[0.0, np.abs(np.diff(pupil))]
    aperture_changes: np.ndarray = np.r_[0.0, np.abs(np.diff(aperture))]
    is_spike: np.ndarray = (pupil_changes > 0.25) | (aperture_changes > 0.3)

    # An 11-sample window extends rejection by five samples on either side.
    # Separate masks let blink filling use a different limit from other tracking loss.
    blink_dilated: np.ndarray = maximum_filter1d(is_blink, 11, mode="constant", cval=0)
    spike_dilated: np.ndarray = maximum_filter1d(is_spike, 11, mode="constant", cval=0)

    signals[blink_dilated | spike_dilated] = np.nan

    # 7. Fill general gaps up to 50 samples and blink gaps up to 400 samples.
    # imputed holds gap-filled data; cleaned below holds the final smoothed data.
    imputed: np.ndarray = signals.copy()
    for column in range(4):
        values: np.ndarray = signals[:, column]
        valid: np.ndarray = np.flatnonzero(~np.isnan(values))
        if valid.size == 0:
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

    # Plot raw data on its own timestamps and cleaned data on the common eye timeline.
    if verbose:
        _, axes = plt.subplots(3, 2, sharex=True, figsize=(16, 9), constrained_layout=True)
        names: tuple[str, ...] = ("Pupil (combined)", "Eyelid aperture (combined)", "Gaze x", "Gaze y")
        units: tuple[str, ...] = ("Diameter [mm]", "Aperture [mm]", "Position [px]", "Position [px]")
        raw: tuple[np.ndarray, ...] = (raw_pupil, raw_aperture, raw_gaze[:, 0], raw_gaze[:, 1])

        for column in range(4):
            axis = axes.flat[column]
            axis.plot(time if column < 2 else gaze_time, raw[column], "k.", markersize=2, label="Raw")
            axis.plot(time, cleaned[:, column], "r-", linewidth=0.5, label="Cleaned")
            axis.set(title=names[column], ylabel=units[column], xlabel="Time [s]")
            axis.legend()
            axis.grid(alpha=0.25)

        # Repeat unpadded blinks beneath each column for visual comparison.
        for axis in axes[2]:
            axis.fill_between(time, unpadded_blinks, step="post", color="#0059b3", alpha=0.25)
            axis.step(time, unpadded_blinks, where="post", color="#0059b3", linewidth=1.5)
            axis.set(
                title="Blinks (no padding)", ylabel="Blink state", xlabel="Time [s]",
                ylim=(-0.1, 1.1), yticks=[0, 1], yticklabels=["No blink", "Blink"],
            )
            axis.grid(axis="x", alpha=0.12)

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
