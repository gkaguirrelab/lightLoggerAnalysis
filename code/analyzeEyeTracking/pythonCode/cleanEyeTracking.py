"""Geoff's eye-cleaning pipeline, with MATLAB reference plots and folder UI."""

import pathlib
import warnings
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import maximum_filter1d



def _mean_available(values: np.ndarray) -> np.ndarray:
    """Average each row while ignoring missing values.

    Args:
        values (np.ndarray):
            Two-dimensional signal array.

    Returns:
        One mean per row, with NaN for rows containing only NaNs.

    Notes:
        Matches MATLAB mean(..., 2, 'omitnan').
    """
    # Count available channels in each row so all-missing rows can stay NaN.
    counts: np.ndarray = np.sum(~np.isnan(values), axis=1)
    return np.divide(
        np.nansum(values, axis=1),
        counts,
        out=np.full(len(values), np.nan),
        where=counts > 0,
    )


def _fill_missing(values: np.ndarray, method: str, max_gap: int) -> np.ndarray:
    """Fill eligible missing runs using sample-index interpolation.

    Args:
        values (np.ndarray):
            One-dimensional signal containing NaNs for missing samples.

        method (str):
            Interpolation method: "pchip" or "nearest".

        max_gap (int):
            Maximum gap measured in sample-index distance.

    Returns:
        A copy of the signal with eligible gaps filled.

    Notes:
        Interior gaps include the distance between bounding valid samples.
        Edge gaps use the distance to the nearest valid sample. PCHIP
        extrapolates at endpoints; nearest ties select the later sample.
    """
    # Allocate a separate output so the supplied signal is preserved.
    result: np.ndarray = values.copy()
    # Locate observed samples that can anchor interpolation or endpoint filling.
    valid: np.ndarray = np.flatnonzero(~np.isnan(values))
    # Leave gaps untouched when there are too few observations to interpolate.
    if valid.size == 0 or (method == "pchip" and valid.size < 2):
        return result
    # Mark the samples that need interpolation.
    missing: np.ndarray = np.isnan(values)
    # Pad the missing mask to detect runs touching either endpoint.
    transitions: np.ndarray = np.diff(np.r_[False, missing, False].astype(np.int8))
    # Find the first missing sample of each run.
    starts: np.ndarray = np.flatnonzero(transitions == 1)
    # Find the exclusive end of each missing run.
    stops: np.ndarray = np.flatnonzero(transitions == -1)
    # Measure how many samples each missing run contains.
    lengths: np.ndarray = stops - starts
    # Interior gap size is the distance between bounding valid samples;
    # an edge gap is measured from the vector endpoint to its valid neighbor.
    gap_sizes: np.ndarray = lengths + ((starts > 0) & (stops < len(values)))
    # Select only missing samples whose run meets the maximum gap limit.
    indices: np.ndarray = np.flatnonzero(missing)[
        np.repeat(gap_sizes <= max_gap, lengths)
    ]
    # Return immediately when every gap exceeds the allowed size.
    if not indices.size:
        return result
    # Select the requested interpolation rule for eligible gaps.
    if method == "pchip":
        # Evaluate the shape-preserving interpolant at the eligible sample indices.
        result[indices] = PchipInterpolator(valid, values[valid], extrapolate=True)(
            indices
        )
    else:
        # A tie goes to the later sample, as in MATLAB nearest interpolation.
        right: np.ndarray = np.minimum(np.searchsorted(valid, indices), len(valid) - 1)
        # Find the candidate immediately before the right-hand neighbor.
        left: np.ndarray = np.maximum(right - 1, 0)
        # Choose the closer neighbor; equal distances select the later sample.
        nearest: np.ndarray = np.where(
            indices - valid[left] < valid[right] - indices, left, right
        )
        # Copy the chosen observed values into the eligible missing positions.
        result[indices] = values[valid[nearest]]
    return result


def _smooth_signals(values: np.ndarray, smoothing_samples: int) -> np.ndarray:
    """Apply median smoothing to physiology and quadratic SG smoothing to gaze.

    Args:
        values (np.ndarray):
            Floating-point (N, 5) array: pupil, left/right aperture, gaze x/y.

        smoothing_samples (int):
            Positive odd window size in samples.

    Returns:
        Smoothed array with the same shape; all-missing windows remain NaN.

    Notes:
        Missing samples are omitted. Gaze fits reduce their polynomial degree
        when fewer than three valid samples are available. Endpoint gaze fits
        use one-sided windows, limited to the recording length.
    """
    # Allocate a separate output so the supplied signal is preserved.
    result: np.ndarray = np.full_like(values, np.nan)
    # Use half the odd window width on each side of an interior sample.
    half: int = smoothing_samples // 2
    # Express window positions relative to the sample being estimated.
    offsets: np.ndarray = np.arange(-half, half + 1, dtype=float)
    for column in range(values.shape[1]):
        # Process one physiology or gaze channel at a time.
        samples: np.ndarray = values[:, column]
        # Pad the recording with NaNs to form centered windows at every sample.
        windows: np.ndarray = np.lib.stride_tricks.sliding_window_view(
            np.pad(samples, (half, half), constant_values=np.nan), smoothing_samples
        )
        # The first three channels are physiology and use a NaN-aware median.
        if column < 3:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="All-NaN slice encountered",
                    category=RuntimeWarning,
                )
                result[:, column] = np.nanmedian(windows, axis=1)
            continue
        patterns: np.ndarray
        groups: np.ndarray
        # Group windows sharing the same valid positions to reuse polynomial fits.
        patterns, groups = np.unique(~np.isnan(windows), axis=0, return_inverse=True)
        for group, present in enumerate(patterns):
            if not present.any():
                continue
            # Reduce the quadratic fit order when fewer valid samples are available.
            degree: int = min(2, int(present.sum()) - 1)
            # Build polynomial terms at the available sample positions.
            design: np.ndarray = np.vander(
                offsets[present], degree + 1, increasing=True
            )
            # Select the windows using this pattern of available samples.
            rows: np.ndarray = np.flatnonzero(groups == group)
            # Evaluate each fitted polynomial at zero offset, the window center.
            result[rows, column] = windows[rows][:, present] @ np.linalg.pinv(design)[0]
        # MATLAB uses a full one-sided window at the recording endpoints.
        width: int = min(smoothing_samples, len(samples))
        for index in sorted(
            set(range(min(half, len(samples))))
            | set(range(max(0, len(samples) - half), len(samples)))
        ):
            # Anchor an endpoint window to the beginning or end of the recording.
            start: int = 0 if index < half else len(samples) - width
            # Extract the one-sided window used for this endpoint estimate.
            edge: np.ndarray = samples[start : start + width]
            # Exclude missing samples from the local polynomial fit.
            present: np.ndarray = ~np.isnan(edge)
            if present.any():
                # Reduce the quadratic fit order when fewer valid samples are available.
                degree = min(2, int(present.sum()) - 1)
                # Build polynomial terms at the available sample positions.
                design = np.vander(
                    (np.arange(start, start + width) - index)[present],
                    degree + 1,
                    increasing=True,
                )
                # Evaluate the one-sided fit at the current endpoint sample.
                result[index, column] = np.linalg.pinv(design)[0] @ edge[present]
    return result


def clean_eye_data(
    raw_neon_path: str | pathlib.Path,
    processing_neon_path: str | pathlib.Path,
    verbose: bool = True,
    pupil_diameter_thresholds: tuple[float, float] = (1.5, 8.5),
    pupil_change_threshold: float = 0.25,
    aperture_change_threshold: float = 0.3,
    *,
    blink_pre_ns: int = 100_000_000,
    blink_post_ns: int = 150_000_000,
    gaze_x_upper_bound: int = 640,
    gaze_y_upper_bound: int = 480,
    gaze_spike_threshold: float = 40.0,
    dilation_samples: int = 5,
    general_max_gap: int = 25,
    blink_max_gap: int = 200,
    smoothing_samples: int = 3,
    plot_label: str | None = None,
    plot_id: str | None = None,
    plot_output_path: str | pathlib.Path | None = None,
    diagnostics: dict[str, np.ndarray] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Clean raw Neon eye measurements and mapped gaze into aligned signals.

    Align eye measurements with the overlapping gaze timeline, then identify
    unreliable samples using physiological limits, gaze bounds, blinks, and
    abrupt changes. Mask those samples, fill gaps within the configured limits,
    and smooth the resulting signals. Combine the cleaned eye channels into
    pupil diameter, mean eyelid aperture, and aperture difference, and return
    them alongside gaze, elapsed time, and the cleaning settings. Optionally
    plot a raw-versus-cleaned comparison and collect diagnostic arrays.

    Args:
        raw_neon_path (str | pathlib.Path):
            Folder containing 3d_eye_states.csv and blinks.csv.

        processing_neon_path (str | pathlib.Path):
            Folder containing alternative_camera_gaze.csv.

        verbose (bool):
            Save and display the six-panel raw/cleaned comparison figure.

        pupil_diameter_thresholds (tuple[float, float]):
            Inclusive pupil diameter bounds in millimeters.

        pupil_change_threshold (float):
            Maximum pupil change per sample in millimeters.

        aperture_change_threshold (float):
            Maximum mean aperture change per sample in millimeters.

        blink_pre_ns (int):
            Blink padding before each event, in nanoseconds.

        blink_post_ns (int):
            Blink padding after each event, in nanoseconds.

        gaze_x_upper_bound (int):
            Inclusive upper bound for the horizontal gaze coordinate (x), in pixels.

        gaze_y_upper_bound (int):
            Inclusive upper bound for the vertical gaze coordinate (y), in pixels.

        gaze_spike_threshold (float):
            Maximum Euclidean gaze change per sample, in pixels.

        dilation_samples (int):
            Positive mask dilation window size, in samples.

        general_max_gap (int):
            Maximum general fill gap, in sample-index distance.

        blink_max_gap (int):
            Maximum blink fill gap, in sample-index distance.

        smoothing_samples (int):
            Positive odd median and quadratic SG window size.

        plot_label (str | None):
            Recording label; defaults to the raw folder name.

        plot_id (str | None):
            Optional recording identifier for the figure title.

        plot_output_path (str | pathlib.Path | None):
            PDF destination; defaults to eye_cleaned_vs_raw.pdf in
            processing_neon_path. Used only when verbose is true.

        diagnostics (dict[str, np.ndarray] | None):
            Optional dictionary populated with intermediate gaze arrays and
            masks for auditing.

    Returns:
        A tuple of the (N, 6) cleaned array and a dictionary of the effective
        cleaning parameters. Array columns are elapsed seconds, pupil diameter,
        mean eyelid aperture, left-minus-right eyelid aperture, gaze x, and gaze y.
        Physiology remains in millimeters and gaze in pixels; missing values
        remain NaN. Time zero is the first retained gaze sample.

    Raises:
        ValueError: Invalid window sizes, negative padding or gaps, invalid or
            unsorted timelines, no overlapping samples, or reversed blink events.

    Notes:
        Defaults reproduce processEyeData.m without downsampling. Eye states are
        interpolated onto the overlapping gaze timeline. Physiology uses blink
        and spike masking, endpoint carry, PCHIP fills, and median smoothing.
        Gaze uses bounds filtering, blink-AND-spike masking, nearest fills, and
        quadratic Savitzky–Golay smoothing. Apertures stay separate until the
        final mean and difference. Interior fill gaps include their bounding
        valid-sample distance; edge gaps use the endpoint-to-valid distance.
        Plotting and diagnostics do not change the returned data or metadata.
    """

    # Dilation and smoothing windows count whole samples and must contain at
    # least one sample for maximum_filter1d and sliding_window_view to work.
    for name, value in (
        ("dilation_samples", dilation_samples),
        ("smoothing_samples", smoothing_samples),
    ):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")

    # An odd width puts the estimated sample at the center with equal numbers
    # of neighbors on either side. It also keeps the SG offsets (-half to
    # +half, inclusive) the same length as each smoothing window.
    if smoothing_samples % 2 == 0:
        raise ValueError("smoothing_samples must be odd.")

    # Padding uses whole nanoseconds; fill limits use whole sample-index
    # distances. Negative padding would shrink or reverse blink intervals,
    # and negative gap limits have no meaningful distance interpretation.
    # Zero allows no extra blink padding or no gap filling in that fill pass.
    for name, value in (
        ("blink_pre_ns", blink_pre_ns),
        ("blink_post_ns", blink_post_ns),
        ("general_max_gap", general_max_gap),
        ("blink_max_gap", blink_max_gap),
    ):
        if not isinstance(value, (int, np.integer)) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer.")

    # Read columns in a fixed order before converting to floating-point arrays.
    # Float timestamps match MATLAB's default numeric representation.
    eye_columns: tuple[str, ...] = (
        "timestamp [ns]",
        "pupil diameter left [mm]",
        "pupil diameter right [mm]",
        "eyelid aperture left [mm]",
        "eyelid aperture right [mm]",
    )
    
    eye_data: np.ndarray = pd.read_csv(
        pathlib.Path(raw_neon_path) / "3d_eye_states.csv",
        usecols=list(eye_columns),
    ).loc[:, list(eye_columns)].to_numpy(dtype=float)
    

    # Read gaze already mapped into the world-camera pixel coordinates.
    gaze_columns: tuple[str, ...] = (
        "timestamp [ns]", "gaze x [px]", "gaze y [px]",
    )
    
    gaze_data: np.ndarray = pd.read_csv(
        pathlib.Path(processing_neon_path) / "alternative_camera_gaze.csv",
        usecols=list(gaze_columns),
    ).loc[:, list(gaze_columns)].to_numpy(dtype=float)
    

    # Read blink intervals used to construct the temporal exclusion mask.
    blink_columns: tuple[str, ...] = (
        "start timestamp [ns]", "end timestamp [ns]",
    )
    
    blink_data: np.ndarray = pd.read_csv(
        pathlib.Path(raw_neon_path) / "blinks.csv",
        usecols=list(blink_columns),
    ).loc[:, list(blink_columns)].to_numpy(dtype=float)

    # Copy mapped gaze so bounds filtering does not overwrite diagnostic input.
    bounded_gaze: np.ndarray = gaze_data[:, 1:].copy()

    # Mark rows where either gaze coordinate falls outside the sensor.
    invalid_gaze_pos_mask: np.ndarray = (
        (bounded_gaze[:, 0] < 0)
        | (bounded_gaze[:, 0] > gaze_x_upper_bound)
        | (bounded_gaze[:, 1] < 0)
        | (bounded_gaze[:, 1] > gaze_y_upper_bound)
    )

    # Invalidate both gaze coordinates when either coordinate is outside bounds.
    bounded_gaze[invalid_gaze_pos_mask] = np.nan

    # Extract the sensor timelines after filtering gaze positions.
    eye_timestamps_ns: np.ndarray = eye_data[:, 0]
    gaze_timestamps_ns: np.ndarray = gaze_data[:, 0]
    raw_eye_measurements: np.ndarray = eye_data[:, 1:]

    # Keep gaze samples only where both sensor timelines have coverage.
    # We find that by finding timestamps between whatever sensor started later
    # as well as whatever sensor ended first
    valid_gaze_timestamp_mask: np.ndarray = (
        gaze_timestamps_ns >= max(gaze_timestamps_ns[0], eye_timestamps_ns[0])
    ) & (gaze_timestamps_ns <= min(gaze_timestamps_ns[-1], eye_timestamps_ns[-1]))

    # We will use this period of overlapping gaze and eye states as the timestamps
    valid_gaze_timestamps_ns: np.ndarray = gaze_timestamps_ns[valid_gaze_timestamp_mask]
    if not len(valid_gaze_timestamps_ns):
        raise ValueError("Eye and gaze timelines have no overlapping gaze samples.")

    # Express time in seconds relative to the first retained gaze sample.
    relative_timestamps_s: np.ndarray = (valid_gaze_timestamps_ns - valid_gaze_timestamps_ns[0]) / 1e9

    # Allocate one row per retained gaze timestamp and one column per eye
    # measurement. Exclude the timestamp column from the measurement names.
    eye_channel_names: tuple[str, ...] = eye_columns[1:]
    interpolated_eye_measurements: np.ndarray = np.empty((len(valid_gaze_timestamps_ns), len(eye_channel_names)), dtype=float)

    for idx, column_name in enumerate(eye_channel_names):
        # raw_eye_measurements has the same channel order as eye_channel_names, so idx
        # selects the measurements for column_name (e.g. left pupil diameter).
        channel_values: np.ndarray = raw_eye_measurements[:, idx]

        # Eye and gaze samples may occur at different times. For each target
        # gaze time t, if there is a matching eye timesstamp, 
        # just use that, otherwise, interpolation based on the surrounding eye times t0 and t1 
        interpolated_eye_measurements[:, idx] = np.interp(valid_gaze_timestamps_ns, eye_timestamps_ns, channel_values)

    # Keep the original mapped gaze for diagnostics before bounds filtering.
    overlapping_gaze: np.ndarray = gaze_data[valid_gaze_timestamp_mask, 1:]

    # Retain the same rows in filtered gaze and its position mask so they
    # remain aligned with the output timeline and diagnostic arrays.
    bounded_gaze = bounded_gaze[valid_gaze_timestamp_mask]
    invalid_gaze_pos_mask = invalid_gaze_pos_mask[valid_gaze_timestamp_mask]

    # Extract the start of each unpadded blink event.
    blink_start_timestamps_ns: np.ndarray = blink_data[:, 0]

    # Extract the corresponding end of each blink event.
    blink_end_timestamps_ns: np.ndarray = blink_data[:, 1]
    if np.any(blink_end_timestamps_ns < blink_start_timestamps_ns):
        raise ValueError("Blink end timestamps must not precede their starts.")

    # Allocate interval-boundary counters, including a terminal sentinel.
    padded_changes: np.ndarray = np.zeros(len(valid_gaze_timestamps_ns) + 1, dtype=np.int64)

    # Enter blink intervals at their first included sample.
    np.add.at(
        padded_changes,
        np.searchsorted(valid_gaze_timestamps_ns, blink_start_timestamps_ns - blink_pre_ns, side="left"),
        1,
    )

    # Leave blink intervals immediately after their last included sample.
    np.add.at(
        padded_changes,
        np.searchsorted(valid_gaze_timestamps_ns, blink_end_timestamps_ns + blink_post_ns, side="right"),
        -1,
    )

    # Convert accumulated interval counts into the padded blink mask.
    padded_blinks: np.ndarray = np.cumsum(padded_changes[:-1]) > 0

    # Take a view of the two pupil channels for physiological bounds filtering.
    pupils: np.ndarray = interpolated_eye_measurements[:, :2]

    # Remove pupil values outside the configured physiological range.
    pupils[
        (pupils < pupil_diameter_thresholds[0])
        | (pupils > pupil_diameter_thresholds[1])
    ] = np.nan

    # Combine available left and right pupil diameters at each timestamp.
    pupil: np.ndarray = _mean_available(pupils)

    # Average apertures temporarily to detect abrupt physiological changes.
    aperture_temp: np.ndarray = _mean_available(interpolated_eye_measurements[:, 2:])

    # The working columns retain each aperture separately through smoothing.
    signals: np.ndarray = np.column_stack((pupil, interpolated_eye_measurements[:, 2:], bounded_gaze))

    # Preserve the pre-mask signals only when a comparison plot is requested.
    before_masking: np.ndarray | None = signals.copy() if verbose else None

    # Flag abrupt pupil or mean-aperture changes relative to the prior sample.
    physio_spike: np.ndarray = (
        np.r_[0.0, np.abs(np.diff(pupil))] > pupil_change_threshold
    ) | (np.r_[0.0, np.abs(np.diff(aperture_temp))] > aperture_change_threshold)

    # Compute horizontal and vertical changes; give the first sample zero change.
    gaze_changes: np.ndarray = np.vstack((np.zeros(2), np.abs(np.diff(bounded_gaze, axis=0))))

    # Compare the Euclidean gaze displacement with the spike threshold.
    gaze_spike: np.ndarray = (
        np.sqrt(np.sum(gaze_changes**2, axis=1)) > gaze_spike_threshold
    )

    # Expand the padded blink mask across the configured sample window.
    blink_dilated: np.ndarray = maximum_filter1d(
        padded_blinks, dilation_samples, mode="constant", cval=0
    )

    # Expand physiology spikes to include their nearby samples.
    physio_spike_dilated: np.ndarray = maximum_filter1d(
        physio_spike, dilation_samples, mode="constant", cval=0
    )

    # Expand gaze spikes before intersecting them with blink periods.
    gaze_spike_dilated: np.ndarray = maximum_filter1d(
        gaze_spike, dilation_samples, mode="constant", cval=0
    )

    # Exclude physiology when either a blink or physiological spike is present.
    physio_mask: np.ndarray = blink_dilated | physio_spike_dilated

    # Exclude gaze only where a blink and gaze spike coincide.
    gaze_mask: np.ndarray = blink_dilated & gaze_spike_dilated

    # Mask pupil and both apertures using the physiological exclusion rule.
    signals[physio_mask, :3] = np.nan

    # Apply the stricter blink-and-spike exclusion rule to the two gaze columns.
    signals[gaze_mask, 3:] = np.nan

    # Match Geoff's endpoint protection before PCHIP: carry the last valid
    # physiology value forward, then the first valid value backward. Both
    # carries are uncapped; all-missing columns remain missing. This prevents
    # extrapolation at either recording edge without changing interior gaps.
    for column in range(3):
        # Locate observed samples that can anchor interpolation or endpoint filling.
        valid: np.ndarray = np.flatnonzero(~np.isnan(signals[:, column]))
        if valid.size:
            # Carry the last valid physiology value through the trailing gap.
            signals[valid[-1] + 1 :, column] = signals[valid[-1], column]

    for column in range(3):
        # Locate observed samples that can anchor interpolation or endpoint filling.
        valid = np.flatnonzero(~np.isnan(signals[:, column]))
        if valid.size:
            # Carry the first valid physiology value through the leading gap.
            signals[: valid[0], column] = signals[valid[0], column]

    # Start with masked signals and replace only eligible missing samples.
    imputed: np.ndarray = signals.copy()
    for column in range(5):
        # Use shape-preserving interpolation for physiology and nearest samples for gaze.
        method: str = "pchip" if column < 3 else "nearest"

        # Compute candidate fills using the more permissive blink gap limit.
        longer_fill: np.ndarray = _fill_missing(
            signals[:, column], method, blink_max_gap
        )

        # Apply the general gap limit before allowing blink-specific exceptions.
        imputed[:, column] = _fill_missing(signals[:, column], method, general_max_gap)

        # Select the blink-related positions allowed to use the longer fill limit.
        overwrite: np.ndarray = blink_dilated if column < 3 else gaze_mask

        # Replace only blink-related positions with the longer-gap candidates.
        imputed[overwrite, column] = longer_fill[overwrite]

    # Smooth physiology and gaze with their respective filters.
    smoothed: np.ndarray = _smooth_signals(imputed, smoothing_samples)

    # Combine the independently cleaned left and right apertures.
    aperture_mean: np.ndarray = _mean_available(smoothed[:, 1:3])

    # Retain aperture asymmetry as left minus right.
    aperture_difference: np.ndarray = smoothed[:, 1] - smoothed[:, 2]

    # Assemble the five cleaned output signals in the documented column order.
    cleaned: np.ndarray = np.column_stack(
        (smoothed[:, 0], aperture_mean, aperture_difference, smoothed[:, 3:])
    )

    # Copy intermediate arrays so later inspection cannot mutate working values.
    if diagnostics is not None:
        diagnostics.update(
            time=relative_timestamps_s.copy(),
            mapped_gaze=overlapping_gaze.copy(),
            bounded_gaze=bounded_gaze.copy(),
            masked_gaze=signals[:, 3:].copy(),
            filled_gaze=imputed[:, 3:].copy(),
            cleaned_gaze=cleaned[:, 3:].copy(),
            out_of_bounds=invalid_gaze_pos_mask.copy(),
            blink_mask=blink_dilated.copy(),
            gaze_spike_mask=gaze_spike_dilated.copy(),
            gaze_mask=gaze_mask.copy(),
        )

    if verbose:
        # Match MATLAB's pupil_raw/aperture_*_raw/gaze_*_raw snapshots. These
        # are the aligned, bounded signals before blink/spike censoring.
        plot_raw: np.ndarray = np.column_stack(
            (
                before_masking[:, 0],
                _mean_available(before_masking[:, 1:3]),
                before_masking[:, 1] - before_masking[:, 2],
                before_masking[:, 3:],
            )
        )

        # MATLAB's last panel shows original blink intervals, without the
        # asymmetric padding or dilation used by the cleaning algorithm.
        raw_blink_changes: np.ndarray = np.zeros(len(valid_gaze_timestamps_ns) + 1, dtype=np.int64)

        # Enter blink intervals at their first included sample.
        np.add.at(
            raw_blink_changes, np.searchsorted(valid_gaze_timestamps_ns, blink_start_timestamps_ns, side="left"), 1
        )

        # Leave blink intervals immediately after their last included sample.
        np.add.at(
            raw_blink_changes, np.searchsorted(valid_gaze_timestamps_ns, blink_end_timestamps_ns, side="right"), -1
        )

        # Convert original interval counts into the blink trace shown in the plot.
        raw_blinks: np.ndarray = np.cumsum(raw_blink_changes[:-1]) > 0
        _plot_matlab_summary(
            relative_timestamps_s,
            plot_raw,
            cleaned,
            raw_blinks,
            plot_label or pathlib.Path(raw_neon_path).name,
            plot_id,
            (
                pathlib.Path(plot_output_path)
                if plot_output_path is not None
                else pathlib.Path(processing_neon_path) / "eye_cleaned_vs_raw.pdf"
            ),
        )

    # Record the effective cleaning settings alongside the returned data.
    metadata: dict[str, object] = {
        "pupil_diameter_thresholds": pupil_diameter_thresholds,
        "pupil_change_threshold": pupil_change_threshold,
        "aperture_change_threshold": aperture_change_threshold,
        "blink_pre_ns": blink_pre_ns,
        "blink_post_ns": blink_post_ns,
        "gaze_x_upper_bound": gaze_x_upper_bound,
        "gaze_y_upper_bound": gaze_y_upper_bound,
        "gaze_spike_threshold": gaze_spike_threshold,
        "dilation_samples": dilation_samples,
        "general_max_gap": general_max_gap,
        "blink_max_gap": blink_max_gap,
        "smoothing_samples": smoothing_samples,
    }

    # Prepend elapsed time and return the settings that produced these signals.
    return np.column_stack((relative_timestamps_s, cleaned)), metadata


def _plot_matlab_summary(
    time: np.ndarray,
    raw: np.ndarray,
    cleaned: np.ndarray,
    raw_blinks: np.ndarray,
    recording_label: str,
    recording_id: str | None,
    output_path: pathlib.Path,
) -> Figure:
    """Save and display the six-panel MATLAB-style raw/cleaned figure.

    Args:
        time (np.ndarray):
            Elapsed sample times in seconds.

        raw (np.ndarray):
            (N, 5) aligned signals before masking.

        cleaned (np.ndarray):
            (N, 5) cleaned signals in output-column order, excluding time.

        raw_blinks (np.ndarray):
            Boolean mask of unpadded blink events.

        recording_label (str):
            Human-readable recording label.

        recording_id (str | None):
            Optional identifier prepended to the figure title.

        output_path (pathlib.Path):
            Destination for the vector PDF.

    Returns:
        The saved Matplotlib figure, displayed without blocking.

    Notes:
        Panels show pupil diameter, mean aperture, aperture difference, gaze x,
        gaze y, and raw blinks. Values are not clipped for display. Custom figure
        attributes preserve the existing notebook caption interface.
    """
    # Label each panel with the physical signal it displays.
    titles: tuple[str, ...] = (
        "Pupil Diameter",
        "Mean Eyelid Aperture (L/R Combined)",
        "Eyelid Aperture Difference (Left - Right)",
        "Gaze X",
        "Gaze Y",
        "Raw Blink Events",
    )
    # Associate each signal panel with its measurement units.
    units: tuple[str, ...] = (
        "Diameter [mm]",
        "Aperture [mm]",
        "Difference [mm]",
        "Position [px]",
        "Position [px]",
        "Blink Active",
    )
    figure: Figure
    axes: np.ndarray
    # Create six aligned panels sharing the same time axis.
    figure, axes = plt.subplots(
        6, 1, figsize=(18, 11), sharex=True, constrained_layout=True
    )
    # Build a stable figure identifier for notebook captions.
    figure_id: str = f"{recording_id}-F1" if recording_id else "F1"
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
    # Create the PDF destination folder if this is its first output.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Preserve line and text quality by saving the comparison as a vector PDF.
    figure.savefig(output_path, format="pdf")
    # Display the figure without stopping a notebook or batch run.
    plt.show(block=False)
    return figure


def main() -> None:
    """Select recording folders, clean eye data, and optionally save a CSV.

    Args:
        None.

    Returns:
        None.

    Notes:
        Runs an interactive Tk folder selector and displays Matplotlib figures.
        The optional CSV contains the data array; metadata is returned by
        clean_eye_data for programmatic callers.
    """
    # Create the folder-selection window before constructing its controls.
    root: tk.Tk = tk.Tk()
    root.title("Eye Tracking · Data Cleaning")
    root.geometry("780x480")
    root.minsize(680, 440)

    # Define shared widget styles for the folder-selection interface.
    style: ttk.Style = ttk.Style(root)
    style.configure("Title.TLabel", font=("Helvetica", 22, "bold"))
    style.configure("Subtitle.TLabel", font=("Helvetica", 11))
    style.configure("Section.TLabel", font=("Helvetica", 11, "bold"))
    style.configure("Folder.TButton", padding=(12, 18))
    style.configure("Run.TButton", padding=(18, 8))

    # Store the selected raw-data folder in a Tk-managed string.
    raw_path: tk.StringVar = tk.StringVar(root, value="")
    # Store the selected mapped-gaze folder independently of the raw folder.
    mapped_path: tk.StringVar = tk.StringVar(root, value="")
    # Track the optional CSV destination chosen by the user.
    output_path: tk.StringVar = tk.StringVar(root, value="")
    # Display whether the selected inputs are ready for processing.
    status: tk.StringVar = tk.StringVar(
        root, value="Choose both input folders to continue."
    )
    # Distinguish pressing Run from closing the selection window.
    run_requested: bool = False

    def select_folder(path: tk.StringVar, title: str) -> None:
        """Select an input folder and update the run-button state.

        Args:
            path (tk.StringVar):
                Tk variable receiving the chosen folder.

            title (str):
                Folder-dialog title.

        Returns:
            None.

        Notes:
            Cancellation preserves the current selection.
        """
        # Open the selection dialog and retain its result for cancellation handling.
        selected: str = filedialog.askdirectory(
            parent=root,
            title=title,
            mustexist=True,
            initialdir=path.get() or str(pathlib.Path.home()),
        )
        if selected:
            path.set(selected)
        # Enable cleaning only after both input folders have been selected.
        ready: bool = bool(raw_path.get() and mapped_path.get())
        run_button.configure(state="normal" if ready else "disabled")
        status.set(
            "Ready to clean and plot."
            if ready
            else "Choose both input folders to continue."
        )

    def select_output() -> None:
        """Choose the optional cleaned-data CSV destination.

        Args:
            None.

        Returns:
            None.

        Notes:
            Cancellation preserves the current destination.
        """
        # Use the previous output selection to initialize the save dialog.
        current: pathlib.Path | None = (
            pathlib.Path(output_path.get()) if output_path.get() else None
        )
        # Open the selection dialog and retain its result for cancellation handling.
        selected: str = filedialog.asksaveasfilename(
            parent=root,
            title="Save cleaned eye data",
            initialdir=(
                str(current.parent)
                if current
                else mapped_path.get() or str(pathlib.Path.home())
            ),
            initialfile=current.name if current else "cleaned_eye_data.csv",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
        )
        if selected:
            output_path.set(selected)

    def run() -> None:
        """Request cleaning and exit the folder-selection event loop.

        Args:
            None.

        Returns:
            None.

        Notes:
            Processing starts after the event loop returns.
        """
        nonlocal run_requested
        # Distinguish pressing Run from closing the selection window.
        run_requested = True
        root.quit()

    # Create the container that holds the folder and output controls.
    content: ttk.Frame = ttk.Frame(root, padding=24)
    content.pack(fill="both", expand=True)
    content.columnconfigure((0, 1), weight=1, uniform="folders")
    content.rowconfigure(3, weight=1)

    ttk.Label(content, text="Eye tracking", style="Title.TLabel").grid(
        row=0, column=0, columnspan=2, sticky="w"
    )
    ttk.Label(
        content,
        text="Select your recordings to clean and compare the signals.",
        style="Subtitle.TLabel",
    ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(6, 20))

    for column, (label, path, description) in enumerate(
        (
            ("RAW", raw_path, "Eye states and blink exports"),
            ("MAPPED", mapped_path, "Alternative camera gaze export"),
        )
    ):
        # Group one input-folder button with its label and selected path.
        frame: ttk.Frame = ttk.Frame(
            content, padding=14, relief="groove", borderwidth=1
        )
        frame.grid(
            row=2,
            column=column,
            rowspan=2,
            sticky="nsew",
            padx=(0, 8) if column == 0 else (8, 0),
        )
        ttk.Label(frame, text=label, style="Section.TLabel").pack(anchor="w")
        ttk.Label(frame, text=description).pack(anchor="w", pady=(4, 12))
        # Bind this folder button to its own path variable and dialog title.
        box: ttk.Button = ttk.Button(
            frame,
            text="Choose folder…",
            style="Folder.TButton",
            command=lambda value=path, name=label: select_folder(
                value, f"Select {name} folder"
            ),
        )
        box.pack(fill="x")
        ttk.Label(frame, textvariable=path, wraplength=290).pack(fill="x", pady=(10, 0))

    # Group the optional CSV destination controls below the inputs.
    output_frame: ttk.Frame = ttk.Frame(content)
    output_frame.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(20, 0))
    output_frame.columnconfigure(0, weight=1)
    ttk.Label(output_frame, text="CSV OUTPUT · Optional", style="Section.TLabel").grid(
        row=0, column=0, columnspan=3, sticky="w", pady=(0, 8)
    )
    # Show the destination read-only so selection goes through the save dialog.
    output_entry: ttk.Entry = ttk.Entry(
        output_frame, textvariable=output_path, state="readonly"
    )
    output_entry.grid(row=1, column=0, sticky="ew")
    output_entry.bind("<Button-1>", lambda event: select_output())
    ttk.Button(output_frame, text="Save as…", command=select_output).grid(
        row=1, column=1, padx=(8, 0)
    )
    ttk.Button(output_frame, text="Clear", command=lambda: output_path.set("")).grid(
        row=1, column=2, padx=(4, 0)
    )
    ttk.Label(output_frame, text="Leave blank to plot without saving a file.").grid(
        row=2, column=0, columnspan=3, sticky="w", pady=(6, 0)
    )

    ttk.Separator(content).grid(
        row=5, column=0, columnspan=2, sticky="ew", pady=(18, 12)
    )
    ttk.Label(content, textvariable=status).grid(row=6, column=0, sticky="w")
    # Keep the processing button disabled until both input folders are ready.
    run_button: ttk.Button = ttk.Button(
        content,
        text="Clean and plot",
        style="Run.TButton",
        command=run,
        state="disabled",
    )
    run_button.grid(row=6, column=1, sticky="e")
    root.protocol("WM_DELETE_WINDOW", root.quit)

    try:
        root.mainloop()
        if not run_requested:
            return
        root.withdraw()
        cleaned_data: np.ndarray
        metadata: dict[str, object]
        # Run cleaning and separate the numerical output from its settings.
        cleaned_data, metadata = clean_eye_data(
            raw_path.get(),
            mapped_path.get(),
            verbose=True,
            plot_output_path=(
                pathlib.Path(output_path.get()).with_name("eye_cleaned_vs_raw.pdf")
                if output_path.get()
                else None
            ),
        )
        if output_path.get():
            output_columns: tuple[str, ...] = (
                "time [s]",
                "pupil diameter [mm]",
                "eyelid aperture [mm]",
                "eyelid aperture difference [mm]",
                "gaze x [px]",
                "gaze y [px]",
            )
            # Label the cleaned array columns before writing the optional CSV.
            cleaned_df: pd.DataFrame = pd.DataFrame(
                cleaned_data, columns=output_columns
            )
            cleaned_df.to_csv(output_path.get(), index=False)
    except Exception as error:
        messagebox.showerror("Eye tracking", str(error), parent=root)
        return
    finally:
        root.destroy()

    plt.show(block=True)


if __name__ == "__main__":
    main()
