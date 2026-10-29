"""Geoff's eye-cleaning pipeline, with MATLAB reference plots and folder UI."""

import pathlib
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import maximum_filter1d



def _fill_missing(values: np.ndarray, method: str, max_gap: int) -> np.ndarray:
    """Fill eligible missing runs using sample-index interpolation.

    Find consecutive missing samples, select gaps within the allowed size,
    and fill those samples from the original observations. Physiology uses
    a shape-preserving PCHIP curve; gaze uses the nearest observed sample.
    Longer gaps stay missing, and the supplied array is not changed.

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
    filled_values: np.ndarray = values.copy()
    missing_sample_mask: np.ndarray = np.isnan(values)
    valid_sample_indices: np.ndarray = np.where(~missing_sample_mask)[0]

    # Nearest filling needs one observed sample; a PCHIP curve needs two.
    minimum_valid_samples: int = 2 if method == "pchip" else 1
    if len(valid_sample_indices) < minimum_valid_samples:
        return filled_values

    # A gap starts where a missing sample follows an observed sample, and ends
    # where a missing sample precedes an observed sample. Treat positions just
    # outside the recording as observed so edge gaps are detected too.
    previous_sample_missing: np.ndarray = np.concatenate(([False], missing_sample_mask[:-1]))
    next_sample_missing: np.ndarray = np.concatenate((missing_sample_mask[1:], [False]))
    gap_start_indices: np.ndarray = np.where(
        missing_sample_mask & ~previous_sample_missing
    )[0]
    gap_end_indices: np.ndarray = np.where(
        missing_sample_mask & ~next_sample_missing
    )[0]

    # Decide which complete gaps are short enough to fill. For an interior
    # gap, count the distance between its bounding observations: two missing
    # samples have a gap distance of three. At an edge, use the missing count.
    fill_sample_mask: np.ndarray = np.zeros(len(values), dtype=bool)
    for gap_start, gap_end in zip(gap_start_indices, gap_end_indices):
        gap_distance: int = int(gap_end - gap_start + 1)
        has_observation_before: bool = gap_start > 0
        has_observation_after: bool = gap_end < len(values) - 1

        if has_observation_before and has_observation_after:
            gap_distance += 1

        if gap_distance <= max_gap:
            fill_sample_mask[gap_start : gap_end + 1] = True

    fill_sample_indices: np.ndarray = np.where(fill_sample_mask)[0]
    if len(fill_sample_indices) == 0:
        return filled_values

    if method == "pchip":
        # Build one shape-preserving curve from the original observations.
        # Evaluate it only inside eligible gaps. Edge gaps use extrapolation.
        interpolator = PchipInterpolator(
            valid_sample_indices, values[valid_sample_indices], extrapolate=True
        )
        filled_values[fill_sample_indices] = interpolator(fill_sample_indices)
    else:
        # Find the observed samples immediately before and after each missing
        # sample. Clip to the first/last observation for gaps at either edge.
        right_neighbor_positions: np.ndarray = np.searchsorted(
            valid_sample_indices, fill_sample_indices
        )
        left_neighbor_positions: np.ndarray = np.clip(
            right_neighbor_positions - 1, 0, len(valid_sample_indices) - 1
        )
        right_neighbor_positions = np.clip(
            right_neighbor_positions, 0, len(valid_sample_indices) - 1
        )
        left_sample_indices: np.ndarray = valid_sample_indices[left_neighbor_positions]
        right_sample_indices: np.ndarray = valid_sample_indices[right_neighbor_positions]

        # Use the closer observation. A tie uses the later sample, matching
        # MATLAB nearest interpolation. New fills never become anchor samples.
        distance_to_left: np.ndarray = np.abs(fill_sample_indices - left_sample_indices)
        distance_to_right: np.ndarray = np.abs(right_sample_indices - fill_sample_indices)
        use_left_sample: np.ndarray = distance_to_left < distance_to_right
        nearest_sample_indices: np.ndarray = np.where(
            use_left_sample, left_sample_indices, right_sample_indices
        )
        filled_values[fill_sample_indices] = values[nearest_sample_indices]

    return filled_values


def _smooth_signals(values: np.ndarray, smoothing_samples: int) -> np.ndarray:
    """Apply median smoothing to physiology and quadratic SG smoothing to gaze.

    Build a neighborhood around each sample and ignore missing observations.
    For pupil and aperture channels, take the neighborhood median. For gaze,
    fit a local polynomial and estimate its value at the current sample.
    Recalculate gaze endpoints with one-sided windows to match MATLAB.
    Return a separate array, preserving the supplied measurements.

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
    signal_columns: tuple[str, ...] = (
        "pupil diameter [mm]",
        "eyelid aperture left [mm]",
        "eyelid aperture right [mm]",
        "gaze x [px]",
        "gaze y [px]",
    )
    gaze_columns: slice = slice(
        signal_columns.index("gaze x [px]"),
        signal_columns.index("gaze y [px]") + 1,
    )
    smoothed_signals: np.ndarray = np.full_like(values, np.nan)
    half_window_size: int = smoothing_samples // 2

    # For a five-sample window, offsets are [-2, -1, 0, 1, 2]. Zero is the
    # sample being estimated; the other positions are its nearby samples.
    window_offsets: np.ndarray = np.arange(
        -half_window_size, half_window_size + 1, dtype=float
    )

    for column, column_name in enumerate(signal_columns):
        channel_values: np.ndarray = values[:, column]

        # Add missing values outside the recording, then take one centered
        # window per sample. Each row contains that sample's neighborhood.
        padded_channel_values: np.ndarray = np.pad(
            channel_values, (half_window_size, half_window_size), constant_values=np.nan
        )
        sample_windows: np.ndarray = np.lib.stride_tricks.sliding_window_view(
            padded_channel_values, smoothing_samples
        )
        valid_window_samples: np.ndarray = ~np.isnan(sample_windows)

        if column not in range(gaze_columns.start, gaze_columns.stop):
            # Pupil and aperture channels use the median of available nearby
            # values. Skip all-missing windows so they stay NaN without warnings.
            has_valid_window_sample: np.ndarray = np.any(valid_window_samples, axis=1)
            smoothed_signals[has_valid_window_sample, column] = np.nanmedian(
                sample_windows[has_valid_window_sample], axis=1
            )
            continue

        # Gaze uses a local quadratic fit (Savitzky-Golay smoothing). Windows
        # with the same missing positions share the same fit weights, so group
        # them to avoid solving an identical fitting problem for every sample.
        valid_sample_patterns: np.ndarray
        window_pattern_indices: np.ndarray
        valid_sample_patterns, window_pattern_indices = np.unique(
            valid_window_samples, axis=0, return_inverse=True
        )

        for pattern_index, valid_position_mask in enumerate(valid_sample_patterns):
            valid_sample_count: int = int(np.sum(valid_position_mask))
            if valid_sample_count == 0:
                continue

            # Three or more observations support a quadratic; two support a
            # line, and one supports only a constant. Ignore missing samples.
            polynomial_degree: int = min(2, valid_sample_count - 1)
            valid_offsets: np.ndarray = window_offsets[valid_position_mask]
            polynomial_terms: np.ndarray = np.vander(
                valid_offsets, polynomial_degree + 1, increasing=True
            )

            # Each row of polynomial_terms is [1, offset, offset**2] for a
            # quadratic. The pseudoinverse gives least-squares fit weights.
            # At offset zero only the constant coefficient contributes, so
            # its weights (row zero) directly estimate the center value.
            center_value_weights: np.ndarray = np.linalg.pinv(polynomial_terms)[0]
            matching_window_indices: np.ndarray = np.where(
                window_pattern_indices == pattern_index
            )[0]
            matching_window_values: np.ndarray = sample_windows[matching_window_indices]
            observed_window_values: np.ndarray = matching_window_values[:, valid_position_mask]
            smoothed_signals[matching_window_indices, column] = (
                observed_window_values @ center_value_weights
            )

        # MATLAB uses a full window at each gaze endpoint instead of a window
        # shortened by padding. Recalculate those samples using the first or
        # last window of the recording, capped by the recording length.
        sample_count: int = len(channel_values)
        endpoint_window_size: int = min(smoothing_samples, sample_count)
        sample_indices: np.ndarray = np.arange(sample_count)
        endpoint_sample_mask: np.ndarray = (
            (sample_indices < half_window_size)
            | (sample_indices >= sample_count - half_window_size)
        )
        endpoint_sample_indices: np.ndarray = np.where(endpoint_sample_mask)[0]

        for sample_index in endpoint_sample_indices:
            if sample_index < half_window_size:
                window_start: int = 0
            else:
                window_start = sample_count - endpoint_window_size

            window_stop: int = window_start + endpoint_window_size
            endpoint_values: np.ndarray = channel_values[window_start:window_stop]
            valid_endpoint_mask: np.ndarray = ~np.isnan(endpoint_values)
            valid_sample_count = int(np.sum(valid_endpoint_mask))
            if valid_sample_count == 0:
                continue

            # Express positions relative to this endpoint sample so evaluating
            # at zero again gives its estimate, even with a one-sided window.
            endpoint_offsets: np.ndarray = np.arange(window_start, window_stop) - sample_index
            polynomial_degree = min(2, valid_sample_count - 1)
            polynomial_terms = np.vander(
                endpoint_offsets[valid_endpoint_mask], polynomial_degree + 1, increasing=True
            )
            center_value_weights = np.linalg.pinv(polynomial_terms)[0]
            smoothed_signals[sample_index, column] = (
                center_value_weights @ endpoint_values[valid_endpoint_mask]
            )

    return smoothed_signals


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
    plot_output_path: str | pathlib.Path | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Clean raw Neon eye measurements and mapped gaze into aligned signals.

    Align eye measurements with the overlapping gaze timeline, then identify
    unreliable samples using physiological limits, gaze bounds, blinks, and
    abrupt changes. Mask those samples, fill gaps within the configured limits,
    and smooth the resulting signals. Combine the cleaned eye channels into
    pupil diameter, mean eyelid aperture, and aperture difference, and return
    them alongside gaze, elapsed time, and the cleaning settings. Optionally
    plot a raw-versus-cleaned comparison.

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

        plot_output_path (str | pathlib.Path | None):
            PDF destination; defaults to eye_cleaned_vs_raw.pdf in
            processing_neon_path. Used only when verbose is true.

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
        Plotting does not change the returned data or metadata.
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

    # Extract a working copy of the mapped gaze coordinates for bounds filtering.
    bounded_gaze: np.ndarray = gaze_data[:, 1:].copy()

    # First, let's immediately set any weird gaze position to NaN. 
    invalid_gaze_pos_mask: np.ndarray = (
        (bounded_gaze[:, 0] < 0)
        | (bounded_gaze[:, 0] > gaze_x_upper_bound)
        | (bounded_gaze[:, 1] < 0)
        | (bounded_gaze[:, 1] > gaze_y_upper_bound)
    )
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

    # Let's iterate now over the eye states and ensure they are interpolated along the gaze 
    # timestamps
    for idx, column_name in enumerate(eye_channel_names):
        # raw_eye_measurements has the same channel order as eye_channel_names, so idx
        # selects the measurements for column_name (e.g. left pupil diameter).
        channel_values: np.ndarray = raw_eye_measurements[:, idx]

        # Eye and gaze samples may occur at different times. For each target
        # gaze time t, if there is a matching eye timesstamp, 
        # just use that, otherwise, interpolation based on the surrounding eye times t0 and t1 
        interpolated_eye_measurements[:, idx] = np.interp(valid_gaze_timestamps_ns, eye_timestamps_ns, channel_values)

    # Look up the adjacent pupil channels by name. A slice keeps this as a
    # view, so bounds filtering also updates the interpolated measurements.
    pupil_columns: slice = slice(
        eye_channel_names.index("pupil diameter left [mm]"),
        eye_channel_names.index("pupil diameter right [mm]") + 1,
    )
    pupils: np.ndarray = interpolated_eye_measurements[:, pupil_columns]

    # Remove pupil values outside the configured physiological range.
    invalid_pupil_diamater_mask: np.ndarray = (pupils < pupil_diameter_thresholds[0]) | (pupils > pupil_diameter_thresholds[1])
    pupils[invalid_pupil_diamater_mask] = np.nan

    # Keep filtered gaze aligned with the retained output timestamps.
    bounded_gaze = bounded_gaze[valid_gaze_timestamp_mask]

    # Extract the start of each unpadded blink event.
    blink_start_timestamps_ns: np.ndarray = blink_data[:, 0]

    # Extract the corresponding end of each blink event.
    blink_end_timestamps_ns: np.ndarray = blink_data[:, 1]
    if np.any(blink_end_timestamps_ns < blink_start_timestamps_ns):
        raise ValueError("Blink end timestamps must not precede their starts.")

    # Extend each blink interval to include the configured time before and
    # after the event. All interval boundaries remain in nanoseconds.
    padded_blink_starts_ns: np.ndarray = blink_start_timestamps_ns - blink_pre_ns
    padded_blink_ends_ns: np.ndarray = blink_end_timestamps_ns + blink_post_ns

    # Allocate an array that will represent a mask of gaze measurements that are during blinks
    padded_blinks: np.ndarray = np.zeros(len(valid_gaze_timestamps_ns), dtype=bool)

    # Iterate over each blink 
    for start_ns, end_ns in zip(padded_blink_starts_ns, padded_blink_ends_ns):

        # Find the gaze_timestamps idxs that are 
        # within this blink timestamp. Update the padded blink mask 
        # with this infrmation
        padded_blinks |= (
            (valid_gaze_timestamps_ns >= start_ns)
            & (valid_gaze_timestamps_ns <= end_ns)
        )

    # Leave the combined diameter as NaN when both eyes are missing.
    combined_pupil: np.ndarray = np.full(len(pupils), np.nan)
    has_valid_pupil: np.ndarray = np.any(~np.isnan(pupils), axis=1)

    # Average left and right diameters, ignoring a missing eye. If only one
    # eye is available, use its diameter. Skip rows with both eyes missing
    # to avoid an empty-mean warning.
    combined_pupil[has_valid_pupil] = np.nanmean(pupils[has_valid_pupil], axis=1)

    # Average apertures temporarily to detect abrupt physiological changes.
    aperture_columns: slice = slice(
        eye_channel_names.index("eyelid aperture left [mm]"),
        eye_channel_names.index("eyelid aperture right [mm]") + 1,
    )
    apertures: np.ndarray = interpolated_eye_measurements[:, aperture_columns]
    combined_aperture: np.ndarray = np.full(len(apertures), np.nan)
    has_valid_aperture: np.ndarray = np.any(~np.isnan(apertures), axis=1)

    # Average left and right apertures, using the available eye if one is
    # missing. Leave rows with both eyes missing as NaN and skip their mean
    # calculation to avoid an empty-mean warning.
    combined_aperture[has_valid_aperture] = np.nanmean(
        apertures[has_valid_aperture], axis=1
    )

    # The working columns retain each aperture separately through smoothing.
    signal_columns: tuple[str, ...] = (
        "pupil diameter [mm]",
        "eyelid aperture left [mm]",
        "eyelid aperture right [mm]",
        "gaze x [px]",
        "gaze y [px]",
    )
    signals: np.ndarray = np.column_stack((combined_pupil, apertures, bounded_gaze))

    # The combined signals have a different column order from the eye CSV.
    # Resolve their indices from their own names for masking and summaries.
    pupil_column: int = signal_columns.index("pupil diameter [mm]")
    left_aperture_column: int = signal_columns.index("eyelid aperture left [mm]")
    right_aperture_column: int = signal_columns.index("eyelid aperture right [mm]")
    signal_aperture_columns: slice = slice(left_aperture_column, right_aperture_column + 1)
    physiology_columns: slice = slice(pupil_column, right_aperture_column + 1)
    signal_gaze_columns: slice = slice(
        signal_columns.index("gaze x [px]"),
        signal_columns.index("gaze y [px]") + 1,
    )

    # Preserve the pre-mask signals only when a comparison plot is requested.
    before_masking: np.ndarray | None = signals.copy() if verbose else None

    # Measure each sample's absolute change from the previous sample.
    # The first sample has no predecessor, so its change is set to zero.
    absolute_pupil_changes: np.ndarray = np.zeros_like(combined_pupil)
    absolute_pupil_changes[1:] = np.abs(np.diff(combined_pupil))

    absolute_aperture_changes: np.ndarray = np.zeros_like(combined_aperture)
    absolute_aperture_changes[1:] = np.abs(np.diff(combined_aperture))

    # Flag a physiological spike when either measurement exceeds its limit.
    pupil_spike: np.ndarray = absolute_pupil_changes > pupil_change_threshold
    aperture_spike: np.ndarray = absolute_aperture_changes > aperture_change_threshold
    physio_spike: np.ndarray = pupil_spike | aperture_spike

    # Measure absolute horizontal and vertical gaze changes separately, again
    # leaving the first sample's changes at zero.
    absolute_gaze_changes: np.ndarray = np.zeros_like(bounded_gaze)
    absolute_gaze_changes[1:] = np.abs(np.diff(bounded_gaze, axis=0))

    # Combine the horizontal and vertical changes into the straight-line
    # distance moved since the previous sample: sqrt(dx**2 + dy**2).
    # Summing across axis=1 combines x and y for each sample, giving one
    # displacement in pixels. Flag a spike only when this distance exceeds
    # gaze_spike_threshold; this measures movement per sample, not velocity.
    gaze_spikes: np.ndarray = (
        np.sqrt(np.sum(absolute_gaze_changes**2, axis=1)) > gaze_spike_threshold
    )

    # For each Boolean mask below, maximum_filter1d examines a moving window
    # of dilation_samples entries and returns its maximum at each position.
    # For booleans, that means True whenever any entry in the window is True,
    # expanding flagged events to include nearby samples. With a width of 5,
    # an isolated True also flags the two samples before and after it.
    # mode="constant", cval=0 treats positions outside the recording as False,
    # so events do not wrap around or reflect at the recording boundaries.
    # This adds padding in samples, on top of the blink padding in nanoseconds.

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
        gaze_spikes, dilation_samples, mode="constant", cval=0
    )

    # Exclude physiology when either a blink or physiological spike is present.
    invalid_physio_mask: np.ndarray = blink_dilated | physio_spike_dilated

    # Exclude gaze only where a blink and gaze spike coincide.
    invalid_gaze_mask: np.ndarray = blink_dilated & gaze_spike_dilated

    # Mask pupil and both apertures using the physiological exclusion rule.
    signals[invalid_physio_mask, physiology_columns] = np.nan

    # Apply the stricter blink-and-spike exclusion rule to the two gaze columns.
    signals[invalid_gaze_mask, signal_gaze_columns] = np.nan

    # Fill missing values at the beginning and end of each physiology channel
    # with its nearest observed value, matching Geoff's endpoint protection.
    # This avoids PCHIP extrapolation at the recording edges. Interior gaps
    # stay unchanged, and these endpoint fills have no gap-length limit.
    for column in range(physiology_columns.start, physiology_columns.stop):
        # This view lets us update the channel directly in signals.
        channel_values: np.ndarray = signals[:, column]
        valid_sample_mask: np.ndarray = ~np.isnan(channel_values)
        # np.where returns the indices where the mask is True. Take the first
        # result because this channel has only one dimension.
        valid_sample_indices: np.ndarray = np.where(valid_sample_mask)[0]

        # A channel with no observations has no value to carry to its edges.
        if valid_sample_indices.size == 0:
            continue

        first_valid_index: int = int(valid_sample_indices[0])
        last_valid_index: int = int(valid_sample_indices[-1])

        # Fill everything after the last observation with that observation.
        channel_values[last_valid_index + 1 :] = channel_values[last_valid_index]

        # Fill everything before the first observation with that observation.
        channel_values[:first_valid_index] = channel_values[first_valid_index]

    # Start with masked signals and replace only eligible missing samples.
    imputed: np.ndarray = signals.copy()
    for column in range(len(signal_columns)):
        # Use shape-preserving interpolation for physiology and nearest samples for gaze.
        method: str = "nearest" if column in range(signal_gaze_columns.start, signal_gaze_columns.stop) else "pchip"

        # Compute candidate fills using the more permissive blink gap limit.
        longer_fill: np.ndarray = _fill_missing(
            signals[:, column], method, blink_max_gap
        )

        # Apply the general gap limit before allowing blink-specific exceptions.
        imputed[:, column] = _fill_missing(signals[:, column], method, general_max_gap)

        # Select the blink-related positions allowed to use the longer fill limit.
        overwrite: np.ndarray = blink_dilated if column < signal_gaze_columns.start else invalid_gaze_mask

        # Replace only blink-related positions with the longer-gap candidates.
        imputed[overwrite, column] = longer_fill[overwrite]

    # Smooth physiology and gaze with their respective filters.
    smoothed_aperture_still_separate: np.ndarray = _smooth_signals(imputed, smoothing_samples)

    # Calculate mean aperture and left-minus-right from the separate eye channels.
    smoothed_apertures: np.ndarray = smoothed_aperture_still_separate[:, signal_aperture_columns]
    has_valid_smoothed_aperture: np.ndarray = np.any(
        ~np.isnan(smoothed_apertures), axis=1
    )
    combined_aperture: np.ndarray = np.full(len(smoothed_aperture_still_separate), np.nan)

    # Average the available eyes. One missing eye uses the other eye's value;
    # two missing eyes stay NaN without an empty-mean warning.
    combined_aperture[has_valid_smoothed_aperture] = np.nanmean(
        smoothed_apertures[has_valid_smoothed_aperture], axis=1
    )
    aperture_difference: np.ndarray = (
        smoothed_aperture_still_separate[:, left_aperture_column] - smoothed_aperture_still_separate[:, right_aperture_column]
    )

    # Build a separate output matrix with pupil, mean aperture, aperture
    # difference, and gaze x/y. Keep the independently smoothed eyes intact.
    smoothed_final_data: np.ndarray = np.column_stack(
        (
            smoothed_aperture_still_separate[:, pupil_column],
            combined_aperture,
            aperture_difference,
            smoothed_aperture_still_separate[:, signal_gaze_columns],
        )
    )

    # The plotting code in the block below is AI-generated and has not been reviewed.
    if verbose:
        # Count observed channels so rows with no observations remain NaN.
        valid_raw_aperture_counts: np.ndarray = np.sum(~np.isnan(before_masking[:, signal_aperture_columns]), axis=1)
        raw_aperture_mean: np.ndarray = np.divide(
            np.nansum(before_masking[:, signal_aperture_columns], axis=1),
            valid_raw_aperture_counts,
            out=np.full(len(before_masking[:, signal_aperture_columns]), np.nan),
            where=valid_raw_aperture_counts > 0,
        )

        # Match MATLAB's pupil_raw/aperture_*_raw/gaze_*_raw snapshots. These
        # are the aligned, bounded signals before blink/spike censoring.
        plot_raw: np.ndarray = np.column_stack(
            (
                before_masking[:, pupil_column],
                raw_aperture_mean,
                before_masking[:, left_aperture_column] - before_masking[:, right_aperture_column],
                before_masking[:, signal_gaze_columns],
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
            smoothed_final_data,
            raw_blinks,
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
    return np.column_stack((relative_timestamps_s, smoothed_final_data)), metadata


def _plot_matlab_summary(
    time: np.ndarray,
    raw: np.ndarray,
    cleaned: np.ndarray,
    raw_blinks: np.ndarray,
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

        output_path (pathlib.Path):
            Destination for the vector PDF.

    Returns:
        The saved Matplotlib figure, displayed without blocking.

    Notes:
        Panels show pupil diameter, mean aperture, aperture difference, gaze x,
        gaze y, and raw blinks. Values are not clipped for display.
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
    figure.suptitle("Eye Data Processing: Raw vs Cleaned", fontsize=12, fontweight="bold")
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
