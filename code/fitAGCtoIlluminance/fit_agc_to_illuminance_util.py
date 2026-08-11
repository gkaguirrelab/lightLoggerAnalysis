from __future__ import annotations

import itertools
import os
import sys
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matlab
import matlab.engine
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize, TwoSlopeNorm
from numpy.lib.npyio import NpzFile
from scipy.io import loadmat, savemat
from scipy.optimize import OptimizeResult, least_squares
from tqdm.auto import tqdm

THIS_FILE: Path = Path(__file__).resolve()
FIT_AGC_TO_ILLUMINANCE_DIR: Path = THIS_FILE.parent
LIGHT_LOGGER_ANALYSIS_ROOT: Path = THIS_FILE.parents[2]
PROJECTS_ROOT: Path = LIGHT_LOGGER_ANALYSIS_ROOT.parent
LIGHT_LOGGER_ROOT: Path = PROJECTS_ROOT / "lightLogger"
PI_UTILITY_DIR: Path = LIGHT_LOGGER_ROOT / "raspberry_pi_firmware" / "utility"
SENSOR_UTILITY_DIR: Path = LIGHT_LOGGER_ANALYSIS_ROOT / "code" / "library" / "sensor_utility"
CACHED_PROCESSING_DATA_DIR: Path = (
    FIT_AGC_TO_ILLUMINANCE_DIR / "cached_processing_data"
)
VIDEO_IO_UTILITY_DIR: Path = (
    LIGHT_LOGGER_ANALYSIS_ROOT / "code" / "library" / "matlabIO" / "python_libraries"
)
MS_MATLAB_CODE_DIR: Path = LIGHT_LOGGER_ANALYSIS_ROOT / "code" / "analyzeMS" / "matlabCode"
AGC_KERNEL_MAT_PATH: Path = (
    LIGHT_LOGGER_ANALYSIS_ROOT / "misc" / "agc_simulation" / "agc_empirical_kernels.mat"
)
FRAME_SATURATION_DATA_PATH: Path = (
    CACHED_PROCESSING_DATA_DIR
    / "frame_saturation_calibration_data.npz"
)
ILLUMINANCE_DIAGNOSTICS_OUTPUT_DIR: Path = (
    CACHED_PROCESSING_DATA_DIR / "illuminance_channel_diagnostics"
)
LINEAR_SCALE_MAT_OUTPUT_PATH: Path = (
    FIT_AGC_TO_ILLUMINANCE_DIR / "camera_agc_illuminance_linear_scale.mat"
)

for import_path in (PI_UTILITY_DIR, SENSOR_UTILITY_DIR, VIDEO_IO_UTILITY_DIR):
    import_path_str: str = os.fspath(import_path)
    if import_path_str not in sys.path:
        sys.path.append(import_path_str)

from Pi_util import parse_chunks
import world_util
import video_io

_MATLAB_ENGINE: Any | None = None

# Sentinel written into the illuminance signal of a recording that could not
# be converted. Such recordings are excluded from every pooled fit.
ILLUMINANCE_ERROR_SIGNAL: float = -1.0

# Column order of the per-recording summary table, declared once so an empty
# run still returns a frame with the expected schema.
SUMMARY_TABLE_COLUMNS: tuple[str, ...] = (
    "subject_id",
    "activity",
    "errored",
    "correlation",
    "shared_lag_s",
    "slope",
    "intercept",
    "camera_score_std",
    "n_paired_samples",
    "error_message",
    "recording_path",
)


def _get_matlab_engine() -> Any:
    """
    Return a cached MATLAB engine with the minispect analysis code on path.
    """

    global _MATLAB_ENGINE

    # Starting MATLAB is expensive, so keep one engine alive and reuse it for
    # repeated illuminance conversions during multi-recording calibration.
    if _MATLAB_ENGINE is None:
        _MATLAB_ENGINE = matlab.engine.start_matlab()

    return _MATLAB_ENGINE


### HELPER FUNCTION
def ms_counts_to_illuminance(
    ms_counts: np.ndarray,
    diagnostics: bool = False,
    diagnostics_output_path: str | os.PathLike[str] | None = None,
    timestamps: np.ndarray | None = None,
) -> np.ndarray:
    """
    Convert minispect AS counts to illuminance using msCounts2Illuminance.m.

    When ``diagnostics`` is true, MATLAB also calculates the raw first-nine
    channel counts, calibrated per-channel values, returned eight-channel
    lux, legacy nine-channel lux, and their numerical differences. Supplying
    ``diagnostics_output_path`` writes those complete vectors to a compressed
    NumPy archive. The returned array is always the production eight-channel
    estimate.
    """
    global _MATLAB_ENGINE 

    ms_counts = np.asarray(ms_counts, dtype=float)
    if ms_counts.ndim != 2:
        raise ValueError("Minispect counts must be a 2-D sample-by-channel matrix.")

    # Keep the Python entry point, but delegate the calibration logic to the
    # MATLAB implementation so coefficient loading and future fixes stay in
    # one canonical location.
    if(_MATLAB_ENGINE is None):
        _MATLAB_ENGINE = _get_matlab_engine()

        # On the first load in of the MALTAB engine, we need to have it tbUse our projcet
        _MATLAB_ENGINE.tbUseProject('lightLoggerAnalysis', nargout=0)

    matlab_engine: object = _MATLAB_ENGINE    
    matlab_counts: matlab.double = matlab.double(ms_counts.tolist())
    illuminance: Any
    if diagnostics:
        matlab_diagnostics: Any
        illuminance, matlab_diagnostics = matlab_engine.msCounts2Illuminance(
            matlab_counts,
            "diagnostics",
            True,
            "verbose",
            False,
            nargout=2,
        )

        # Save all comparison vectors immediately so later filtering and
        # plotting cannot alter the values being tested.
        if diagnostics_output_path is not None:
            _write_illuminance_channel_diagnostics(
                diagnostics_output_path,
                matlab_diagnostics,
                timestamps,
            )
    else:
        illuminance = matlab_engine.msCounts2Illuminance(
            matlab_counts,
            nargout=1,
        )

    return np.asarray(illuminance, dtype=float).reshape(-1)


def _matlab_struct_field(matlab_struct: Any, field_name: str) -> Any:
    """Return one field from a MATLAB scalar struct returned by the engine."""

    if isinstance(matlab_struct, dict):
        return matlab_struct[field_name]

    try:
        return matlab_struct[field_name]
    except (KeyError, TypeError):
        return getattr(matlab_struct, field_name)


def _write_illuminance_channel_diagnostics(
    output_path: str | os.PathLike[str],
    matlab_diagnostics: Any,
    timestamps: np.ndarray | None,
) -> None:
    """Write explicit eight- and nine-channel calculations to an archive."""

    available: bool = bool(
        np.asarray(
            _matlab_struct_field(matlab_diagnostics, "available")
        ).squeeze()
    )
    if not available:
        raise RuntimeError(
            "MATLAB could not calculate the nine-channel illuminance diagnostics."
        )

    raw_counts: np.ndarray = np.asarray(
        _matlab_struct_field(matlab_diagnostics, "raw_counts"),
        dtype=float,
    )
    illuminance_by_channel: np.ndarray = np.asarray(
        _matlab_struct_field(matlab_diagnostics, "illuminance_by_channel"),
        dtype=float,
    )
    lux_8_channels: np.ndarray = np.asarray(
        _matlab_struct_field(matlab_diagnostics, "lux_8_channels"),
        dtype=float,
    ).reshape(-1)
    lux_9_channels: np.ndarray = np.asarray(
        _matlab_struct_field(matlab_diagnostics, "lux_9_channels"),
        dtype=float,
    ).reshape(-1)
    difference_lux: np.ndarray = np.asarray(
        _matlab_struct_field(
            matlab_diagnostics,
            "difference_8_minus_9_lux",
        ),
        dtype=float,
    ).reshape(-1)
    maximum_absolute_difference_lux: float = float(
        np.asarray(
            _matlab_struct_field(
                matlab_diagnostics,
                "maximum_absolute_difference_lux",
            )
        ).squeeze()
    )
    maximum_relative_difference: float = float(
        np.asarray(
            _matlab_struct_field(
                matlab_diagnostics,
                "maximum_relative_difference",
            )
        ).squeeze()
    )
    exactly_identical: bool = bool(
        np.asarray(
            _matlab_struct_field(matlab_diagnostics, "exactly_identical")
        ).squeeze()
    )
    numerically_identical: bool = bool(
        np.asarray(
            _matlab_struct_field(matlab_diagnostics, "numerically_identical")
        ).squeeze()
    )
    comparison_tolerance_lux: float = float(
        np.asarray(
            _matlab_struct_field(
                matlab_diagnostics,
                "comparison_tolerance_lux",
            )
        ).squeeze()
    )

    # A missing timestamp vector is represented explicitly rather than by a
    # fabricated sample rate. Normal fit calls always provide real AS times.
    if timestamps is None:
        timestamp_array: np.ndarray = np.full(
            lux_8_channels.shape,
            np.nan,
            dtype=float,
        )
    else:
        timestamp_array = np.asarray(timestamps, dtype=float).reshape(-1)
        if timestamp_array.size != lux_8_channels.size:
            raise ValueError(
                "Diagnostic timestamps and illuminance vectors must have "
                "the same length."
            )

    destination: Path = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Store each option as a self-contained calculation: the raw inputs, each
    # independently calibrated channel, and the resulting across-channel mean.
    # Explicit arrays avoid requiring later analysis to infer the 8-channel
    # option by slicing arrays labeled for the 9-channel option.
    np.savez_compressed(
        destination,
        timestamp=timestamp_array,
        raw_counts_8_channels=raw_counts[:, :8],
        converted_lux_per_channel_8_channels=illuminance_by_channel[:, :8],
        average_lux_8_channels=lux_8_channels,
        raw_counts_9_channels=raw_counts[:, :9],
        converted_lux_per_channel_9_channels=illuminance_by_channel[:, :9],
        average_lux_9_channels=lux_9_channels,
        difference_8_minus_9_lux=difference_lux,
        maximum_absolute_difference_lux=maximum_absolute_difference_lux,
        maximum_relative_difference=maximum_relative_difference,
        exactly_identical=exactly_identical,
        numerically_identical=numerically_identical,
        comparison_tolerance_lux=comparison_tolerance_lux,
    )


def _load_nine_channel_illuminance(
    diagnostic_path: str | os.PathLike[str],
) -> np.ndarray:
    """Load the complete legacy nine-channel average from a diagnostic archive."""

    with np.load(diagnostic_path) as diagnostic_data:
        return np.asarray(
            diagnostic_data["average_lux_9_channels"],
            dtype=float,
        ).copy()


def _activity_point_cloud_from_model(
    recording_result: dict[str, Any],
    model_result: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Build the camera-score and filtered-lux vectors used by the activity plot."""

    as_time: np.ndarray = np.asarray(recording_result["as_time"], dtype=float)
    camera_time: np.ndarray = np.asarray(
        recording_result["camera_time"],
        dtype=float,
    )
    camera_score: np.ndarray = np.asarray(
        recording_result["camera_score"],
        dtype=float,
    )

    # Apply this option's shared lag in physical lux, then align the camera
    # score to the same minispect timestamps used by the final scatter plot.
    filtered_illuminance: np.ndarray = np.interp(
        as_time,
        model_result["time"],
        model_result["filtered_illuminance"],
        left=np.nan,
        right=np.nan,
    )
    camera_score_on_as_time: np.ndarray = np.interp(
        as_time,
        camera_time,
        camera_score,
        left=np.nan,
        right=np.nan,
    )
    valid: np.ndarray = (
        np.isfinite(camera_score_on_as_time)
        & np.isfinite(filtered_illuminance)
    )

    return camera_score_on_as_time[valid], filtered_illuminance[valid]


def _append_activity_point_cloud_diagnostics(
    diagnostic_path: str | os.PathLike[str],
    camera_score_8_channels: np.ndarray,
    filtered_illuminance_8_channels: np.ndarray,
    shared_lag_8_channels: float,
    correlation_8_channels: float,
    camera_score_9_channels: np.ndarray,
    filtered_illuminance_9_channels: np.ndarray,
    shared_lag_9_channels: float,
    correlation_9_channels: float,
    minimum_correlation: float,
) -> None:
    """Append both activity-plot point clouds to one recording archive."""

    destination: Path = Path(diagnostic_path)
    with np.load(destination) as diagnostic_data:
        payload: dict[str, np.ndarray] = {
            key: diagnostic_data[key].copy()
            for key in diagnostic_data.files
        }

    payload.update(
        {
            "activity_point_cloud_camera_score_8_channels": np.asarray(
                camera_score_8_channels,
                dtype=float,
            ),
            "activity_point_cloud_filtered_illuminance_8_channels": np.asarray(
                filtered_illuminance_8_channels,
                dtype=float,
            ),
            "activity_point_cloud_shared_lag_seconds_8_channels": np.asarray(
                shared_lag_8_channels,
                dtype=float,
            ),
            "activity_point_cloud_correlation_8_channels": np.asarray(
                correlation_8_channels,
                dtype=float,
            ),
            "activity_point_cloud_included_8_channels": np.asarray(
                np.isfinite(correlation_8_channels)
                and correlation_8_channels >= minimum_correlation,
                dtype=bool,
            ),
            "activity_point_cloud_camera_score_9_channels": np.asarray(
                camera_score_9_channels,
                dtype=float,
            ),
            "activity_point_cloud_filtered_illuminance_9_channels": np.asarray(
                filtered_illuminance_9_channels,
                dtype=float,
            ),
            "activity_point_cloud_shared_lag_seconds_9_channels": np.asarray(
                shared_lag_9_channels,
                dtype=float,
            ),
            "activity_point_cloud_correlation_9_channels": np.asarray(
                correlation_9_channels,
                dtype=float,
            ),
            "activity_point_cloud_included_9_channels": np.asarray(
                np.isfinite(correlation_9_channels)
                and correlation_9_channels >= minimum_correlation,
                dtype=bool,
            ),
            "activity_point_cloud_minimum_correlation": np.asarray(
                minimum_correlation,
                dtype=float,
            ),
        }
    )

    # Replace the archive atomically so an interrupted write cannot leave a
    # partially valid diagnostic file.
    temporary_path: Path = destination.with_suffix(".tmp.npz")
    np.savez_compressed(temporary_path, **payload)
    os.replace(temporary_path, destination)


### HELPER FUNCTION
def safe_ms_counts_to_illuminance(
    ms_counts: np.ndarray,
    diagnostics: bool = False,
    diagnostics_output_path: str | os.PathLike[str] | None = None,
    timestamps: np.ndarray | None = None,
) -> tuple[np.ndarray, bool, str]:
    """
    Convert AS counts to illuminance, flagging failures instead of raising.

    msCounts2Illuminance takes log10 of the first eight spectral channels, so
    any sample that is not strictly positive in those channels cannot be
    converted. Rather than aborting the whole calibration run, the
    unconvertible portion is filled with
    ILLUMINANCE_ERROR_SIGNAL and the recording is marked as errored so the
    caller can exclude it from every pooled fit.

    When ``diagnostics`` is true, the MATLAB conversion calculates an
    additional comparison against the legacy nine-channel calculation and
    writes it when an output path is supplied. Returns the illuminance signal,
    whether the recording errored, and a message describing the failure.
    """

    ms_counts = np.asarray(ms_counts, dtype=float)

    if ms_counts.ndim != 2:
        raise ValueError("Minispect counts must be a 2-D sample-by-channel matrix.")
    if ms_counts.shape[1] < 8:
        raise ValueError(
            "Minispect counts must contain at least the first 8 spectral channels."
        )

    n_samples: int = ms_counts.shape[0]
    illuminance: np.ndarray = np.full(
        n_samples,
        ILLUMINANCE_ERROR_SIGNAL,
        dtype=float,
    )

    # Clear and NIR do not contribute to illuminance, so only AS_0-AS_7 can
    # make a sample unconvertible.
    spectral_counts: np.ndarray = ms_counts[:, :8]
    invalid_samples: np.ndarray = ~np.all(
        np.isfinite(spectral_counts) & (spectral_counts > 0),
        axis=1,
    )

    if not np.any(invalid_samples):
        # The counts look convertible, but MATLAB remains the authority on
        # whether they actually are, so a failure there is still caught.
        try:
            return (
                ms_counts_to_illuminance(
                    ms_counts,
                    diagnostics,
                    diagnostics_output_path,
                    timestamps,
                ),
                False,
                "",
            )
        except Exception as conversion_error:
            return illuminance, True, str(conversion_error).strip()

    first_invalid_index: int = int(np.argmax(invalid_samples))
    n_invalid: int = int(np.sum(invalid_samples))
    message: str = (
        f"{n_invalid} of {n_samples} minispect samples are not strictly "
        f"positive; first at sample index {first_invalid_index}."
    )

    # Convert whatever prefix is valid so the failure point stays visible in
    # the data, and flag everything from the first bad sample onward.
    if first_invalid_index >= 1:
        try:
            illuminance[:first_invalid_index] = ms_counts_to_illuminance(
                ms_counts[:first_invalid_index],
                diagnostics,
                diagnostics_output_path,
                (
                    None
                    if timestamps is None
                    else np.asarray(timestamps)[:first_invalid_index]
                ),
            )
        except Exception as conversion_error:
            message = (
                f"{message} Converting the valid prefix also failed: "
                f"{str(conversion_error).strip()}"
            )

    return illuminance, True, message


def flatten_ms_chunks(ms_chunks: Sequence[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    """
    Flatten minispect AS values and timestamps from parsed recording chunks.
    """

    if len(ms_chunks) == 0:
        raise RuntimeError("No minispect chunks were parsed.")

    as_pairs: list[tuple[np.ndarray, np.ndarray]] = list(
        filter(
            lambda as_pair: len(as_pair[0]) > 0,
            map(
                lambda chunk: (chunk["M"]["v"]["AS"], chunk["M"]["t"]["AS"]),
                ms_chunks,
            ),
        )
    )

    if len(as_pairs) == 0:
        raise RuntimeError("No AS data were found in the minispect chunks.")

    mismatched_chunks: list[tuple[np.ndarray, np.ndarray]] = list(
        filter(
            lambda as_pair: len(as_pair[0]) != len(as_pair[1]),
            as_pairs,
        )
    )

    if len(mismatched_chunks) > 0:
        raise ValueError(
            "AS value and timestamp lengths do not match within a minispect chunk."
        )

    as_value_chunks: tuple[np.ndarray, ...]
    as_time_chunks: tuple[np.ndarray, ...]
    as_value_chunks, as_time_chunks = zip(*as_pairs)

    as_values: np.ndarray = np.vstack(as_value_chunks)
    as_time: np.ndarray = np.concatenate(as_time_chunks).astype(float)

    return as_values, as_time


def extract_agc_columns(metadata: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract analog gain, digital gain, and exposure from known AGC schemas.
    """

    legacy_columns: tuple[str, str, str] = ("Again", "Dgain", "exposure")
    modern_columns: tuple[str, str, str] = (
        "cameraAgain",
        "AGCDgain",
        "cameraExposure",
    )

    if all(column in metadata.columns for column in legacy_columns):
        analog_gain_column: str = "Again"
        digital_gain_column: str = "Dgain"
        exposure_column: str = "exposure"
    elif all(column in metadata.columns for column in modern_columns):
        analog_gain_column = "cameraAgain"
        digital_gain_column = "AGCDgain"
        exposure_column = "cameraExposure"
    else:
        raise ValueError(
            "World metadata must contain either "
            "['timestamp', 'Again', 'Dgain', 'exposure'] or "
            "['timestamp', 'cameraAgain', 'AGCDgain', 'cameraExposure', "
            "'AGCAgain', 'AGCExposure']."
        )

    analog_gain: np.ndarray = metadata[analog_gain_column].to_numpy(dtype=float)
    digital_gain: np.ndarray = metadata[digital_gain_column].to_numpy(dtype=float)
    exposure: np.ndarray = metadata[exposure_column].to_numpy(dtype=float)

    return analog_gain, digital_gain, exposure


### HELPER FUNCTION
def apply_empirical_kernel(
    ms_time: np.ndarray,
    ms_signal: np.ndarray,
    kernel_mat_path: str,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply the AGC-derived empirical kernel to a minispect signal.
    """

    mat: dict[str, Any] = loadmat(kernel_mat_path)

    kernel_time: np.ndarray = np.squeeze(mat["commonKernelTime"]).astype(float)
    kernel: np.ndarray = np.squeeze(mat["meanKernel"]).astype(float)

    # Put minispect signal on the kernel's 4-Hz timebase.
    dt_kernel: np.float64 = np.median(np.diff(kernel_time))
    t_uniform: np.ndarray = np.arange(ms_time[0], ms_time[-1], dt_kernel)
    ms_uniform: np.ndarray = np.interp(t_uniform, ms_time, ms_signal)

    # Convert continuous kernel values to discrete convolution weights.
    kernel_weights: np.ndarray = kernel / np.sum(kernel)

    # Pad with the initial signal value so the filter does not assume
    # that the minispect signal was zero before the recording began.
    n_kernel: int = len(kernel_weights)
    padded_signal: np.ndarray = np.pad(ms_uniform, (n_kernel - 1, 0), mode="edge")

    filtered: np.ndarray = np.convolve(padded_signal, kernel_weights, mode="valid")

    return t_uniform, filtered

### HELPER FUNCTION
def characterize_empirical_kernel(
    camera_time: np.ndarray,
    camera_score: np.ndarray,
    ms_time: np.ndarray,
    ms_signal: np.ndarray,
    kernel_mat_path: str,
    lag_range: tuple[float, float] = (-10, 10),
    n_lags: int = 201,
    fixed_lag: float | None = None,
) -> dict[str, Any]:
    """
    Apply the AGC-derived empirical kernel and find the temporal lag that
    maximizes correlation with the brightness-oriented camera AGC product.

    Supplying fixed_lag reports that lag instead of the best-correlating
    one. The kernel shape is already identical across recordings, so fixing
    the lag makes the whole empirical model identical as well.
    """

    empirical_time: np.ndarray
    unshifted_filtered_illuminance: np.ndarray
    empirical_time, unshifted_filtered_illuminance = apply_empirical_kernel(
        ms_time, ms_signal, kernel_mat_path
    )

    # The camera score is the AGC product, which decreases as illumination
    # increases. Negate only its standardized temporal trace so a well-aligned
    # illuminance prediction retains a positive correlation quality metric.
    camera_z: np.ndarray = -zscore_signal(camera_score)
    unshifted_prediction: np.ndarray = zscore_signal(
        unshifted_filtered_illuminance
    )

    # Put observed camera score onto the empirical kernel timebase.
    camera_on_empirical_time: np.ndarray = np.interp(
        empirical_time,
        camera_time,
        camera_z,
        left=np.nan,
        right=np.nan
    )

    candidate_lags: np.ndarray = np.linspace(lag_range[0], lag_range[1], n_lags)
    correlations: np.ndarray = np.full(len(candidate_lags), np.nan)

    for lag_idx, lag in enumerate(candidate_lags):

        # Positive lag means the camera response follows the minispect signal.
        predicted_camera: np.ndarray = np.interp(
            empirical_time - lag,
            empirical_time,
            unshifted_prediction,
            left=np.nan,
            right=np.nan
        )

        valid: np.ndarray = (
            np.isfinite(camera_on_empirical_time) &
            np.isfinite(predicted_camera)
        )

        if np.sum(valid) < 3:
            continue

        correlations[lag_idx] = np.corrcoef(
            camera_on_empirical_time[valid],
            predicted_camera[valid]
        )[0, 1]

    best_idx: int | None = _select_lag_index(
        correlations,
        candidate_lags,
        fixed_lag,
    )

    if best_idx is None:
        return {
            "time": empirical_time,
            "camera": camera_on_empirical_time,
            "illuminance": unshifted_prediction,
            "prediction": np.full(empirical_time.shape, np.nan, dtype=float),
            "filtered_illuminance": np.full(
                empirical_time.shape,
                np.nan,
                dtype=float,
            ),
            "unshifted_filtered_illuminance": unshifted_filtered_illuminance,
            "best_lag": np.nan,
            "best_correlation": np.nan,
            "rmse": np.nan,
            "candidate_lags": candidate_lags,
            "correlations": correlations,
        }

    best_lag: np.float64 = candidate_lags[best_idx]
    best_correlation: np.float64 = correlations[best_idx]

    best_prediction: np.ndarray = np.interp(
        empirical_time - best_lag,
        empirical_time,
        unshifted_prediction,
        left=np.nan,
        right=np.nan
    )
    filtered_illuminance: np.ndarray = np.interp(
        empirical_time - best_lag,
        empirical_time,
        unshifted_filtered_illuminance,
        left=np.nan,
        right=np.nan,
    )

    valid = (
        np.isfinite(camera_on_empirical_time) &
        np.isfinite(best_prediction)
    )

    rmse: np.float64 = np.sqrt(np.mean(
        (camera_on_empirical_time[valid] - best_prediction[valid]) ** 2
    ))

    return {
        "time": empirical_time,
        "camera": camera_on_empirical_time,
        "illuminance": unshifted_prediction,
        "prediction": best_prediction,
        # Apply the same lag in physical units so correlation, scatter, and
        # calibration all describe exactly the same temporal model.
        "filtered_illuminance": filtered_illuminance,
        "unshifted_filtered_illuminance": unshifted_filtered_illuminance,
        "best_lag": best_lag,
        "best_correlation": best_correlation,
        "rmse": rmse,
        "candidate_lags": candidate_lags,
        "correlations": correlations
    }


def characterize_shared_empirical_kernel(
    recordings: Sequence[dict[str, Any]],
    kernel_mat_path: str,
    lag_range: tuple[float, float] = (-10, 10),
    n_lags: int = 201,
    fixed_lag: float | None = None,
) -> dict[str, Any]:
    """
    Apply one empirical AGC kernel and one shared lag to all recordings.

    The empirical kernel already fixes the temporal response shape. When no
    lag is supplied, this function chooses the single lag that maximizes the
    mean correlation across recordings. Supplying ``fixed_lag`` forces that
    lag instead. Either way, every returned recording uses the same kernel
    and lag.
    """

    # First calculate each recording's complete lag-correlation curve. The
    # temporary per-recording optimum is discarded after the shared lag is
    # selected; it is needed only to avoid duplicating kernel convolution.
    provisional_results: list[dict[str, Any]] = [
        characterize_empirical_kernel(
            recording["camera_time"],
            recording["camera_score"],
            recording["as_time"],
            recording["illuminance"],
            kernel_mat_path,
            lag_range=lag_range,
            n_lags=n_lags,
        )
        for recording in recordings
    ]

    candidate_lags: np.ndarray = provisional_results[0]["candidate_lags"]
    correlations: np.ndarray = np.vstack(
        [result["correlations"] for result in provisional_results]
    )
    mean_correlation_by_lag: np.ndarray = np.nanmean(correlations, axis=0)
    shared_lag_index: int | None = _select_lag_index(
        mean_correlation_by_lag,
        candidate_lags,
        fixed_lag,
    )

    if shared_lag_index is None:
        raise RuntimeError(
            "No candidate lag produced overlapping empirical-kernel and "
            "camera-score samples."
        )

    shared_lag: float = float(candidate_lags[shared_lag_index])
    per_recording_results: list[dict[str, Any]] = []

    # Rebuild only the inexpensive shifted arrays at the chosen lag. This
    # keeps each recording's physical and normalized predictions aligned.
    for provisional_result in provisional_results:
        empirical_time: np.ndarray = provisional_result["time"]
        unshifted_prediction: np.ndarray = provisional_result["illuminance"]
        unshifted_filtered_illuminance: np.ndarray = provisional_result[
            "unshifted_filtered_illuminance"
        ]
        prediction: np.ndarray = np.interp(
            empirical_time - shared_lag,
            empirical_time,
            unshifted_prediction,
            left=np.nan,
            right=np.nan,
        )
        filtered_illuminance: np.ndarray = np.interp(
            empirical_time - shared_lag,
            empirical_time,
            unshifted_filtered_illuminance,
            left=np.nan,
            right=np.nan,
        )
        camera: np.ndarray = provisional_result["camera"]
        valid: np.ndarray = np.isfinite(camera) & np.isfinite(prediction)
        correlation: float = (
            float(np.corrcoef(camera[valid], prediction[valid])[0, 1])
            if np.sum(valid) >= 3
            else np.nan
        )
        rmse: float = (
            float(np.sqrt(np.mean((camera[valid] - prediction[valid]) ** 2)))
            if np.any(valid)
            else np.nan
        )

        per_recording_results.append(
            {
                "time": empirical_time,
                "camera": camera,
                "illuminance": unshifted_prediction,
                "prediction": prediction,
                "filtered_illuminance": filtered_illuminance,
                "best_lag": shared_lag,
                "best_correlation": correlation,
                "rmse": rmse,
                "candidate_lags": candidate_lags,
                "correlations": provisional_result["correlations"],
            }
        )

    return {
        "shared_lag": shared_lag,
        "candidate_lags": candidate_lags,
        "mean_correlation_by_lag": mean_correlation_by_lag,
        "correlations": correlations,
        "per_recording_results": per_recording_results,
    }

### HELPER FUNCTION 
def zscore_signal(x: np.ndarray) -> np.ndarray:
    """
    Normalize a signal to mean 0 and standard deviation 1. This allows us to compare temporal shape even though the camera
    score and minispect channel have different units / scales.
    """
    x = np.asarray(x, dtype=float)

    mean_x: np.float64 = np.nanmean(x)
    std_x: np.float64 = np.nanstd(x)

    if std_x == 0:
        raise ValueError("Cannot z-score a signal with zero variance.")

    return (x - mean_x) / std_x

def _select_lag_index(
    correlations: np.ndarray,
    candidate_lags: np.ndarray,
    fixed_lag: float | None,
) -> int | None:
    """
    Choose which candidate lag to report for one recording.

    With no fixed lag this is the best-correlating lag. With a fixed lag the
    nearest candidate is reported instead, so every recording is forced onto
    the same temporal model rather than each choosing its own alignment.
    """

    if fixed_lag is not None:
        return int(np.argmin(np.abs(candidate_lags - fixed_lag)))

    # Every lag can be NaN when two signals never overlap, and nanargmax
    # raises rather than returning a sentinel in that case.
    if not np.any(np.isfinite(correlations)):
        return None

    return int(np.nanargmax(correlations))


def _sample_calibration_points(
    recording_result: dict[str, Any],
    sample_interval_seconds: float = 1.0,
    filtered_illuminance_key: str = "filtered_illuminance",
) -> dict[str, np.ndarray]:
    """
    Sample a selected filtered illuminance signal and camera score sparsely.
    """

    as_time: np.ndarray = np.asarray(recording_result["as_time"], dtype=float)
    filtered_illuminance: np.ndarray = np.asarray(
        recording_result[filtered_illuminance_key],
        dtype=float,
    )
    camera_score_on_as_time: np.ndarray = np.asarray(
        recording_result["camera_score_on_as_time"],
        dtype=float,
    )

    valid: np.ndarray = (
        np.isfinite(as_time)
        & np.isfinite(filtered_illuminance)
        & np.isfinite(camera_score_on_as_time)
    )

    if np.sum(valid) < 2:
        return {
            "time": np.array([], dtype=float),
            "filtered_illuminance": np.array([], dtype=float),
            "camera_score": np.array([], dtype=float),
        }

    valid_time: np.ndarray = as_time[valid]
    valid_illuminance: np.ndarray = filtered_illuminance[valid]
    valid_camera_score: np.ndarray = camera_score_on_as_time[valid]
    time_order: np.ndarray = np.argsort(valid_time)
    valid_time = valid_time[time_order]
    valid_illuminance = valid_illuminance[time_order]
    valid_camera_score = valid_camera_score[time_order]

    sample_time: np.ndarray = np.arange(
        valid_time[0],
        valid_time[-1],
        sample_interval_seconds,
    )

    return {
        "time": sample_time,
        "filtered_illuminance": np.interp(
            sample_time,
            valid_time,
            valid_illuminance,
        ),
        "camera_score": np.interp(
            sample_time,
            valid_time,
            valid_camera_score,
        ),
    }


def fit_steady_state_illuminance_calibration(
    results_by_recording: Sequence[dict[str, Any]],
    sample_interval_seconds: float = 1.0,
    filtered_illuminance_key: str = "filtered_illuminance",
) -> dict[str, Any]:
    """
    Fit one pooled robust illuminance-from-camera-score calibration line.

    Recordings are sampled at a fixed interval so long recordings do not
    dominate the fit simply because they contain more sensor samples.
    """

    sampled_results: list[dict[str, np.ndarray]] = [
        _sample_calibration_points(
            recording_result,
            sample_interval_seconds,
            filtered_illuminance_key,
        )
        for recording_result in results_by_recording
    ]

    filtered_illuminance: np.ndarray = np.concatenate(
        [sampled_result["filtered_illuminance"] for sampled_result in sampled_results]
    )
    camera_score: np.ndarray = np.concatenate(
        [sampled_result["camera_score"] for sampled_result in sampled_results]
    )

    valid: np.ndarray = np.isfinite(filtered_illuminance) & np.isfinite(camera_score)
    if np.sum(valid) < 2:
        raise RuntimeError(
            "At least two valid paired calibration samples are required "
            "to fit camera AGC product against filtered illuminance."
        )

    line_fit: dict[str, Any] = fit_illuminance_on_camera_score(
        camera_score[valid],
        filtered_illuminance[valid],
    )
    slope: float = line_fit["slope"]
    intercept: float = line_fit["intercept"]
    fitted_illuminance: np.ndarray = (
        intercept + slope * camera_score[valid]
    )
    residuals: np.ndarray = filtered_illuminance[valid] - fitted_illuminance
    rmse: np.float64 = np.sqrt(np.mean(residuals ** 2))
    correlation: np.float64 = np.corrcoef(
        filtered_illuminance[valid],
        camera_score[valid],
    )[0, 1]

    return {
        "sample_interval_seconds": sample_interval_seconds,
        "filtered_illuminance_key": filtered_illuminance_key,
        "sampled_results": sampled_results,
        "filtered_illuminance": filtered_illuminance[valid],
        "camera_score": camera_score[valid],
        "slope": slope,
        "intercept": intercept,
        "rmse": rmse,
        "correlation": correlation,
        "fit_method": line_fit["fit_method"],
        "residual_scale": line_fit["residual_scale"],
        "implied_illuminance_equation": (
            "illuminance = intercept + slope * camera_score"
        ),
    }


def _draw_figure_now(figure: Any) -> None:
    """
    Render a figure as soon as it is built rather than deferring the draw.

    A notebook normally holds every figure until the whole cell finishes, so
    a long run looks like it has hung with nothing to show. Drawing each
    figure the moment it is complete makes progress visible while the run is
    still going.
    """

    try:
        from IPython import get_ipython
        from IPython.display import display

        if get_ipython() is not None:
            display(figure)
            # The figure has already been shown, so close it to stop the
            # notebook from drawing a second copy when the cell finishes.
            plt.close(figure)
            return
    except ImportError:
        pass

    # Outside a notebook, nudge the event loop so an interactive backend
    # paints now. Headless backends have nothing to paint into, and warn
    # rather than raise, so both outcomes are swallowed here.
    try:
        figure.canvas.draw_idle()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plt.pause(0.001)
    except Exception:
        pass


def _set_legend_top_right(ax: Axes, **kwargs: Any) -> Axes:
    """
    Place an axis legend inside the top-right corner when labels exist.
    """

    handles, labels = ax.get_legend_handles_labels()
    if len(handles) == 0:
        return ax

    ax.legend(
        handles,
        labels,
        loc="upper right",
        **kwargs,
    )
    return ax


def _plot_camera_settings(
    ax: Axes,
    camera_time_relative: np.ndarray,
    analog_gain: np.ndarray,
    digital_gain: np.ndarray,
    exposure: np.ndarray,
) -> Axes:
    """
    Plot camera gain and exposure settings with separate y-axes.
    """

    exposure_axis: Axes = ax.twinx()

    ax.plot(camera_time_relative, analog_gain, label="Analog gain", color="tab:blue")
    ax.plot(camera_time_relative, digital_gain, label="Digital gain", color="tab:cyan")
    exposure_axis.plot(
        camera_time_relative,
        exposure,
        label="Exposure",
        color="tab:orange",
    )

    ax.set_xlabel("Time from recording start (s)")
    ax.set_ylabel("Gain")
    exposure_axis.set_ylabel("Exposure")
    ax.set_ylim(1, 11)
    exposure_axis.set_ylim(25, 8500)
    ax.set_title("Camera Settings")
    ax.grid()

    handles_left, labels_left = ax.get_legend_handles_labels()
    handles_right, labels_right = exposure_axis.get_legend_handles_labels()
    ax.legend(
        handles_left + handles_right,
        labels_left + labels_right,
        loc="upper right",
        fontsize=8,
    )
    return ax


def _plot_camera_light_score(
    ax: Axes,
    camera_time_relative: np.ndarray,
    camera_score: np.ndarray,
) -> Axes:
    """
    Plot the camera AGC product on an existing axis.
    """

    ax.plot(camera_time_relative, camera_score, label="Camera AGC product")
    ax.set_xlabel("Time from recording start (s)")
    ax.set_ylabel("Again x Dgain x exposure")
    ax.set_title("Camera AGC product")
    ax.grid()
    _set_legend_top_right(ax, fontsize=8)
    return ax


def _plot_video_saturation(
    ax: Axes,
    saturation_time_relative: np.ndarray,
    video_saturation: np.ndarray,
) -> Axes:
    """Plot the saturated-pixel proportion for each decoded world-video frame."""

    ax.plot(
        saturation_time_relative,
        video_saturation,
        label="Saturated pixels",
    )
    ax.set_xlabel("Time from video start (s)")
    ax.set_ylabel("Saturated pixel proportion")
    ax.set_title("World video saturation")
    ax.set_ylim(-0.1, 1.1)
    ax.grid()
    _set_legend_top_right(ax, fontsize=8)
    return ax


def _plot_all_video_saturation(
    results_by_recording: Sequence[dict[str, Any]],
) -> Axes:
    """Plot percent saturation over time for every processed recording."""

    figure: Any
    ax: Axes
    figure, ax = plt.subplots(figsize=(12, 7))

    # Convert the stored 0-1 proportions to percentages and retain a separate
    # labeled line for every supplied recording.
    for recording_result in results_by_recording:
        recording_label: str = _recording_name_from_path(
            recording_result["recording_path"]
        )
        saturation_percent: np.ndarray = (
            np.asarray(recording_result["video_saturation"], dtype=float)
            * 100.0
        )
        ax.plot(
            recording_result["saturation_time_relative"],
            saturation_percent,
            label=recording_label,
        )

    ax.set_title("Saturation over Time")
    ax.set_xlabel("Time since start (s)")
    ax.set_ylabel("Percent saturation")
    ax.set_ylim(-10.0, 110.0)
    ax.grid()
    _set_legend_top_right(ax, fontsize=9)
    figure.tight_layout()
    _draw_figure_now(figure)
    return ax


def _video_saturation_at_times_percent(
    sample_time_relative: np.ndarray,
    saturation_time_relative: np.ndarray,
    video_saturation: np.ndarray,
) -> np.ndarray:
    """Return nearest-frame spatial saturation percentages at sample times."""

    # Validate the frame-level time and saturation vectors before indexing so
    # malformed video metadata cannot silently color the wrong observations.
    sample_times: np.ndarray = np.asarray(sample_time_relative, dtype=float)
    frame_times: np.ndarray = np.asarray(saturation_time_relative, dtype=float)
    saturation_values: np.ndarray = np.asarray(video_saturation, dtype=float)
    if frame_times.ndim != 1 or saturation_values.ndim != 1:
        raise ValueError("Video saturation times and values must be one-dimensional.")
    if frame_times.size != saturation_values.size:
        raise ValueError(
            "Video saturation times and values must have the same length."
        )
    if frame_times.size == 0:
        return np.full(sample_times.shape, np.nan, dtype=float)

    # Locate the frame immediately to either side of every sample, then retain
    # whichever frame timestamp is nearest. Clipping maps tiny endpoint timing
    # differences to the first or last decoded frame.
    insertion_indices: np.ndarray = np.searchsorted(frame_times, sample_times)
    right_indices: np.ndarray = np.clip(
        insertion_indices,
        0,
        frame_times.size - 1,
    )
    left_indices: np.ndarray = np.clip(
        insertion_indices - 1,
        0,
        frame_times.size - 1,
    )
    use_right: np.ndarray = (
        np.abs(frame_times[right_indices] - sample_times)
        < np.abs(sample_times - frame_times[left_indices])
    )
    nearest_indices: np.ndarray = np.where(
        use_right,
        right_indices,
        left_indices,
    )

    return 100.0 * saturation_values[nearest_indices]


def _plot_raw_as_channels(
    ax: Axes,
    as_time_relative: np.ndarray,
    as_values: np.ndarray,
) -> Axes:
    """
    Plot each flattened minispect AS channel on an existing axis.
    """

    for column_number in range(as_values.shape[1]):
        ax.plot(
            as_time_relative,
            as_values[:, column_number],
            label=f"AS {column_number}",
        )

    ax.set_xlabel("Time from AS recording start (s)")
    ax.set_ylabel("AS sensor value")
    ax.set_title("Flattened minispect AS channels")
    _set_legend_top_right(ax, fontsize=7)
    ax.grid()
    return ax


def _plot_illuminance(
    ax: Axes,
    as_time_relative: np.ndarray,
    illuminance: np.ndarray,
) -> Axes:
    """
    Plot calibrated minispect illuminance on an existing axis.
    """

    ax.plot(as_time_relative, illuminance, label="Illuminance")
    ax.set_xlabel("Time from minispect recording start (s)")
    ax.set_ylabel("Illuminance (lux)")
    ax.set_title("Calibrated minispect illuminance")
    ax.grid()
    _set_legend_top_right(ax, fontsize=8)
    return ax


def fit_illuminance_on_camera_score(
    camera_score: np.ndarray,
    illuminance: np.ndarray,
) -> dict[str, Any]:
    """
    Robustly fit illuminance against the camera AGC product using Huber loss.

    The affine fit is oriented as
    illuminance = intercept + slope * camera_score. Because camera score is
    the product of analog gain, digital gain, and exposure, its slope may be
    negative. Huber loss retains quadratic weighting for ordinary residuals
    while reducing the influence of unusually large residuals.

    Camera-score standard deviation is the sample standard deviation of the
    valid paired scores. It is not divided by the square root of the sample
    count and therefore is not a standard error.
    """

    camera_score = np.asarray(camera_score, dtype=float)
    illuminance = np.asarray(illuminance, dtype=float)
    valid: np.ndarray = np.isfinite(camera_score) & np.isfinite(illuminance)
    valid_camera_score: np.ndarray = camera_score[valid]
    valid_illuminance: np.ndarray = illuminance[valid]
    n_samples: int = int(np.sum(valid))
    camera_score_std: float = (
        float(np.std(valid_camera_score, ddof=1)) if n_samples >= 2 else np.nan
    )

    empty_fit: dict[str, Any] = {
        "slope": np.nan,
        "intercept": np.nan,
        "correlation": np.nan,
        "camera_score_std": camera_score_std,
        "n_samples": n_samples,
        "fit_method": "Huber affine regression",
        "residual_scale": np.nan,
    }

    # Retain low-variance recordings so their potentially unstable slopes can
    # be diagnosed using camera_score_std in the returned summary.
    if n_samples < 2:
        return empty_fit

    # Center and scale the AGC product before optimization so its much larger
    # exposure-driven values remain numerically well-conditioned.
    score_center: float = float(np.median(valid_camera_score))
    centered_camera_score: np.ndarray = valid_camera_score - score_center
    score_scale: float = float(np.max(np.abs(centered_camera_score)))
    if score_scale == 0:
        return empty_fit

    scaled_camera_score: np.ndarray = centered_camera_score / score_scale
    initial_coefficients: np.ndarray = np.polyfit(
        scaled_camera_score,
        valid_illuminance,
        deg=1,
    )
    initial_scaled_slope: float = float(initial_coefficients[0])
    initial_scaled_intercept: float = float(initial_coefficients[1])

    # Estimate the ordinary residual scale robustly. SciPy's f_scale sets the
    # transition between squared-error and linearly growing Huber penalties.
    initial_residuals: np.ndarray = (
        initial_scaled_intercept
        + initial_scaled_slope * scaled_camera_score
        - valid_illuminance
    )
    residual_median: float = float(np.median(initial_residuals))
    median_absolute_deviation: float = float(
        np.median(np.abs(initial_residuals - residual_median))
    )
    residual_scale: float = 1.4826 * median_absolute_deviation

    # When most residuals are identical, MAD can be zero even if a few extreme
    # outliers remain. Use a numerical floor based on the illuminance scale so
    # those outliers do not define their own rejection threshold.
    if residual_scale == 0:
        illuminance_scale: float = max(
            float(np.max(np.abs(valid_illuminance))),
            1.0,
        )
        residual_scale = np.sqrt(np.finfo(float).eps) * illuminance_scale

    def _affine_residuals(parameters: np.ndarray) -> np.ndarray:
        """Return affine illuminance residuals for slope and intercept."""

        scaled_slope: float = float(parameters[0])
        scaled_intercept: float = float(parameters[1])
        return (
            scaled_intercept
            + scaled_slope * scaled_camera_score
            - valid_illuminance
        )

    # Leave the slope unconstrained because the AGC product is expected to
    # decrease as environmental illuminance increases.
    fit_result: OptimizeResult = least_squares(
        _affine_residuals,
        x0=np.array(
            [initial_scaled_slope, initial_scaled_intercept],
            dtype=float,
        ),
        loss="huber",
        f_scale=residual_scale,
        x_scale="jac",
    )
    robust_scaled_slope: float = float(fit_result.x[0])
    robust_scaled_intercept: float = float(fit_result.x[1])

    slope: float = robust_scaled_slope / score_scale
    intercept: float = (
        robust_scaled_intercept - slope * score_center
    )

    correlation: np.float64 = np.nan
    if np.ptp(valid_camera_score) > 0 and np.ptp(valid_illuminance) > 0:
        correlation = np.corrcoef(valid_camera_score, valid_illuminance)[0, 1]

    return {
        "slope": slope,
        "intercept": intercept,
        "correlation": float(correlation),
        "camera_score_std": camera_score_std,
        "n_samples": n_samples,
        "fit_method": "Huber affine regression",
        "residual_scale": residual_scale,
    }


def _plot_absolute_scatter(
    ax: Axes,
    scatter_camera_score: np.ndarray,
    scatter_illuminance: np.ndarray,
    title: str = "Camera AGC product vs filtered illuminance",
) -> Axes:
    """
    Plot paired camera AGC products against filtered illuminance values.

    A robust Huber line is fit through the paired samples in the plotted
    orientation, so the slope is reported in lux per unit camera score.
    """

    ax.scatter(
        scatter_camera_score,
        scatter_illuminance,
        s=12,
        alpha=0.5,
        label="Paired samples",
    )

    camera_score: np.ndarray = np.asarray(scatter_camera_score, dtype=float)
    illuminance: np.ndarray = np.asarray(scatter_illuminance, dtype=float)
    line_fit: dict[str, Any] = fit_illuminance_on_camera_score(
        camera_score,
        illuminance,
    )
    slope: float = line_fit["slope"]
    intercept: float = line_fit["intercept"]

    if np.isfinite(slope):
        valid: np.ndarray = np.isfinite(camera_score) & np.isfinite(illuminance)

        # Draw the affine fit across the observed AGC-product range.
        score_fit: np.ndarray = np.linspace(
            float(np.min(camera_score[valid])),
            float(np.max(camera_score[valid])),
            200,
        )
        ax.plot(
            score_fit,
            intercept + slope * score_fit,
            color="black",
            linewidth=2,
            label=(
                "Robust fit: illuminance = "
                f"{intercept:.4g} + {slope:.4g} x score"
            ),
        )

    ax.set_xlabel("Camera AGC product")
    ax.set_ylabel("Filtered illuminance (lux)")
    ax.set_title(title)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid()
    _set_legend_top_right(ax, fontsize=8)
    return ax


def _plot_calibration_fit(
    ax: Axes,
    calibration_result: dict[str, Any],
) -> Axes:
    """
    Plot pooled sparse calibration points and their linear fit.
    """

    filtered_illuminance: np.ndarray = calibration_result["filtered_illuminance"]
    camera_score: np.ndarray = calibration_result["camera_score"]
    slope: float = calibration_result["slope"]
    intercept: float = calibration_result["intercept"]

    ax.scatter(
        camera_score,
        filtered_illuminance,
        s=12,
        alpha=0.45,
        label="Sparse paired samples",
    )

    valid: np.ndarray = np.isfinite(filtered_illuminance) & np.isfinite(camera_score)
    if np.any(valid):
        # Draw the fitted relationship in the same orientation as both the
        # per-recording scatter and the implied-illuminance equation.
        camera_score_fit: np.ndarray = np.linspace(
            np.nanmin(camera_score[valid]),
            np.nanmax(camera_score[valid]),
            200,
        )
        illuminance_fit: np.ndarray = intercept + slope * camera_score_fit
        ax.plot(
            camera_score_fit,
            illuminance_fit,
            color="black",
            linewidth=2,
            label=(
                "Robust fit: illuminance = "
                f"{intercept:.4g} + {slope:.4g} x score"
            ),
        )

    ax.set_xlabel("Camera AGC product")
    ax.set_ylabel("Filtered illuminance (lux)")
    ax.set_title("Robust steady-state AGC calibration")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid()
    _set_legend_top_right(ax, fontsize=8)
    return ax


def _plot_temporal_filter_inputs(
    ax: Axes,
    model_result: dict[str, Any],
) -> Axes:
    """
    Plot normalized camera and illuminance traces before kernel filtering.
    """

    relative_ms_time: np.ndarray = (
        model_result["time"] - model_result["time"][0]
    )
    ax.plot(
        relative_ms_time,
        model_result["illuminance"],
        label="Illuminance (z-scored)",
    )
    ax.plot(
        relative_ms_time,
        model_result["camera"],
        label="Brightness-oriented camera AGC product",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized value")
    ax.set_title("Normalized camera and illuminance")
    _set_legend_top_right(ax, fontsize=8)
    ax.grid()
    return ax


def _plot_temporal_filter_prediction(
    ax: Axes,
    model_result: dict[str, Any],
) -> Axes:
    """
    Plot observed camera score against the shared empirical-kernel prediction.
    """

    relative_ms_time: np.ndarray = (
        model_result["time"] - model_result["time"][0]
    )
    ax.plot(
        relative_ms_time,
        model_result["camera"],
        label="Observed brightness-oriented AGC product",
    )
    ax.plot(
        relative_ms_time,
        model_result["prediction"],
        label=f"Empirical kernel (lag={model_result['best_lag']:.2f} s)",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized value")
    ax.set_title("Empirical-kernel prediction")
    _set_legend_top_right(ax, fontsize=8)
    ax.grid()
    return ax


def _plot_filtered_illuminance_prediction(
    ax: Axes,
    model_result: dict[str, Any],
) -> Axes:
    """
    Plot the empirical-kernel output in physical illuminance units.
    """

    ax.plot(
        model_result["time"] - model_result["time"][0],
        model_result["filtered_illuminance"],
        label="Filtered illuminance",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Illuminance (lux)")
    ax.set_title("Empirical filtered illuminance")
    _set_legend_top_right(ax, fontsize=8)
    ax.grid()
    return ax


def _plot_filter_lag_search(
    ax: Axes,
    model_result: dict[str, Any],
) -> Axes:
    """
    Plot one recording's empirical-kernel correlation versus lag.
    """

    ax.plot(model_result["candidate_lags"], model_result["correlations"])
    ax.axvline(
        model_result["best_lag"],
        linestyle="--",
        label=f"Shared lag = {model_result['best_lag']:.2f} s",
    )
    ax.set_xlabel("Camera lag relative to minispect (s)")
    ax.set_ylabel("Correlation")
    ax.set_title("Empirical-kernel lag search")
    _set_legend_top_right(ax, fontsize=8)
    ax.grid()
    return ax


def _plot_shared_lag_summary(
    ax: Axes,
    shared_model_result: dict[str, Any],
) -> Axes:
    """
    Plot the aggregate correlation used to choose the shared model lag.
    """

    ax.plot(
        shared_model_result["candidate_lags"],
        shared_model_result["mean_correlation_by_lag"],
    )
    ax.axvline(
        shared_model_result["shared_lag"],
        linestyle="--",
        label=f"Shared lag = {shared_model_result['shared_lag']:.2f} s",
    )
    ax.set_xlabel("Camera lag relative to minispect (s)")
    ax.set_ylabel("Mean correlation")
    ax.set_title("Shared empirical-model lag")
    _set_legend_top_right(ax, fontsize=8)
    ax.grid()
    return ax


def _recording_identity_from_path(recording_path: str) -> tuple[str, str]:
    """
    Split a recording path into its subject ID and activity name.
    """

    normalized_path: str = recording_path.rstrip("/")
    terminal_name: str = os.path.basename(normalized_path)
    recording_name: str = os.path.basename(os.path.dirname(normalized_path))
    participant_name: str = os.path.basename(os.path.dirname(os.path.dirname(normalized_path)))

    if terminal_name != "GKA":
        raise ValueError(
            "Recording paths must have format "
            ".../FLIC_<NUMERICALID>/<RECORDING_NAME>/GKA."
        )

    if not participant_name.startswith("FLIC_"):
        raise ValueError(
            "Recording paths must have format "
            ".../FLIC_<NUMERICALID>/<RECORDING_NAME>/GKA."
        )

    participant_id: str = participant_name.removeprefix("FLIC_")
    if not participant_id.isdigit() or len(recording_name) == 0:
        raise ValueError(
            "Recording paths must have format "
            ".../FLIC_<NUMERICALID>/<RECORDING_NAME>/GKA."
        )

    return participant_id, recording_name


def _recording_name_from_path(recording_path: str) -> str:
    """
    Return a compact recording label as NUMERICALID_RECORDING_NAME.
    """

    subject_id: str
    activity: str
    subject_id, activity = _recording_identity_from_path(recording_path)

    return f"{subject_id}_{activity}"


def _processed_world_video_path(recording_path: str) -> str:
    """Construct the processed W.avi path corresponding to a raw recording."""

    # Inputs point to raw GKA chunk directories. The processed world video
    # mirrors that directory tree under FLIC_processing and is named W.avi.
    if "FLIC_raw" not in recording_path:
        raise ValueError(
            "Recording paths must point into FLIC_raw so the corresponding "
            "FLIC_processing video path can be constructed."
        )

    processed_recording_path: str = recording_path.replace(
        "FLIC_raw",
        "FLIC_processing",
    )
    return os.path.join(processed_recording_path, "W.avi")


def _padded_limits(
    value_arrays: Sequence[np.ndarray],
    padding_fraction: float = 0.05,
) -> tuple[float, float] | None:
    """
    Return pooled min/max of several arrays, padded by a fraction of range.
    """

    finite_values: list[np.ndarray] = []
    for value_array in value_arrays:
        flat_values: np.ndarray = np.asarray(value_array, dtype=float).ravel()
        finite_values.append(flat_values[np.isfinite(flat_values)])

    pooled_values: np.ndarray = (
        np.concatenate(finite_values) if len(finite_values) > 0
        else np.array([], dtype=float)
    )

    if pooled_values.size == 0:
        return None

    minimum_value: float = float(np.min(pooled_values))
    maximum_value: float = float(np.max(pooled_values))
    value_range: float = maximum_value - minimum_value

    # A flat signal has no range to scale, so pad relative to its magnitude
    # instead and fall back to a fixed margin for an all-zero signal.
    if value_range == 0:
        padding: float = (
            padding_fraction * abs(maximum_value) if maximum_value != 0 else 1.0
        )
    else:
        padding = padding_fraction * value_range

    return (minimum_value - padding, maximum_value + padding)


def _relative_time(time_values: np.ndarray) -> np.ndarray:
    """
    Return a time vector shifted so that it starts at zero.
    """

    time_values = np.asarray(time_values, dtype=float)

    if time_values.size == 0:
        return time_values

    return time_values - time_values[0]


def compute_shared_axis_limits(
    results_by_recording: Sequence[dict[str, Any]],
    padding_fraction: float = 0.05,
) -> dict[str, dict[str, tuple[float, float] | None]]:
    """
    Compute per-panel x and y axis limits pooled across every recording.

    Panels are grouped by the quantity they display rather than pooled into
    one global range, because the dashboard mixes seconds, lux, camera light
    score, z-scores, and correlations across its panels. The three normalized
    panels share one range so their traces stay directly comparable, and the
    time axes span the longest recording so shorter ones read as shorter.
    """

    camera_score_arrays: list[np.ndarray] = []
    as_value_arrays: list[np.ndarray] = []
    illuminance_arrays: list[np.ndarray] = []
    scatter_illuminance_arrays: list[np.ndarray] = []
    normalized_arrays: list[np.ndarray] = []
    correlation_arrays: list[np.ndarray] = []

    camera_time_arrays: list[np.ndarray] = []
    as_time_arrays: list[np.ndarray] = []
    saturation_time_arrays: list[np.ndarray] = []
    scatter_camera_score_arrays: list[np.ndarray] = []
    filter_time_arrays: list[np.ndarray] = []
    lag_arrays: list[np.ndarray] = []

    for recording_result in results_by_recording:
        camera_score_arrays.append(recording_result["camera_score"])
        as_value_arrays.append(recording_result["as_values"])
        illuminance_arrays.append(recording_result["illuminance"])
        scatter_illuminance_arrays.append(recording_result["scatter_illuminance"])

        camera_time_arrays.append(recording_result["camera_time_relative"])
        as_time_arrays.append(recording_result["as_time_relative"])
        saturation_time_arrays.append(
            recording_result["saturation_time_relative"]
        )
        scatter_camera_score_arrays.append(recording_result["scatter_camera_score"])

        model_result: dict[str, Any] = recording_result["model_result"]
        normalized_arrays.extend(
            (
                model_result["illuminance"],
                model_result["camera"],
                model_result["prediction"],
            )
        )

        # The filter panels plot time relative to their own first sample, so
        # pool the same relative vectors that those panels actually draw.
        filter_time_arrays.append(_relative_time(model_result["time"]))
        lag_arrays.append(model_result["candidate_lags"])
        correlation_arrays.append(model_result["correlations"])

    return {
        "y": {
            "camera_light_score": _padded_limits(
                camera_score_arrays,
                padding_fraction,
            ),
            "as_channels": _padded_limits(as_value_arrays, padding_fraction),
            "illuminance": _padded_limits(illuminance_arrays, padding_fraction),
            "scatter_illuminance": _padded_limits(
                scatter_illuminance_arrays,
                padding_fraction,
            ),
            "normalized": _padded_limits(normalized_arrays, padding_fraction),
            "correlation": _padded_limits(correlation_arrays, padding_fraction),
        },
        "x": {
            "camera_time": _padded_limits(camera_time_arrays, padding_fraction),
            "as_time": _padded_limits(as_time_arrays, padding_fraction),
            "saturation_time": _padded_limits(
                saturation_time_arrays,
                padding_fraction,
            ),
            "scatter_camera_score": _padded_limits(
                scatter_camera_score_arrays,
                padding_fraction,
            ),
            "filter_time": _padded_limits(filter_time_arrays, padding_fraction),
            "lag": _padded_limits(lag_arrays, padding_fraction),
        },
    }


def build_recording_summary_table(
    results_by_recording: Sequence[dict[str, Any]],
) -> pd.DataFrame:
    """
    Build a per-recording summary of model correlation and calibration fit.

    One row per recording is identified by subject ID and activity. Every
    metric comes from the shared empirical-kernel model: correlation measures
    its temporal prediction, while slope and intercept describe the robust
    affine physical-unit relationship

        illuminance = intercept + slope * camera_score.

    Camera-score standard deviation is calculated from the paired samples
    used by that fit so low-variance recordings can be identified directly.

    Recordings that failed illuminance conversion still get a row, with NaN
    metrics and the failure reason, so an excluded recording is visible in
    the summary rather than silently absent.
    """

    summary_rows: list[dict[str, Any]] = []

    for recording_result in results_by_recording:
        subject_id: str
        activity: str
        subject_id, activity = _recording_identity_from_path(
            recording_result["recording_path"]
        )

        if recording_result.get("errored", False):
            summary_rows.append(
                {
                    "subject_id": subject_id,
                    "activity": activity,
                    "errored": True,
                    "correlation": np.nan,
                    "shared_lag_s": np.nan,
                    "slope": np.nan,
                    "intercept": np.nan,
                    "camera_score_std": np.nan,
                    "n_paired_samples": 0,
                    "error_message": recording_result.get("error_message", ""),
                    "recording_path": recording_result["recording_path"],
                }
            )
            continue

        model_result: dict[str, Any] = recording_result["model_result"]
        line_fit: dict[str, Any] = recording_result["absolute_fit"]

        summary_rows.append(
            {
                "subject_id": subject_id,
                "activity": activity,
                "errored": False,
                "correlation": float(model_result["best_correlation"]),
                "shared_lag_s": float(model_result["best_lag"]),
                "slope": line_fit["slope"],
                "intercept": line_fit["intercept"],
                "camera_score_std": line_fit["camera_score_std"],
                "n_paired_samples": line_fit["n_samples"],
                "error_message": "",
                "recording_path": recording_result["recording_path"],
            }
        )

    # Naming the columns explicitly keeps the schema stable even when no
    # recording produced a row.
    summary_table: pd.DataFrame = pd.DataFrame(
        summary_rows,
        columns=list(SUMMARY_TABLE_COLUMNS),
    )

    if len(summary_table) > 0:
        # Best-fitting recordings first. Errored rows carry a NaN correlation
        # and sort to the bottom rather than disappearing.
        summary_table = summary_table.sort_values(
            ["correlation", "subject_id", "activity"],
            ascending=[False, True, True],
            na_position="last",
            ignore_index=True,
        )

    return summary_table


def _plot_empirical_kernel_summary(
    results_by_recording: Sequence[dict[str, Any]],
    summary_table: pd.DataFrame,
) -> np.ndarray:
    """
    Plot compact empirical-kernel summaries across every recording.
    """

    figure, axes = plt.subplots(1, 3, figsize=(21, 7))
    axes_flat: np.ndarray = axes.ravel()
    figure.suptitle("Empirical Kernel Model Fit - All Recordings", fontsize=16)

    # SCATTER AND PER-RECORDING CALIBRATION LINES
    # -------------------------------------------------------------
    for recording_result in results_by_recording:
        label: str = _recording_name_from_path(recording_result["recording_path"])
        line_fit: dict[str, Any] = recording_result["absolute_fit"]

        scatter_camera_score: np.ndarray = recording_result[
            "scatter_camera_score"
        ]
        scatter_illuminance: np.ndarray = recording_result[
            "scatter_illuminance"
        ]

        points: Any = axes_flat[0].scatter(
            scatter_camera_score,
            scatter_illuminance,
            s=8,
            alpha=0.35,
            label=label,
        )

        if np.isfinite(line_fit["slope"]) and scatter_camera_score.size > 0:
            score_fit: np.ndarray = np.linspace(
                float(np.min(scatter_camera_score)),
                float(np.max(scatter_camera_score)),
                200,
            )
            axes_flat[0].plot(
                score_fit,
                line_fit["intercept"] + line_fit["slope"] * score_fit,
                color=points.get_facecolor()[0],
                linewidth=2,
            )

    axes_flat[0].set_xlabel("Camera AGC product")
    axes_flat[0].set_ylabel("Empirical-kernel illuminance (lux)")
    axes_flat[0].set_title("Robust AGC product vs empirical-kernel illuminance")
    axes_flat[0].set_xlim(left=0)
    axes_flat[0].set_ylim(bottom=0)

    # MODEL CORRELATION BY RECORDING
    # -------------------------------------------------------------
    recording_labels: list[str] = [
        f"{row.subject_id}_{row.activity}" for row in summary_table.itertuples()
    ]
    label_positions: np.ndarray = np.arange(len(recording_labels))

    axes_flat[1].bar(
        label_positions,
        summary_table["correlation"],
        label="Empirical kernel",
    )
    axes_flat[1].set_xticks(label_positions)
    axes_flat[1].set_xticklabels(recording_labels, rotation=45, ha="right", fontsize=8)
    axes_flat[1].set_ylabel("Brightness-oriented AGC-product correlation")
    axes_flat[1].set_title("Correlation by recording")

    # CALIBRATION SLOPE COMPARISON
    # -------------------------------------------------------------
    axes_flat[2].bar(
        label_positions,
        summary_table["slope"],
        label="Empirical kernel",
    )
    axes_flat[2].set_xticks(label_positions)
    axes_flat[2].set_xticklabels(recording_labels, rotation=45, ha="right", fontsize=8)
    axes_flat[2].set_ylabel("Slope (lux per AGC-product unit)")
    axes_flat[2].set_title("Calibration slope by recording")

    # The per-recording numbers are returned as a DataFrame rather than
    # rendered here, so this figure stays plots only.
    for ax in axes_flat:
        ax.grid()
        _set_legend_top_right(ax, fontsize=7)

    figure.tight_layout(rect=(0, 0, 1, 0.93))
    _draw_figure_now(figure)
    return axes


def _plot_recording_diagnostics(
    recording_result: dict[str, Any],
    axis_limits: dict[str, dict[str, tuple[float, float] | None]] | None = None,
) -> np.ndarray:
    """
    Plot all per-recording diagnostics in one subplot dashboard figure.

    The first absolute scatter uses limits pooled across recordings for direct
    comparison, while the adjacent scatter keeps video-local limits so its
    within-recording structure remains visible. Camera settings keeps its own
    fixed gain and exposure scaling.
    """

    figure, axes = plt.subplots(3, 4, figsize=(22, 14))
    axes_flat: np.ndarray = axes.ravel()
    recording_name: str = _recording_name_from_path(recording_result["recording_path"])
    figure.suptitle(recording_name, fontsize=16)

    _plot_camera_settings(
        axes_flat[0],
        recording_result["camera_time_relative"],
        recording_result["analog_gain"],
        recording_result["digital_gain"],
        recording_result["exposure"],
    )
    _plot_camera_light_score(
        axes_flat[1],
        recording_result["camera_time_relative"],
        recording_result["camera_score"],
    )
    _plot_raw_as_channels(
        axes_flat[2],
        recording_result["as_time_relative"],
        recording_result["as_values"],
    )
    _plot_illuminance(
        axes_flat[3],
        recording_result["as_time_relative"],
        recording_result["illuminance"],
    )
    _plot_absolute_scatter(
        axes_flat[4],
        recording_result["scatter_camera_score"],
        recording_result["scatter_illuminance"],
        title="AGC product vs filtered illuminance (global limits)",
    )
    _plot_absolute_scatter(
        axes_flat[5],
        recording_result["scatter_camera_score"],
        recording_result["scatter_illuminance"],
        title="AGC product vs filtered illuminance (video-local limits)",
    )
    _plot_temporal_filter_inputs(axes_flat[6], recording_result["model_result"])
    _plot_temporal_filter_prediction(axes_flat[7], recording_result["model_result"])
    _plot_filtered_illuminance_prediction(
        axes_flat[8],
        recording_result["model_result"],
    )
    _plot_filter_lag_search(axes_flat[9], recording_result["model_result"])
    _plot_video_saturation(
        axes_flat[10],
        recording_result["saturation_time_relative"],
        recording_result["video_saturation"],
    )
    axes_flat[11].axis("off")

    # Apply the pooled limits after the panels are built so each helper stays
    # usable on its own. Index 0 is camera settings, which is left alone
    # because its twin axes already use fixed gain and exposure scaling.
    if axis_limits is not None:
        y_limit_keys: dict[int, str] = {
            1: "camera_light_score",
            2: "as_channels",
            3: "illuminance",
            4: "scatter_illuminance",
            6: "normalized",
            7: "normalized",
            8: "scatter_illuminance",
            9: "correlation",
        }
        x_limit_keys: dict[int, str] = {
            1: "camera_time",
            2: "as_time",
            3: "as_time",
            4: "scatter_camera_score",
            6: "filter_time",
            7: "filter_time",
            8: "filter_time",
            9: "lag",
            10: "saturation_time",
        }

        y_limits: dict[str, tuple[float, float] | None] = axis_limits.get("y", {})
        x_limits: dict[str, tuple[float, float] | None] = axis_limits.get("x", {})

        for panel_index, limit_key in y_limit_keys.items():
            panel_limits: tuple[float, float] | None = y_limits.get(limit_key)
            if panel_limits is not None:
                axes_flat[panel_index].set_ylim(*panel_limits)

        for panel_index, limit_key in x_limit_keys.items():
            panel_limits = x_limits.get(limit_key)
            if panel_limits is not None:
                axes_flat[panel_index].set_xlim(*panel_limits)

    # Camera AGC products and illuminance are nonnegative measurements. Their
    # upper limits remain global or local as titled; the affine fit itself is
    # drawn only across the observed score range.
    for panel_index in (4, 5):
        axes_flat[panel_index].set_xlim(left=0)
        axes_flat[panel_index].set_ylim(bottom=0)

    figure.tight_layout(rect=(0, 0, 1, 0.96))
    _draw_figure_now(figure)
    return axes


def _gain_log_variation(gain: np.ndarray) -> float:
    """Return sample variation in a gain's multiplicative, logarithmic scale."""

    gain = np.asarray(gain, dtype=float)
    valid_gain: np.ndarray = gain[np.isfinite(gain) & (gain > 0)]

    # Two observations are required for a sample standard deviation. A fixed
    # or undersampled gain contributes no evidence of active variation.
    if valid_gain.size < 2:
        return 0.0

    return float(np.std(np.log(valid_gain), ddof=1))


def _primary_gain_control(
    analog_gain: np.ndarray,
    digital_gain: np.ndarray,
) -> str:
    """Classify a recording by which gain control varies more proportionally."""

    # Camera gain enters the light score multiplicatively, so logarithmic
    # variation compares proportional analog and digital changes in one scale.
    analog_variation: float = _gain_log_variation(analog_gain)
    digital_variation: float = _gain_log_variation(digital_gain)

    return "analog" if analog_variation >= digital_variation else "digital"


def _distinct_activity_colors(
    number_of_colors: int,
) -> list[tuple[float, float, float]]:
    """Return maximally separated, high-contrast colors for activity points."""

    if number_of_colors <= 0:
        return []

    # Build a broad set of saturated colors while excluding very dark, very
    # light, and gray candidates that are difficult to see on a white plot.
    channel_levels: np.ndarray = np.linspace(0.05, 0.95, 7)
    luminance_weights: np.ndarray = np.array([0.2126, 0.7152, 0.0722])
    candidate_colors: np.ndarray = np.array(
        [
            color
            for color in itertools.product(channel_levels, repeat=3)
            if max(color) - min(color) >= 0.45
            and 0.12 <= np.dot(color, luminance_weights) <= 0.82
        ],
        dtype=float,
    )

    if number_of_colors > candidate_colors.shape[0]:
        raise ValueError(
            f"Requested {number_of_colors} activity colors, but only "
            f"{candidate_colors.shape[0]} distinguishable candidates exist."
        )

    # Begin with a vivid blue, then repeatedly choose the candidate farthest
    # from its nearest selected color. This directly avoids adjacent shades.
    initial_color: np.ndarray = np.array([0.05, 0.35, 0.95])
    initial_distances: np.ndarray = np.sum(
        (candidate_colors - initial_color) ** 2,
        axis=1,
    )
    selected_indices: list[int] = [int(np.argmin(initial_distances))]
    available: np.ndarray = np.ones(candidate_colors.shape[0], dtype=bool)
    available[selected_indices[0]] = False

    while len(selected_indices) < number_of_colors:
        distance_squared: np.ndarray = np.sum(
            (
                candidate_colors[:, np.newaxis, :]
                - candidate_colors[selected_indices][np.newaxis, :, :]
            )
            ** 2,
            axis=2,
        )
        nearest_selected_distance: np.ndarray = np.min(
            distance_squared,
            axis=1,
        )
        nearest_selected_distance[~available] = -1.0
        next_index: int = int(np.argmax(nearest_selected_distance))
        selected_indices.append(next_index)
        available[next_index] = False

    selected_colors: np.ndarray = candidate_colors[selected_indices]
    return [
        tuple(float(channel) for channel in color)
        for color in selected_colors
    ]


def _plot_high_correlation_calibration(
    results_by_recording: Sequence[dict[str, Any]],
    minimum_correlation: float = 0.9,
) -> np.ndarray:
    """
    Plot two correlation-filtered views of the calibration data.

    The left panel gives each activity a distinct color. In the right panel,
    recordings driven primarily by analog-gain variation are red, while those
    driven primarily by digital-gain variation are black. Both panels display
    the same robust pooled affine line.
    """

    figure: Any
    axes: np.ndarray
    figure, axes = plt.subplots(1, 2, figsize=(20, 8))
    axes_flat: np.ndarray = axes.ravel()

    # FILTER RECORDINGS BY SHARED-MODEL CORRELATION
    # -------------------------------------------------------------
    selected_recordings: list[dict[str, Any]] = [
        recording_result
        for recording_result in results_by_recording
        if np.isfinite(recording_result["model_result"]["best_correlation"])
        and recording_result["model_result"]["best_correlation"]
        >= minimum_correlation
    ]
    pooled_camera_scores: list[np.ndarray] = []
    pooled_illuminance: list[np.ndarray] = []
    activity_colors: list[tuple[float, float, float]] = _distinct_activity_colors(
        len(selected_recordings)
    )
    activity_markers: tuple[str, ...] = (
        "o",
        "s",
        "^",
        "D",
        "v",
        "P",
        "X",
        "<",
        ">",
        "*",
    )

    # PLOT CORRELATION-FILTERED ACTIVITIES BY NAME AND PRIMARY GAIN CONTROL
    # -------------------------------------------------------------
    for recording_index, recording_result in enumerate(selected_recordings):
        primary_gain: str = _primary_gain_control(
            recording_result["analog_gain"],
            recording_result["digital_gain"],
        )
        agc_type_color: str = "red" if primary_gain == "analog" else "black"
        activity_color: tuple[float, float, float] = activity_colors[recording_index]
        activity_marker: str = activity_markers[
            recording_index % len(activity_markers)
        ]
        recording_label: str = _recording_name_from_path(
            recording_result["recording_path"]
        )
        camera_score: np.ndarray = np.asarray(
            recording_result["scatter_camera_score"],
            dtype=float,
        )
        filtered_illuminance: np.ndarray = np.asarray(
            recording_result["scatter_illuminance"],
            dtype=float,
        )
        valid: np.ndarray = (
            np.isfinite(camera_score) & np.isfinite(filtered_illuminance)
        )

        axes_flat[0].scatter(
            camera_score[valid],
            filtered_illuminance[valid],
            color=activity_color,
            marker=activity_marker,
            s=18,
            alpha=0.8,
            edgecolors="white",
            linewidths=0.25,
            label=recording_label,
        )
        axes_flat[1].scatter(
            camera_score[valid],
            filtered_illuminance[valid],
            color=agc_type_color,
            s=10,
            alpha=0.35,
            label=recording_label,
        )
        pooled_camera_scores.append(camera_score[valid])
        pooled_illuminance.append(filtered_illuminance[valid])

    # FIT ONE ROBUST LINE ACROSS THE DISPLAYED ACTIVITIES
    # -------------------------------------------------------------
    nonempty_scores: list[np.ndarray] = [
        score for score in pooled_camera_scores if score.size > 0
    ]
    nonempty_illuminance: list[np.ndarray] = [
        illuminance for illuminance in pooled_illuminance if illuminance.size > 0
    ]

    if len(nonempty_scores) > 0:
        combined_camera_score: np.ndarray = np.concatenate(nonempty_scores)
        combined_illuminance: np.ndarray = np.concatenate(nonempty_illuminance)
        line_fit: dict[str, Any] = fit_illuminance_on_camera_score(
            combined_camera_score,
            combined_illuminance,
        )
        slope: float = line_fit["slope"]
        intercept: float = line_fit["intercept"]

        if np.isfinite(slope):
            score_fit: np.ndarray = np.linspace(
                float(np.min(combined_camera_score)),
                float(np.max(combined_camera_score)),
                200,
            )
            for ax in axes_flat:
                ax.plot(
                    score_fit,
                    intercept + slope * score_fit,
                    color="tab:blue",
                    linewidth=2.5,
                    label=(
                        "Shared calibration fit: "
                        f"illuminance = {intercept:.4g} "
                        f"+ {slope:.4g} x score"
                    ),
                )

    # FORMAT THE TWO MATCHED CALIBRATION PANELS
    # -------------------------------------------------------------
    figure.suptitle(
        "Camera AGC product vs filtered illuminance "
        f"(shared-model correlation >= {minimum_correlation:.2f})",
        fontsize=16,
    )
    axes_flat[0].set_title("Activities by Name")
    axes_flat[1].set_title("Activities by AGC Type")

    for ax in axes_flat:
        ax.set_xlabel("Camera AGC product")
        ax.set_ylabel("Filtered illuminance (lux)")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid()
    if len(selected_recordings) > 0:
        _set_legend_top_right(axes_flat[0], fontsize=8)
        _set_legend_top_right(axes_flat[1], fontsize=8)

    figure.tight_layout(rect=(0, 0, 1, 0.95))
    _draw_figure_now(figure)
    return axes


def _fit_shifted_reciprocal_illuminance(
    camera_agc_product: np.ndarray,
    illuminance: np.ndarray,
) -> dict[str, float | int]:
    """
    Robustly fit an affine shifted reciprocal model using Huber loss.

    The model is ``illuminance = intercept + numerator / (product + offset)``.
    As in the original reciprocal fit, residuals are measured in untransformed
    lux and the intercept is unconstrained. The numerator and denominator
    offset are constrained to be nonnegative.
    """

    camera_agc_product = np.asarray(camera_agc_product, dtype=float)
    illuminance = np.asarray(illuminance, dtype=float)
    if camera_agc_product.shape != illuminance.shape:
        raise ValueError(
            "Camera AGC product and illuminance must have matching shapes."
        )

    # Log-scale plotting and reciprocal modeling both require positive paired
    # observations, so remove unusable values before estimating parameters.
    valid: np.ndarray = (
        np.isfinite(camera_agc_product)
        & (camera_agc_product > 0)
        & np.isfinite(illuminance)
        & (illuminance > 0)
    )
    valid_product: np.ndarray = camera_agc_product[valid]
    valid_illuminance: np.ndarray = illuminance[valid]
    n_samples: int = int(valid_product.size)
    empty_fit: dict[str, float | int] = {
        "intercept": np.nan,
        "numerator": np.nan,
        "denominator_offset": np.nan,
        "r_squared": np.nan,
        "residual_scale": np.nan,
        "n_samples": n_samples,
    }
    if n_samples < 2 or np.ptp(valid_product) == 0:
        return empty_fit

    # Scale the large AGC products before nonlinear optimization. In scaled
    # coordinates the model is intercept + beta / (product + gamma).
    product_scale: float = float(np.median(valid_product))
    scaled_product: np.ndarray = valid_product / product_scale
    initial_coefficients: np.ndarray = np.polyfit(
        1.0 / scaled_product,
        valid_illuminance,
        deg=1,
    )
    initial_beta: float = max(float(initial_coefficients[0]), 0.0)
    initial_intercept: float = float(initial_coefficients[1])
    initial_gamma: float = 0.0
    initial_prediction: np.ndarray = (
        initial_intercept + initial_beta / scaled_product
    )

    # Match the original fit by estimating the Huber transition from residuals
    # in lux, before adding the denominator offset to the model.
    initial_residual: np.ndarray = initial_prediction - valid_illuminance
    residual_median: float = float(np.median(initial_residual))
    median_absolute_deviation: float = float(
        np.median(np.abs(initial_residual - residual_median))
    )
    residual_scale: float = 1.4826 * median_absolute_deviation
    if residual_scale == 0:
        residual_scale = np.sqrt(np.finfo(float).eps) * max(
            float(np.max(np.abs(valid_illuminance))),
            1.0,
        )

    def _shifted_reciprocal_residuals(
        parameters: np.ndarray,
    ) -> np.ndarray:
        """Return raw-lux residuals for the shifted reciprocal model."""

        intercept: float = float(parameters[0])
        beta: float = float(parameters[1])
        gamma: float = float(parameters[2])
        return (
            intercept
            + beta / (scaled_product + gamma)
            - valid_illuminance
        )

    fit_result: OptimizeResult = least_squares(
        _shifted_reciprocal_residuals,
        x0=np.array(
            [initial_intercept, initial_beta, initial_gamma],
            dtype=float,
        ),
        bounds=(
            np.array([-np.inf, 0.0, 0.0], dtype=float),
            np.array([np.inf, np.inf, np.inf], dtype=float),
        ),
        loss="huber",
        f_scale=residual_scale,
        x_scale="jac",
    )
    fitted_intercept: float = float(fit_result.x[0])
    fitted_beta: float = float(fit_result.x[1])
    fitted_gamma: float = float(fit_result.x[2])
    fitted_illuminance: np.ndarray = (
        fitted_intercept
        + fitted_beta / (scaled_product + fitted_gamma)
    )
    residual_sum_squares: float = float(
        np.sum((valid_illuminance - fitted_illuminance) ** 2)
    )
    total_sum_squares: float = float(
        np.sum((valid_illuminance - np.mean(valid_illuminance)) ** 2)
    )
    r_squared: float = (
        1.0 - residual_sum_squares / total_sum_squares
        if total_sum_squares > 0
        else np.nan
    )
    return {
        "intercept": fitted_intercept,
        "numerator": fitted_beta * product_scale,
        "denominator_offset": fitted_gamma * product_scale,
        "r_squared": r_squared,
        "residual_scale": residual_scale,
        "n_samples": n_samples,
    }


def _plot_frame_saturation_calibration(
    results_by_recording: Sequence[dict[str, Any]],
    minimum_correlation: float = 0.9,
    maximum_saturation_percent: float | None = None,
    initial_samples_to_exclude: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Plot calibration points by diagnostics and selected activity groups.

    Points above ``maximum_saturation_percent`` are excluded from both the
    scatters and robust contextual shifted-reciprocal AGC fit. Both axes
    use logarithmic scales, so nonpositive camera products and illuminance
    values are omitted. The first ``initial_samples_to_exclude`` observations
    from each recording are also omitted from both plotting and fitting. The
    derivative is calculated within each complete video before these filters
    are applied. A saturation threshold of ``None`` retains all points.
    The displayed untransformed x/y vectors are also written to a MATLAB
    struct beside the MATLAB fitting routine.
    """

    if (
        maximum_saturation_percent is not None
        and (
            not np.isfinite(maximum_saturation_percent)
            or maximum_saturation_percent < 0.0
            or maximum_saturation_percent > 100.0
        )
    ):
        raise ValueError(
            "maximum_saturation_percent must be between 0 and 100."
        )
    if isinstance(initial_samples_to_exclude, bool) or not isinstance(
        initial_samples_to_exclude,
        (int, np.integer),
    ):
        raise TypeError("initial_samples_to_exclude must be an integer.")
    if initial_samples_to_exclude < 0:
        raise ValueError(
            "initial_samples_to_exclude must be greater than or equal to 0."
        )

    figure: Any
    axes: np.ndarray
    figure, axes = plt.subplots(3, 2, figsize=(20, 22))
    axes_flat: np.ndarray = axes.ravel()
    diagnostic_axes_flat: np.ndarray = axes_flat[0:2]
    activity_axes_flat: np.ndarray = axes_flat[2:4]
    high_correlation_axes_flat: np.ndarray = axes_flat[4:6]

    # Calculate the camera-score change separately within each complete video.
    # This preserves video boundaries and avoids threshold-induced jumps.
    point_camera_score_arrays: list[np.ndarray] = [
        np.asarray(
            recording_result["scatter_camera_score"],
            dtype=float,
        )
        for recording_result in results_by_recording
    ]
    point_saturation_arrays: list[np.ndarray] = [
        np.asarray(
            recording_result["scatter_saturation_percent"],
            dtype=float,
        )
        for recording_result in results_by_recording
    ]
    point_derivative_arrays: list[np.ndarray] = []
    for camera_score in point_camera_score_arrays:
        if camera_score.size == 0:
            camera_score_derivative: np.ndarray = np.array([], dtype=float)
        elif camera_score.size == 1:
            camera_score_derivative = np.zeros(1, dtype=float)
        else:
            camera_score_derivative = np.gradient(camera_score)
        point_derivative_arrays.append(camera_score_derivative)

    # Build one mask per recording. All panels use these exact observations,
    # while each panel's marker color is controlled by only its own metric.
    display_masks: list[np.ndarray] = []
    retained_camera_score_arrays: list[np.ndarray] = []
    retained_illuminance_arrays: list[np.ndarray] = []
    retained_saturation_arrays: list[np.ndarray] = []
    retained_derivative_arrays: list[np.ndarray] = []
    for recording_index, recording_result in enumerate(results_by_recording):
        recording_label: str = _recording_name_from_path(
            recording_result["recording_path"]
        )
        camera_score: np.ndarray = point_camera_score_arrays[recording_index]
        filtered_illuminance: np.ndarray = np.asarray(
            recording_result["scatter_illuminance"],
            dtype=float,
        )
        saturation_percent: np.ndarray = point_saturation_arrays[
            recording_index
        ]
        camera_score_derivative: np.ndarray = point_derivative_arrays[
            recording_index
        ]
        if not (
            camera_score.shape
            == filtered_illuminance.shape
            == saturation_percent.shape
            == camera_score_derivative.shape
        ):
            raise ValueError(
                "Cached calibration vectors must have matching shapes for "
                f"{recording_label}."
            )

        display_valid: np.ndarray = (
            np.isfinite(camera_score)
            & (camera_score > 0)
            & np.isfinite(filtered_illuminance)
            & (filtered_illuminance > 0)
            & np.isfinite(saturation_percent)
            & np.isfinite(camera_score_derivative)
        )
        if maximum_saturation_percent is not None:
            display_valid &= (
                saturation_percent <= maximum_saturation_percent
            )

        # Apply the index-based start exclusion independently to every video.
        # Clipping naturally handles recordings shorter than the requested N.
        excluded_sample_count: int = min(
            initial_samples_to_exclude,
            display_valid.size,
        )
        display_valid[:excluded_sample_count] = False
        display_masks.append(display_valid)
        if np.any(display_valid):
            retained_camera_score_arrays.append(
                camera_score[display_valid]
            )
            retained_illuminance_arrays.append(
                filtered_illuminance[display_valid]
            )
            retained_saturation_arrays.append(
                saturation_percent[display_valid]
            )
            retained_derivative_arrays.append(
                camera_score_derivative[display_valid]
            )

    # Scale saturation across all displayed videos so equal values receive the
    # same color regardless of recording identity.
    if len(retained_saturation_arrays) == 0:
        saturation_minimum: float = 0.0
        saturation_maximum: float = (
            maximum_saturation_percent
            if maximum_saturation_percent is not None
            else 100.0
        )
    else:
        finite_saturation: np.ndarray = np.concatenate(
            retained_saturation_arrays
        )
        saturation_minimum = float(np.min(finite_saturation))
        saturation_maximum = float(np.max(finite_saturation))
        if saturation_minimum == saturation_maximum:
            saturation_padding: float = max(0.5, saturation_minimum * 0.05)
            saturation_minimum = max(0.0, saturation_minimum - saturation_padding)
            saturation_upper_bound: float = (
                maximum_saturation_percent
                if maximum_saturation_percent is not None
                else 100.0
            )
            saturation_maximum = min(
                saturation_upper_bound,
                saturation_maximum + saturation_padding,
            )

    saturation_normalization: Normalize = Normalize(
        vmin=saturation_minimum,
        vmax=saturation_maximum,
    )
    saturation_colormap: Colormap = plt.get_cmap("viridis")

    # Center the derivative color scale at zero with symmetric robust bounds.
    # The 99th percentile prevents a few large changes from flattening the
    # visible color variation among the remaining observations.
    if len(retained_derivative_arrays) == 0:
        derivative_color_limit: float = 1.0
    else:
        retained_derivative: np.ndarray = np.concatenate(
            retained_derivative_arrays
        )
        derivative_color_limit = float(
            np.percentile(np.abs(retained_derivative), 99.0)
        )
        if (
            not np.isfinite(derivative_color_limit)
            or derivative_color_limit == 0.0
        ):
            derivative_color_limit = 1.0

    derivative_normalization: TwoSlopeNorm = TwoSlopeNorm(
        vmin=-derivative_color_limit,
        vcenter=0.0,
        vmax=derivative_color_limit,
    )
    derivative_colormap: Colormap = plt.get_cmap("RdBu_r")

    selected_recordings: list[dict[str, Any]] = [
        recording_result
        for recording_result in results_by_recording
        if np.isfinite(recording_result["model_result"]["best_correlation"])
        and recording_result["model_result"]["best_correlation"]
        >= minimum_correlation
    ]
    selected_recording_indices: dict[int, int] = {
        id(recording_result): recording_index
        for recording_index, recording_result in enumerate(selected_recordings)
    }
    high_correlation_colors: list[tuple[float, float, float]] = (
        _distinct_activity_colors(len(selected_recordings))
    )
    activity_markers: tuple[str, ...] = (
        "o",
        "s",
        "^",
        "D",
        "v",
        "P",
        "X",
        "<",
        ">",
        "*",
    )

    # All panels use the same display mask for a given recording. The first
    # row uses continuous diagnostic colors, the second row highlights focal
    # activity groups, and the third row shows the correlation-qualified
    # calibration subset by recording and primary gain control.
    read_legend_labels: set[str] = set()
    sit_biopond_legend_labels: set[str] = set()
    agc_type_legend_labels: set[str] = set()
    high_correlation_camera_scores: list[np.ndarray] = []
    high_correlation_illuminance: list[np.ndarray] = []
    for recording_index, recording_result in enumerate(results_by_recording):
        camera_score: np.ndarray = point_camera_score_arrays[recording_index]
        filtered_illuminance: np.ndarray = np.asarray(
            recording_result["scatter_illuminance"],
            dtype=float,
        )
        saturation_percent: np.ndarray = point_saturation_arrays[
            recording_index
        ]
        camera_score_derivative: np.ndarray = point_derivative_arrays[
            recording_index
        ]
        display_valid: np.ndarray = display_masks[recording_index]
        diagnostic_axes_flat[0].scatter(
            camera_score[display_valid],
            filtered_illuminance[display_valid],
            c=saturation_percent[display_valid],
            cmap=saturation_colormap,
            norm=saturation_normalization,
            marker="o",
            s=78,
            alpha=1.0,
            edgecolors="none",
            linewidths=0.0,
        )
        diagnostic_axes_flat[1].scatter(
            camera_score[display_valid],
            filtered_illuminance[display_valid],
            c=camera_score_derivative[display_valid],
            cmap=derivative_colormap,
            norm=derivative_normalization,
            marker="o",
            s=78,
            alpha=1.0,
            edgecolors="none",
            linewidths=0.0,
        )

        if np.any(display_valid):
            subject_id: str
            activity: str
            subject_id, activity = _recording_identity_from_path(
                recording_result["recording_path"]
            )

            # Highlight subjects 1029 and 1034 among the read recordings.
            # Every non-read activity remains a black background observation.
            if activity != "read":
                read_color: str = "black"
                read_label: str = "Non-read activities"
                read_zorder: int = 1
            elif subject_id in {"1029", "1034"}:
                read_color = "tab:red"
                read_label = "1029 and 1034 read"
                read_zorder = 3
            else:
                read_color = "tab:blue"
                read_label = "Other read activities"
                read_zorder = 2
            read_display_label: str | None = (
                read_label
                if read_label not in read_legend_labels
                else None
            )
            read_legend_labels.add(read_label)
            activity_axes_flat[0].scatter(
                camera_score[display_valid],
                filtered_illuminance[display_valid],
                color=read_color,
                label=read_display_label,
                marker="o",
                s=78,
                alpha=1.0,
                edgecolors="none",
                linewidths=0.0,
                zorder=read_zorder,
            )

            # Apply the same grouping to sitBiopond, with subject 1034 as the
            # red focal recording and every other sitBiopond recording blue.
            if activity != "sitBiopond":
                sit_color: str = "black"
                sit_label: str = "Non-sitBiopond activities"
                sit_zorder: int = 1
            elif subject_id == "1034":
                sit_color = "tab:red"
                sit_label = "1034 sitBiopond"
                sit_zorder = 3
            else:
                sit_color = "tab:blue"
                sit_label = "Other sitBiopond activities"
                sit_zorder = 2
            sit_display_label: str | None = (
                sit_label
                if sit_label not in sit_biopond_legend_labels
                else None
            )
            sit_biopond_legend_labels.add(sit_label)
            activity_axes_flat[1].scatter(
                camera_score[display_valid],
                filtered_illuminance[display_valid],
                color=sit_color,
                label=sit_display_label,
                marker="o",
                s=78,
                alpha=1.0,
                edgecolors="none",
                linewidths=0.0,
                zorder=sit_zorder,
            )

            selected_index: int | None = selected_recording_indices.get(
                id(recording_result)
            )
            if selected_index is not None:
                if (
                    "analog_gain" in recording_result
                    and "digital_gain" in recording_result
                ):
                    primary_gain: str | None = _primary_gain_control(
                        recording_result["analog_gain"],
                        recording_result["digital_gain"],
                    )
                else:
                    primary_gain = None

                if primary_gain is None:
                    agc_type_color = "0.45"
                    agc_type_label = "AGC type unavailable"
                elif primary_gain == "analog":
                    agc_type_color = "red"
                    agc_type_label = "Analog-gain primary"
                else:
                    agc_type_color = "black"
                    agc_type_label = "Digital-gain primary"
                agc_type_display_label: str | None = (
                    agc_type_label
                    if agc_type_label not in agc_type_legend_labels
                    else None
                )
                agc_type_legend_labels.add(agc_type_label)
                recording_label: str = _recording_name_from_path(
                    recording_result["recording_path"]
                )
                high_correlation_axes_flat[0].scatter(
                    camera_score[display_valid],
                    filtered_illuminance[display_valid],
                    color=high_correlation_colors[selected_index],
                    marker=activity_markers[
                        selected_index % len(activity_markers)
                    ],
                    s=48,
                    alpha=0.8,
                    edgecolors="white",
                    linewidths=0.25,
                    label=recording_label,
                )
                high_correlation_axes_flat[1].scatter(
                    camera_score[display_valid],
                    filtered_illuminance[display_valid],
                    color=agc_type_color,
                    s=24,
                    alpha=0.35,
                    label=agc_type_display_label,
                )
                high_correlation_camera_scores.append(
                    camera_score[display_valid]
                )
                high_correlation_illuminance.append(
                    filtered_illuminance[display_valid]
                )

    nonempty_high_correlation_scores: list[np.ndarray] = [
        camera_score
        for camera_score in high_correlation_camera_scores
        if camera_score.size > 0
    ]
    nonempty_high_correlation_illuminance: list[np.ndarray] = [
        illuminance
        for illuminance in high_correlation_illuminance
        if illuminance.size > 0
    ]
    if len(nonempty_high_correlation_scores) > 0:
        combined_high_correlation_camera_score: np.ndarray = np.concatenate(
            nonempty_high_correlation_scores
        )
        combined_high_correlation_illuminance: np.ndarray = np.concatenate(
            nonempty_high_correlation_illuminance
        )
        line_fit: dict[str, Any] = fit_illuminance_on_camera_score(
            combined_high_correlation_camera_score,
            combined_high_correlation_illuminance,
        )
        slope: float = line_fit["slope"]
        intercept: float = line_fit["intercept"]
        if np.isfinite(slope):
            score_fit: np.ndarray = np.geomspace(
                float(np.min(combined_high_correlation_camera_score)),
                float(np.max(combined_high_correlation_camera_score)),
                200,
            )
            illuminance_fit: np.ndarray = intercept + slope * score_fit
            positive_fit: np.ndarray = (
                np.isfinite(illuminance_fit) & (illuminance_fit > 0)
            )
            for ax in high_correlation_axes_flat:
                ax.plot(
                    score_fit[positive_fit],
                    illuminance_fit[positive_fit],
                    color="tab:blue",
                    linewidth=2.5,
                    label=(
                        "Shared calibration fit: "
                        f"illuminance = {intercept:.4g} "
                        f"+ {slope:.4g} x score"
                    ),
                )
    # Fit one shifted reciprocal model to the correlation-qualified activities.
    # This retains the original raw-lux Huber objective and affine intercept;
    # only the added denominator offset differs from the original model.
    pooled_camera_scores: list[np.ndarray] = []
    pooled_illuminance: list[np.ndarray] = []
    for recording_index, recording_result in enumerate(results_by_recording):
        correlation: float = float(
            recording_result["model_result"]["best_correlation"]
        )
        if not np.isfinite(correlation) or correlation < minimum_correlation:
            continue

        camera_score: np.ndarray = point_camera_score_arrays[recording_index]
        filtered_illuminance: np.ndarray = np.asarray(
            recording_result["scatter_illuminance"],
            dtype=float,
        )
        fit_valid: np.ndarray = display_masks[recording_index]
        if np.any(fit_valid):
            pooled_camera_scores.append(camera_score[fit_valid])
            pooled_illuminance.append(filtered_illuminance[fit_valid])

    if len(pooled_camera_scores) > 0:
        combined_camera_score: np.ndarray = np.concatenate(
            pooled_camera_scores
        )
        combined_illuminance: np.ndarray = np.concatenate(pooled_illuminance)
        shifted_reciprocal_fit: dict[str, float | int] = (
            _fit_shifted_reciprocal_illuminance(
                combined_camera_score,
                combined_illuminance,
            )
        )
        fitted_intercept: float = float(
            shifted_reciprocal_fit["intercept"]
        )
        numerator: float = float(shifted_reciprocal_fit["numerator"])
        denominator_offset: float = float(
            shifted_reciprocal_fit["denominator_offset"]
        )
        if (
            np.isfinite(fitted_intercept)
            and np.isfinite(numerator)
            and np.isfinite(denominator_offset)
        ):
            score_fit: np.ndarray = np.geomspace(
                float(np.min(combined_camera_score)),
                float(np.max(combined_camera_score)),
                200,
            )
            illuminance_fit: np.ndarray = (
                fitted_intercept
                + numerator / (score_fit + denominator_offset)
            )
            positive_fit: np.ndarray = (
                np.isfinite(illuminance_fit) & (illuminance_fit > 0)
            )
            for ax in diagnostic_axes_flat:
                ax.plot(
                    score_fit[positive_fit],
                    illuminance_fit[positive_fit],
                    color="tab:blue",
                    linewidth=2.5,
                    zorder=4,
                )

    figure_title: str = "Camera AGC Product vs Filtered Illuminance Diagnostics"
    if maximum_saturation_percent is not None:
        figure_title += f" (saturation <= {maximum_saturation_percent:g}%)"
    figure.suptitle(figure_title, fontsize=16)
    diagnostic_axes_flat[0].set_title("Frame Spatial Saturation")
    diagnostic_axes_flat[1].set_title("Camera AGC Product First Derivative")
    activity_axes_flat[0].set_title("Read Activities")
    activity_axes_flat[1].set_title("SitBiopond Activities")
    high_correlation_axes_flat[0].set_title(
        f"Correlation >= {minimum_correlation:.2f}: Activities by Name"
    )
    high_correlation_axes_flat[1].set_title(
        f"Correlation >= {minimum_correlation:.2f}: Activities by AGC Type"
    )
    if len(read_legend_labels) > 0:
        _set_legend_top_right(activity_axes_flat[0], fontsize=9)
    if len(sit_biopond_legend_labels) > 0:
        _set_legend_top_right(activity_axes_flat[1], fontsize=9)
    if len(selected_recordings) > 0:
        _set_legend_top_right(high_correlation_axes_flat[0], fontsize=8)
    if len(agc_type_legend_labels) > 0:
        _set_legend_top_right(high_correlation_axes_flat[1], fontsize=9)
    displayed_camera_score: np.ndarray = (
        np.concatenate(retained_camera_score_arrays)
        if len(retained_camera_score_arrays) > 0
        else np.array([], dtype=float)
    )
    displayed_illuminance: np.ndarray = (
        np.concatenate(retained_illuminance_arrays)
        if len(retained_illuminance_arrays) > 0
        else np.array([], dtype=float)
    )
    # Export exactly the displayed point cloud in its untransformed units.
    # MATLAB receives one scalar struct with no fit or recording metadata.
    savemat(
        os.fspath(LINEAR_SCALE_MAT_OUTPUT_PATH),
        {
            "linear_scale_data": {
                "x_linear_scale": displayed_camera_score,
                "y_linear_scale": displayed_illuminance,
            }
        },
        do_compression=True,
        oned_as="column",
    )

    for ax in axes_flat:
        ax.set_xlabel("Camera AGC product (log scale)")
        ax.set_ylabel("Filtered illuminance (lux, log scale)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        if displayed_camera_score.size > 0:
            score_minimum: float = float(np.min(displayed_camera_score))
            score_maximum: float = float(np.max(displayed_camera_score))
            ax.set_xlim(
                score_minimum / 1.05,
                score_maximum * 1.05,
            )
        if displayed_illuminance.size > 0:
            illuminance_minimum: float = float(
                np.min(displayed_illuminance)
            )
            illuminance_maximum: float = float(
                np.max(displayed_illuminance)
            )
            ax.set_ylim(
                illuminance_minimum / 1.05,
                illuminance_maximum * 1.05,
            )
        ax.grid()

    saturation_mappable: ScalarMappable = ScalarMappable(
        norm=saturation_normalization,
        cmap=saturation_colormap,
    )
    saturation_mappable.set_array(
        np.concatenate(retained_saturation_arrays)
        if len(retained_saturation_arrays) > 0
        else np.array([], dtype=float)
    )
    saturation_colorbar: Any = figure.colorbar(
        saturation_mappable,
        ax=diagnostic_axes_flat[0],
        pad=0.02,
    )
    saturation_colorbar.set_label("Frame spatial saturation (%)")

    derivative_mappable: ScalarMappable = ScalarMappable(
        norm=derivative_normalization,
        cmap=derivative_colormap,
    )
    derivative_mappable.set_array(
        np.concatenate(retained_derivative_arrays)
        if len(retained_derivative_arrays) > 0
        else np.array([], dtype=float)
    )
    derivative_colorbar: Any = figure.colorbar(
        derivative_mappable,
        ax=diagnostic_axes_flat[1],
        pad=0.02,
    )
    derivative_colorbar.set_label(
        "Camera AGC product first derivative (per sample)"
    )

    figure.tight_layout(rect=(0, 0, 1, 0.95))
    _draw_figure_now(figure)
    return axes, activity_axes_flat


def _save_frame_saturation_calibration_data(
    results_by_recording: Sequence[dict[str, Any]],
    output_path: str | os.PathLike[str],
) -> Path:
    """Save all processed inputs needed to recreate the saturation plot."""

    destination: Path = Path(output_path).expanduser()
    if destination.suffix.lower() != ".npz":
        destination = destination.with_suffix(".npz")
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Store each recording under indexed keys so vectors may retain different
    # lengths without object arrays or pickle-based serialization.
    archive_values: dict[str, np.ndarray] = {
        "format_version": np.asarray([2], dtype=int),
        "recording_count": np.asarray([len(results_by_recording)], dtype=int),
        "camera_score_definition": np.asarray(
            "Again * Dgain * exposure"
        ),
    }
    for recording_index, recording_result in enumerate(results_by_recording):
        key_suffix: str = str(recording_index)
        archive_values[f"recording_path_{key_suffix}"] = np.asarray(
            recording_result["recording_path"]
        )
        archive_values[f"camera_score_{key_suffix}"] = np.asarray(
            recording_result["scatter_camera_score"],
            dtype=float,
        )
        archive_values[f"filtered_illuminance_{key_suffix}"] = np.asarray(
            recording_result["scatter_illuminance"],
            dtype=float,
        )
        archive_values[f"frame_saturation_percent_{key_suffix}"] = np.asarray(
            recording_result["scatter_saturation_percent"],
            dtype=float,
        )
        archive_values[f"model_correlation_{key_suffix}"] = np.asarray(
            [recording_result["model_result"]["best_correlation"]],
            dtype=float,
        )

    np.savez_compressed(destination, **archive_values)
    return destination


def plot_frame_saturation_from_processed_data(
    processed_data_path: str | os.PathLike[str] = FRAME_SATURATION_DATA_PATH,
    maximum_saturation_percent: float | None = 40.0,
    minimum_correlation: float = 0.9,
    initial_samples_to_exclude: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Plot frame saturation calibration from previously processed recording data.

    The plotted linear-scale values are written to
    ``code/fitAGCtoIlluminance/camera_agc_illuminance_linear_scale.mat``
    as the fields ``x_linear_scale`` and ``y_linear_scale`` of
    ``linear_scale_data``.

    Args:
        processed_data_path: NumPy archive written by
            :func:`fit_agc_to_illuminance`.
        maximum_saturation_percent: Inclusive maximum frame saturation to
            display and fit. Use ``None`` to retain every processed point.
        minimum_correlation: Minimum recording-level model correlation used
            when fitting the contextual robust calibration line.
        initial_samples_to_exclude: Number of observations to omit from the
            start of every recording before plotting and fitting.

    Returns:
        Dashboard axes and the activity-highlight row axes.
    """

    source_path: Path = Path(processed_data_path).expanduser()
    if not source_path.is_file():
        raise FileNotFoundError(
            f"Processed frame-saturation data does not exist:\n{source_path}"
        )

    # Reconstruct the minimal per-recording structure expected by the plotting
    # helper. Loading the numeric archive requires no raw recordings or MATLAB.
    results_by_recording: list[dict[str, Any]] = []
    archive: NpzFile
    with np.load(source_path, allow_pickle=False) as archive:
        format_version: int = int(archive["format_version"][0])
        if format_version not in (1, 2):
            raise ValueError(
                f"Unsupported frame-saturation data version: {format_version}."
            )

        recording_count: int = int(archive["recording_count"][0])
        for recording_index in range(recording_count):
            key_suffix: str = str(recording_index)
            cached_camera_score: np.ndarray = np.asarray(
                archive[f"camera_score_{key_suffix}"],
                dtype=float,
            )

            # Version 1 stored the reciprocal brightness score. Recover the
            # AGC product exactly so old caches remain immediately usable.
            if format_version == 1:
                reciprocal_valid: np.ndarray = (
                    np.isfinite(cached_camera_score)
                    & (cached_camera_score > 0)
                )
                camera_score: np.ndarray = np.full(
                    cached_camera_score.shape,
                    np.nan,
                    dtype=float,
                )
                camera_score[reciprocal_valid] = (
                    1.0 / cached_camera_score[reciprocal_valid]
                )
            else:
                camera_score = cached_camera_score

            recording_result: dict[str, Any] = {
                "recording_path": str(
                    archive[f"recording_path_{key_suffix}"].item()
                ),
                "scatter_camera_score": camera_score,
                "scatter_illuminance": np.asarray(
                    archive[f"filtered_illuminance_{key_suffix}"],
                    dtype=float,
                ),
                "scatter_saturation_percent": np.asarray(
                    archive[f"frame_saturation_percent_{key_suffix}"],
                    dtype=float,
                ),
                "model_result": {
                    "best_correlation": float(
                        archive[f"model_correlation_{key_suffix}"][0]
                    )
                },
            }
            results_by_recording.append(recording_result)

    return _plot_frame_saturation_calibration(
        results_by_recording,
        minimum_correlation=minimum_correlation,
        maximum_saturation_percent=maximum_saturation_percent,
        initial_samples_to_exclude=initial_samples_to_exclude,
    )


def _plot_channel_count_activity_comparison(
    results_by_recording: Sequence[dict[str, Any]],
    minimum_correlation: float = 0.9,
) -> np.ndarray:
    """Compare the activity point cloud using eight versus nine channels."""

    figure: Any
    axes: np.ndarray
    figure, axes = plt.subplots(1, 2, figsize=(20, 8))
    axes_flat: np.ndarray = axes.ravel()

    activity_colors: list[tuple[float, float, float]] = _distinct_activity_colors(
        len(results_by_recording)
    )
    activity_markers: tuple[str, ...] = (
        "o",
        "s",
        "^",
        "D",
        "v",
        "P",
        "X",
        "<",
        ">",
        "*",
    )
    variant_settings: tuple[
        tuple[str, str, str, str],
        tuple[str, str, str, str],
    ] = (
        (
            "8 Channels",
            "model_result",
            "scatter_camera_score",
            "scatter_illuminance",
        ),
        (
            "9 Channels",
            "model_result_9_channels",
            "scatter_camera_score_9_channels",
            "scatter_illuminance_9_channels",
        ),
    )
    displayed_camera_scores: list[np.ndarray] = []
    displayed_illuminance: list[np.ndarray] = []

    # Build each panel independently because changing the channel average can
    # change both the shared lag and which recordings pass the threshold.
    for ax, (
        title,
        model_key,
        camera_score_key,
        illuminance_key,
    ) in zip(axes_flat, variant_settings):
        pooled_camera_scores: list[np.ndarray] = []
        pooled_illuminance: list[np.ndarray] = []

        for recording_index, recording_result in enumerate(results_by_recording):
            model_result: dict[str, Any] = recording_result[model_key]
            correlation: float = float(model_result["best_correlation"])
            if not np.isfinite(correlation) or correlation < minimum_correlation:
                continue

            camera_score: np.ndarray = np.asarray(
                recording_result[camera_score_key],
                dtype=float,
            )
            filtered_illuminance: np.ndarray = np.asarray(
                recording_result[illuminance_key],
                dtype=float,
            )
            valid: np.ndarray = (
                np.isfinite(camera_score) & np.isfinite(filtered_illuminance)
            )
            camera_score = camera_score[valid]
            filtered_illuminance = filtered_illuminance[valid]
            if camera_score.size == 0:
                continue

            recording_label: str = _recording_name_from_path(
                recording_result["recording_path"]
            )
            ax.scatter(
                camera_score,
                filtered_illuminance,
                color=activity_colors[recording_index],
                marker=activity_markers[
                    recording_index % len(activity_markers)
                ],
                s=18,
                alpha=0.8,
                edgecolors="white",
                linewidths=0.25,
                label=recording_label,
            )
            pooled_camera_scores.append(camera_score)
            pooled_illuminance.append(filtered_illuminance)
            displayed_camera_scores.append(camera_score)
            displayed_illuminance.append(filtered_illuminance)

        # Fit the same robust affine calibration used by the original
        # activity panel, but fit it separately for each option.
        if len(pooled_camera_scores) > 0:
            combined_camera_score: np.ndarray = np.concatenate(
                pooled_camera_scores
            )
            combined_illuminance: np.ndarray = np.concatenate(
                pooled_illuminance
            )
            line_fit: dict[str, Any] = fit_illuminance_on_camera_score(
                combined_camera_score,
                combined_illuminance,
            )
            slope: float = float(line_fit["slope"])
            intercept: float = float(line_fit["intercept"])
            if np.isfinite(slope):
                score_fit: np.ndarray = np.linspace(
                    float(np.min(combined_camera_score)),
                    float(np.max(combined_camera_score)),
                    200,
                )
                fitted_illuminance: np.ndarray = (
                    intercept + slope * score_fit
                )
                ax.plot(
                    score_fit,
                    fitted_illuminance,
                    color="tab:blue",
                    linewidth=2.5,
                    label=(
                        "Shared calibration fit: "
                        f"illuminance = {intercept:.4g} "
                        f"+ {slope:.4g} x score"
                    ),
                )
                displayed_camera_scores.append(score_fit)
                displayed_illuminance.append(fitted_illuminance)

        ax.set_title(title)
        ax.set_xlabel("Camera AGC product")
        ax.set_ylabel("Filtered illuminance (lux)")
        ax.grid()
        if len(pooled_camera_scores) > 0:
            _set_legend_top_right(ax, fontsize=8)

    # Use common nonnegative limits so apparent differences between panels
    # come from the channel calculation rather than independent autoscaling.
    maximum_camera_score: float = (
        max(float(np.max(values)) for values in displayed_camera_scores)
        if len(displayed_camera_scores) > 0
        else 1.0
    )
    maximum_illuminance: float = (
        max(float(np.max(values)) for values in displayed_illuminance)
        if len(displayed_illuminance) > 0
        else 1.0
    )
    for ax in axes_flat:
        ax.set_xlim(0.0, maximum_camera_score * 1.05)
        ax.set_ylim(0.0, maximum_illuminance * 1.05)

    figure.suptitle(
        "Activity point clouds: 8-channel vs 9-channel illuminance "
        f"(correlation >= {minimum_correlation:.2f})",
        fontsize=16,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    _draw_figure_now(figure)
    return axes


def _plot_all_recordings_diagnostics(
    results_by_recording: Sequence[dict[str, Any]],
    shared_model_result: dict[str, Any],
    calibration_result: dict[str, Any],
) -> np.ndarray:
    """
    Plot combined diagnostics across all recordings in one dashboard figure.
    """

    figure, axes = plt.subplots(4, 3, figsize=(18, 17))
    axes_flat: np.ndarray = axes.ravel()
    figure.suptitle("All Videos", fontsize=16)

    for recording_result in results_by_recording:
        label: str = _recording_name_from_path(recording_result["recording_path"])
        axes_flat[0].plot(
            recording_result["camera_time_relative"],
            recording_result["camera_score"],
            label=label,
        )
        axes_flat[1].plot(
            recording_result["as_time_relative"],
            recording_result["illuminance"],
            label=label,
        )
        axes_flat[2].scatter(
            recording_result["scatter_camera_score"],
            recording_result["scatter_illuminance"],
            s=8,
            alpha=0.35,
            label=label,
        )
        model_result: dict[str, Any] = recording_result["model_result"]
        relative_ms_time: np.ndarray = (
            model_result["time"] - model_result["time"][0]
        )
        axes_flat[3].plot(
            relative_ms_time,
            model_result["camera"],
            alpha=0.8,
            label=f"{label} camera",
        )
        axes_flat[3].plot(
            relative_ms_time,
            model_result["illuminance"],
            alpha=0.6,
            linestyle="--",
            label=f"{label} illum",
        )
        axes_flat[4].plot(
            relative_ms_time,
            model_result["prediction"],
            label=label,
        )
        axes_flat[5].plot(
            model_result["candidate_lags"],
            model_result["correlations"],
            label=label,
        )
        saturation_percent: np.ndarray = (
            np.asarray(recording_result["video_saturation"], dtype=float)
            * 100.0
        )
        axes_flat[9].plot(
            recording_result["saturation_time_relative"],
            saturation_percent,
            label=label,
        )

    axes_flat[0].set_title("Camera AGC product")
    axes_flat[0].set_xlabel("Time (s)")
    axes_flat[0].set_ylabel("Again x Dgain x exposure")
    axes_flat[1].set_title("Calibrated illuminance")
    axes_flat[1].set_xlabel("Time (s)")
    axes_flat[1].set_ylabel("Illuminance (lux)")
    axes_flat[2].set_title("Camera AGC product vs illuminance")
    axes_flat[2].set_xlabel("Camera AGC product")
    axes_flat[2].set_ylabel("Filtered illuminance (lux)")
    axes_flat[2].set_xlim(left=0)
    axes_flat[2].set_ylim(bottom=0)
    axes_flat[3].set_title("Normalized inputs")
    axes_flat[3].set_xlabel("Time (s)")
    axes_flat[3].set_ylabel("Z-score")
    axes_flat[4].set_title("Shared empirical-model predictions")
    axes_flat[4].set_xlabel("Time (s)")
    axes_flat[4].set_ylabel("Predicted brightness-oriented AGC product")
    axes_flat[5].set_title("Lag correlation by recording")
    axes_flat[5].set_xlabel("Lag (s)")
    axes_flat[5].set_ylabel("Correlation")
    _plot_shared_lag_summary(axes_flat[6], shared_model_result)
    _plot_calibration_fit(axes_flat[7], calibration_result)

    # The final panel makes poor model fits immediately visible without
    # requiring the user to inspect the CSV first.
    recording_labels: list[str] = [
        _recording_name_from_path(result["recording_path"])
        for result in results_by_recording
    ]
    recording_correlations: list[float] = [
        float(result["model_result"]["best_correlation"])
        for result in results_by_recording
    ]
    axes_flat[8].bar(recording_labels, recording_correlations, label="Correlation")
    axes_flat[8].tick_params(axis="x", rotation=45, labelsize=7)
    axes_flat[8].set_ylabel("Brightness-oriented AGC-product correlation")
    axes_flat[8].set_title("Empirical-model correlation")
    axes_flat[9].set_title("World video saturation")
    axes_flat[9].set_xlabel("Time since video start (s)")
    axes_flat[9].set_ylabel("Saturated pixels (%)")
    axes_flat[9].set_ylim(-10.0, 110.0)
    axes_flat[10].axis("off")
    axes_flat[11].axis("off")

    for ax in axes_flat:
        ax.grid()
        _set_legend_top_right(ax, fontsize=7)

    figure.tight_layout(rect=(0, 0, 1, 0.96))
    _draw_figure_now(figure)
    return axes


### MAIN
def fit_agc_to_illuminance(
    recording_paths: Sequence[str],
    fixed_lag: float | None = None,
    illuminance_diagnostics: bool = False,
    illuminance_diagnostics_output_dir: (
        str | os.PathLike[str] | None
    ) = None,
    frame_saturation_data_output_path: (
        str | os.PathLike[str] | None
    ) = None,
) -> pd.DataFrame:
    """
    Fit camera AGC settings to environmental illuminance.

    GOAL
    ----
    Derive a calibration function relating world-camera AGC settings
    (analog gain, digital gain, and exposure) to environmental illuminance
    measured simultaneously by the minispectrometer.

    For each supplied environmental recording:

    1. Load world-camera AGC metadata.
    2. Calculate a camera AGC score:

        score = Again * Dgain * exposure

    3. Measure saturation in the corresponding processed world video.
    4. Load the corresponding minispectrometer data.
    5. Flatten the AS measurements across minispect chunks.
    6. Convert minispectrometer measurements to illuminance.
    7. Align camera and minispectrometer timebases.

    All recordings are then modeled with ONE temporal response:

    8. Apply the same empirical AGC kernel to every recording.
    9. Fit one shared lag across all recordings, or use ``fixed_lag``.
    10. Pool sparse filtered illuminance / camera-score observations.
    11. Robustly fit illuminance = intercept + slope * camera_score.
    12. Summarize each recording by subject ID and activity and plot the
        shared empirical-kernel model across all recordings, including video
        saturation time series.
    13. Build one calibration-selection dashboard with saturation-colored
        points, AGC-derivative-colored points, focal read / sitBiopond
        highlights, and the correlation-qualified calibration subset.
    14. Save the unfiltered point data so thresholded versions can be plotted
        without reprocessing.

    Diagnostic runs additionally repeat the shared temporal model with the
    legacy nine-channel average, save both activity point clouds per recording,
    and finish with a matched 8-versus-9-channel comparison figure.

    Finally, characterize the relationship:

        filtered_illuminance ~= f(camera_score)

    and use the fitted relationship to derive a model that estimates
    environmental illuminance from the recent history of camera AGC settings.
    The steady-state estimate is:

        illuminance = intercept + slope * camera_score

    The function accepts an iterable of raw recording chunk paths beneath
    ``FLIC_raw`` so individual recordings can easily be included or excluded
    from calibration. Set ``illuminance_diagnostics`` to compare the returned
    eight-channel illuminance against the legacy nine-channel calculation for
    every recording. Full vectors are written beneath
    ``illuminance_diagnostics_output_dir``; when omitted, the default is
    ``code/fitAGCtoIlluminance/cached_processing_data/illuminance_channel_diagnostics``.
    Plot-ready frame saturation data are saved to
    ``frame_saturation_data_output_path``; when omitted, the default is
    ``code/fitAGCtoIlluminance/cached_processing_data/frame_saturation_calibration_data.npz``.
    Use :func:`plot_frame_saturation_from_processed_data` to apply a
    saturation threshold later without rerunning this workflow.
    Every supplied recording path is processed and included in the saturation
    comparison plots.

    Returns a per-recording summary DataFrame sorted by model correlation,
    best fit first. The shared lag and pooled calibration are stored in
    ``DataFrame.attrs``. Recordings that failed illuminance conversion appear
    last with NaN metrics. The caller is responsible for saving the table.
    """

    results_by_recording: list[dict[str, Any]] = []

    # LOOP THROUGH RECORDINGS
    # =================================================================
    # Process every supplied recording independently so diagnostics can reveal
    # recording-specific issues without imposing an internal activity filter.

    for recording_path in tqdm(
        recording_paths,
        desc="Loading recordings",
        unit="recording",
    ):

        # Validate the required participant/activity/GKA structure before any
        # recording data are loaded or converted.
        _recording_identity_from_path(recording_path)

        if not os.path.isdir(recording_path):
            raise FileNotFoundError(
                f"Recording directory does not exist:\n{recording_path}"
            )

        # CAMERA / AGC
        # =============================================================
        # Load the world-camera AGC metadata needed to reconstruct the
        # camera-side brightness score.
        metadata: pd.DataFrame = world_util.world_metadata_from_chunks(
            recording_path,
            convert_to_seconds=True,
            verbose=False,
        )

        # Pull the camera timing and exposure-control terms into arrays so
        # the later alignment, filtering, and plotting steps operate on a
        # consistent NumPy representation.
        camera_time: np.ndarray = metadata["timestamp"].to_numpy(dtype=float)

        analog_gain: np.ndarray
        digital_gain: np.ndarray
        exposure: np.ndarray
        analog_gain, digital_gain, exposure = extract_agc_columns(metadata)

        # Put camera time onto a relative timebase for plotting.
        camera_time_relative: np.ndarray = camera_time - camera_time[0]


        # CALCULATE CAMERA LIGHT SCORE
        # =============================================================
        # More gain or longer exposure means the camera needed more help to
        # see the scene, so their product increases as illumination decreases.
        # Camera AGC score:
        #
        #     Dgain * Again * exposure
        #
        # Invalid / missing AGC values remain NaN.

        agc_product: np.ndarray = analog_gain * digital_gain * exposure
        camera_score: np.ndarray = np.full(agc_product.shape, np.nan, dtype=float)
        valid_camera: np.ndarray = (
            np.isfinite(agc_product) & (agc_product > 0)
        )
        camera_score[valid_camera] = agc_product[valid_camera]

        # MEASURE PROCESSED WORLD-VIDEO SATURATION
        # =============================================================
        # The raw GKA path maps to a mirrored processed directory containing
        # the decoded world-camera video used for frame-level saturation.
        world_video_path: str = _processed_world_video_path(recording_path)
        if not os.path.isfile(world_video_path):
            raise FileNotFoundError(
                f"Processed world video does not exist:\n{world_video_path}"
            )

        video_saturation: np.ndarray = video_io.mesaure_video_saturation(
            world_video_path,
            verbose=True,
        )
        world_video_frame_rate: float = float(
            video_io.inspect_video_FPS(world_video_path)
        )
        if world_video_frame_rate <= 0:
            raise RuntimeError(
                f"Processed world video has invalid FPS "
                f"{world_video_frame_rate}:\n{world_video_path}"
            )

        # Convert frame indices to elapsed seconds while retaining one
        # saturation proportion for every decoded video frame.
        saturation_time_relative: np.ndarray = (
            np.arange(video_saturation.size, dtype=float)
            / world_video_frame_rate
        )

        # MINISPECT
        # =============================================================
        # Minispect chunks provide the external light measurement. Loading
        # only "M" avoids spending memory on unrelated world-camera frames.

        ms_chunks: list[dict[str, Any]] = parse_chunks(
            recording_path,
            convert_time_units=True,
            convert_to_float=True,
            chunk_ranges={
                "M": (0, None)
            }
        )

        # COMBINE / FLATTEN AS MINISPECT CHUNKS
        # =============================================================
        # Flattening happens in one helper so the nested chunk structure does
        # not leak into the higher-level calibration workflow.
        as_values: np.ndarray
        as_time: np.ndarray
        as_values, as_time = flatten_ms_chunks(ms_chunks)

        # Put AS time onto a relative timebase for plotting.
        as_time_relative: np.ndarray = as_time - as_time[0]

        # OPTIONAL DATAFRAME FOR EASIER INSPECTION
        # =============================================================
        # The DataFrame is not required for modeling, but gives a compact
        # tabular view and summary statistics for debugging sensor values.

        as_column_names: list[str] = [f"AS_{ii}" for ii in range(as_values.shape[1])]
        as_df: pd.DataFrame = pd.DataFrame(as_values, columns=as_column_names)
        as_df.insert(0, "timestamp", as_time)

        # CONVERT MINISPECT COUNTS TO ILLUMINANCE
        # =============================================================
        # Store the calibrated light signal for the second pass. The shared
        # temporal model is fit only after every recording has illuminance.
        # The production conversion uses the eight filtered channels F1-F8.
        # Diagnostic mode additionally calculates the legacy F1-F8 + Clear
        # result without changing the returned illuminance.
        # A recording whose counts cannot be converted is flagged rather than
        # allowed to abort the run, and is then held out of every pooled fit
        # so it cannot influence the other recordings.
        diagnostic_output_path: Path | None = None
        if illuminance_diagnostics:
            diagnostic_output_dir: Path = (
                Path(illuminance_diagnostics_output_dir)
                if illuminance_diagnostics_output_dir is not None
                else ILLUMINANCE_DIAGNOSTICS_OUTPUT_DIR
            )
            diagnostic_output_path = (
                diagnostic_output_dir
                / f"{_recording_name_from_path(recording_path)}.npz"
            )

        illuminance: np.ndarray
        illuminance_errored: bool
        illuminance_error_message: str
        (
            illuminance,
            illuminance_errored,
            illuminance_error_message,
        ) = safe_ms_counts_to_illuminance(
            as_values,
            diagnostics=illuminance_diagnostics,
            diagnostics_output_path=diagnostic_output_path,
            timestamps=as_time,
        )

        if illuminance_errored:
            print(
                f"\n[WARNING] {_recording_name_from_path(recording_path)} could "
                f"not be converted to illuminance and is excluded from all "
                f"pooled fits.\n          {illuminance_error_message}\n"
            )

        # Retain the nine-channel average in memory only during diagnostic
        # runs so the alternate shared model can be fit after all recordings
        # have been loaded.
        illuminance_9_channels: np.ndarray | None = None
        if (
            illuminance_diagnostics
            and not illuminance_errored
            and diagnostic_output_path is not None
        ):
            illuminance_9_channels = _load_nine_channel_illuminance(
                diagnostic_output_path
            )

        # =============================================================
        # SAVE PER-RECORDING INPUTS FOR THE SHARED FIT
        # =============================================================
        # At this point each recording has camera AGC data and calibrated MS
        # illuminance, but no temporal model yet. The shared fit below adds
        # filtered illuminance, lags, predictions, and diagnostic axes.

        recording_result: dict[str, Any] = {
            "recording_path": recording_path,

            "camera_metadata": metadata,

            "camera_time": camera_time,
            "camera_time_relative": camera_time_relative,

            "analog_gain": analog_gain,
            "digital_gain": digital_gain,
            "exposure": exposure,

            "camera_score": camera_score,

            "world_video_path": world_video_path,
            "video_saturation": video_saturation,
            "saturation_time_relative": saturation_time_relative,

            "as_time": as_time,
            "as_time_relative": as_time_relative,
            "as_values": as_values,
            "as_dataframe": as_df,

            "illuminance": illuminance,
            "illuminance_9_channels": illuminance_9_channels,
            "illuminance_diagnostic_path": diagnostic_output_path,

            "errored": illuminance_errored,
            "error_message": illuminance_error_message,
        }

        results_by_recording.append(
            recording_result
        )

    if len(results_by_recording) == 0:
        return build_recording_summary_table(results_by_recording)

    # HOLD OUT RECORDINGS THAT COULD NOT BE CONVERTED
    # =================================================================
    # Errored recordings keep their error signal and still appear in the
    # summary table, but every fit, pooled axis limit, and dashboard below
    # is built only from the usable recordings.
    usable_recordings: list[dict[str, Any]] = [
        recording_result for recording_result in results_by_recording
        if not recording_result["errored"]
    ]
    errored_recordings: list[dict[str, Any]] = [
        recording_result for recording_result in results_by_recording
        if recording_result["errored"]
    ]

    if len(usable_recordings) == 0:
        summary_table: pd.DataFrame = build_recording_summary_table(
            results_by_recording
        )

        print(
            f"\n[ERROR] All {len(errored_recordings)} recording(s) failed "
            "illuminance conversion. No model was fit.\n"
        )
        print(
            summary_table[
                ["subject_id", "activity", "errored", "error_message"]
            ].to_string(index=False)
        )

        return summary_table

    # FIT ONE EMPIRICAL TEMPORAL MODEL ACROSS ALL RECORDINGS
    # =================================================================
    # The AGC response is a property of the camera system, so every recording
    # uses the same empirical kernel and the same lag. If fixed_lag is None,
    # one lag is selected by mean correlation across all usable recordings.
    shared_model_result: dict[str, Any] = characterize_shared_empirical_kernel(
        usable_recordings,
        os.fspath(AGC_KERNEL_MAT_PATH),
        lag_range=(-10, 10),
        n_lags=201,
        fixed_lag=fixed_lag,
    )

    # APPLY THE SHARED MODEL TO EVERY RECORDING
    # =================================================================
    # Correlation, physical-unit scatter samples, line fits, and diagnostics
    # all use the exact same lag-aligned empirical-kernel output.

    for recording_result, model_result in zip(
        tqdm(
            usable_recordings,
            desc="Applying shared empirical model",
            unit="recording",
        ),
        shared_model_result["per_recording_results"],
    ):
        as_time = recording_result["as_time"]
        camera_time = recording_result["camera_time"]
        camera_score = recording_result["camera_score"]

        # Put the lag-aligned kernel output back on the original minispect
        # timebase so calibration sampling and raw-signal plots share a clock.
        filtered_illuminance: np.ndarray = np.interp(
            as_time,
            model_result["time"],
            model_result["filtered_illuminance"],
            left=np.nan,
            right=np.nan,
        )

        # Put the camera AGC product onto the minispect timebase.
        # camera_time and as_time are already expressed in the same
        # recording-clock units, so interpolate directly using them.
        camera_score_on_as_time: np.ndarray = np.interp(
            as_time,
            camera_time,
            camera_score,
            left=np.nan,
            right=np.nan,
        )

        valid_scatter: np.ndarray = (
            np.isfinite(camera_score_on_as_time)
            & np.isfinite(filtered_illuminance)
        )
        scatter_camera_score: np.ndarray = camera_score_on_as_time[valid_scatter]
        scatter_illuminance: np.ndarray = filtered_illuminance[valid_scatter]

        # Match each calibration observation to the nearest decoded video
        # frame. Both elapsed times are measured from the camera recording
        # start, so each point receives that frame's spatial saturation.
        scatter_time_relative: np.ndarray = (
            as_time[valid_scatter] - camera_time[0]
        )
        scatter_saturation_percent: np.ndarray = (
            _video_saturation_at_times_percent(
                scatter_time_relative,
                recording_result["saturation_time_relative"],
                recording_result["video_saturation"],
            )
        )

        # Fit the robust affine calibration once so the figure and returned
        # summary report the same slope and intercept.
        absolute_fit: dict[str, Any] = fit_illuminance_on_camera_score(
            scatter_camera_score,
            scatter_illuminance,
        )

        recording_result["filtered_illuminance"] = filtered_illuminance
        recording_result["camera_score_on_as_time"] = camera_score_on_as_time
        recording_result["scatter_camera_score"] = scatter_camera_score
        recording_result["scatter_illuminance"] = scatter_illuminance
        recording_result["scatter_time_relative"] = scatter_time_relative
        recording_result["scatter_saturation_percent"] = (
            scatter_saturation_percent
        )
        recording_result["model_result"] = model_result
        recording_result["absolute_fit"] = absolute_fit

    # BUILD THE ALTERNATE NINE-CHANNEL ACTIVITY POINT CLOUDS
    # =================================================================
    # Diagnostic mode repeats the complete shared-lag model with the legacy
    # nine-channel illuminance. This produces the true x/y vectors that the
    # final activity plot would use under either channel-selection option.
    if illuminance_diagnostics:
        nine_channel_model_inputs: list[dict[str, Any]] = []
        for recording_result in usable_recordings:
            illuminance_9_channels = recording_result[
                "illuminance_9_channels"
            ]
            if illuminance_9_channels is None:
                raise RuntimeError(
                    "Nine-channel illuminance is missing for diagnostic "
                    f"recording {recording_result['recording_path']}."
                )

            nine_channel_model_input: dict[str, Any] = recording_result.copy()
            nine_channel_model_input["illuminance"] = illuminance_9_channels
            nine_channel_model_inputs.append(nine_channel_model_input)

        shared_model_result_9_channels: dict[str, Any] = (
            characterize_shared_empirical_kernel(
                nine_channel_model_inputs,
                os.fspath(AGC_KERNEL_MAT_PATH),
                lag_range=(-10, 10),
                n_lags=201,
                fixed_lag=fixed_lag,
            )
        )
        minimum_activity_correlation: float = 0.9

        for recording_result, model_result_9_channels in zip(
            usable_recordings,
            shared_model_result_9_channels["per_recording_results"],
        ):
            camera_score_9_channels: np.ndarray
            filtered_illuminance_9_channels: np.ndarray
            (
                camera_score_9_channels,
                filtered_illuminance_9_channels,
            ) = _activity_point_cloud_from_model(
                recording_result,
                model_result_9_channels,
            )

            # Keep the alternate model and its point vectors in memory so the
            # final side-by-side figure uses exactly the arrays saved below.
            recording_result["model_result_9_channels"] = (
                model_result_9_channels
            )
            recording_result["scatter_camera_score_9_channels"] = (
                camera_score_9_channels
            )
            recording_result["scatter_illuminance_9_channels"] = (
                filtered_illuminance_9_channels
            )

            diagnostic_path: Path | None = recording_result[
                "illuminance_diagnostic_path"
            ]
            if diagnostic_path is None:
                raise RuntimeError(
                    "Diagnostic output path is missing for "
                    f"{recording_result['recording_path']}."
                )

            _append_activity_point_cloud_diagnostics(
                diagnostic_path,
                recording_result["scatter_camera_score"],
                recording_result["scatter_illuminance"],
                float(shared_model_result["shared_lag"]),
                float(recording_result["model_result"]["best_correlation"]),
                camera_score_9_channels,
                filtered_illuminance_9_channels,
                float(shared_model_result_9_channels["shared_lag"]),
                float(model_result_9_channels["best_correlation"]),
                minimum_activity_correlation,
            )

    # BUILD PER-RECORDING DASHBOARDS ON COMMON AXES
    # =================================================================
    # Plotting waits until every recording has been modeled so each panel can
    # use x and y limits pooled across all recordings. Comparing dashboards
    # side by side is only meaningful when the axes do not rescale per
    # recording.
    shared_axis_limits: dict[
        str,
        dict[str, tuple[float, float] | None],
    ] = compute_shared_axis_limits(
        usable_recordings,
        padding_fraction=0.05,
    )
    shared_model_result["shared_axis_limits"] = shared_axis_limits

    for recording_result in tqdm(
        usable_recordings,
        desc="Plotting recording dashboards",
        unit="recording",
    ):
        recording_result["diagnostic_axes"] = _plot_recording_diagnostics(
            recording_result,
            shared_axis_limits,
        )

    # FINALIZE
    # =================================================================
    # These closing stages used to run silently after the last progress bar
    # finished, which made a slow combined plot look like a hang. Each stage
    # reports itself, and tqdm.write is used for messages so they print
    # above the bar instead of corrupting it.
    number_of_finalization_steps: int = 6 if illuminance_diagnostics else 5
    finalization_progress: Any = tqdm(
        total=number_of_finalization_steps,
        desc="Finalizing",
        unit="step",
    )

    # Pool one sample per second from the shared empirical model and fit the
    # steady-state illuminance implied by the camera AGC settings.
    finalization_progress.set_postfix_str("steady-state calibration")
    calibration_result: dict[str, Any] = fit_steady_state_illuminance_calibration(
        usable_recordings,
        sample_interval_seconds=1.0,
        filtered_illuminance_key="filtered_illuminance",
    )
    shared_model_result["calibration_result"] = calibration_result
    finalization_progress.update(1)

    # Per-recording dashboards are created above. This dashboard shows the
    # across-recording relationships, shared empirical fit, and video
    # saturation time series in one figure.
    finalization_progress.set_postfix_str("combined dashboard")
    shared_model_result["diagnostic_axes"] = _plot_all_recordings_diagnostics(
        usable_recordings,
        shared_model_result,
        calibration_result,
    )
    finalization_progress.update(1)

    # Render point-selection diagnostics, focal activity highlights, and the
    # correlation-qualified calibration subset in one dashboard.
    finalization_progress.set_postfix_str("calibration selection diagnostics")
    (
        shared_model_result["frame_saturation_axes"],
        shared_model_result["activity_highlight_axes"],
    ) = (
        _plot_frame_saturation_calibration(
            usable_recordings,
            minimum_correlation=0.9,
        )
    )
    finalization_progress.update(1)

    # Persist unfiltered, plot-ready point vectors once. Threshold exploration
    # can then reload this small archive without repeating video or sensor work.
    finalization_progress.set_postfix_str("saving frame-saturation data")
    resolved_frame_saturation_data_path: Path = (
        Path(frame_saturation_data_output_path)
        if frame_saturation_data_output_path is not None
        else FRAME_SATURATION_DATA_PATH
    )
    saved_frame_saturation_data_path: Path = (
        _save_frame_saturation_calibration_data(
            usable_recordings,
            resolved_frame_saturation_data_path,
        )
    )
    shared_model_result["frame_saturation_data_path"] = (
        saved_frame_saturation_data_path
    )
    finalization_progress.update(1)

    # End diagnostic runs with a direct activity-by-activity comparison of
    # the two channel-count options. Normal production runs omit this figure.
    if illuminance_diagnostics:
        finalization_progress.set_postfix_str("8-vs-9 activity comparison")
        shared_model_result["channel_count_comparison_axes"] = (
            _plot_channel_count_activity_comparison(
                usable_recordings,
                minimum_correlation=0.9,
            )
        )
        finalization_progress.update(1)

    # One row per recording keyed by subject ID and activity, holding the
    # shared-model correlation and the center-panel calibration line. Errored
    # recordings appear in the table so an exclusion is
    # visible, but only the usable rows carry fitted numbers.
    finalization_progress.set_postfix_str("summary table")
    summary_table: pd.DataFrame = build_recording_summary_table(results_by_recording)
    finalization_progress.update(1)

    tqdm.write(
        summary_table[
            [
                "subject_id",
                "activity",
                "errored",
                "correlation",
                "slope",
                "intercept",
                "camera_score_std",
            ]
        ].to_string(index=False)
    )

    if len(errored_recordings) > 0:
        tqdm.write(
            f"\n{len(errored_recordings)} of {len(results_by_recording)} "
            f"recording(s) were excluded from all fits:"
        )
        for recording_result in errored_recordings:
            tqdm.write(
                f"  {_recording_name_from_path(recording_result['recording_path'])}: "
                f"{recording_result['error_message']}"
            )

    finalization_progress.close()

    # ATTACH POOLED RESULTS TO THE RETURNED TABLE
    # =================================================================
    # DataFrame.attrs keeps the one shared lag and pooled calibration
    # available programmatically without repeating them in every CSV row.
    summary_table.attrs["model"] = "empirical AGC kernel"
    summary_table.attrs["camera_score_definition"] = (
        "Again * Dgain * exposure"
    )
    summary_table.attrs["shared_lag_s"] = shared_model_result["shared_lag"]
    summary_table.attrs["pooled_slope"] = calibration_result["slope"]
    summary_table.attrs["pooled_intercept"] = calibration_result["intercept"]
    summary_table.attrs["calibration_fit_method"] = calibration_result["fit_method"]
    summary_table.attrs["implied_illuminance_equation"] = (
        calibration_result["implied_illuminance_equation"]
    )
    summary_table.attrs["frame_saturation_data_path"] = os.fspath(
        saved_frame_saturation_data_path
    )

    tqdm.write(
        "\nPooled robust steady-state calibration:\n"
        f"illuminance = {calibration_result['intercept']:.6g} + "
        f"{calibration_result['slope']:.6g} * camera_score\n"
    )


    # =================================================================
    # RETURN THE SUMMARY TABLE
    # =================================================================
    # Sorted best fit first; persistence is left to the caller.

    return summary_table


def main() -> None:
    """Run the AGC-to-illuminance workflow for command-line recording paths."""

    recording_paths: list[str] = sys.argv[1:]
    if len(recording_paths) == 0:
        raise ValueError(
            "Provide one or more recording paths with format "
            ".../FLIC_<NUMERICALID>/<RECORDING_NAME>/GKA."
        )

    summary_table: pd.DataFrame = fit_agc_to_illuminance(recording_paths)
    print(f"Finished processing {len(summary_table)} recording(s).")


if __name__ == "__main__":
    main()
