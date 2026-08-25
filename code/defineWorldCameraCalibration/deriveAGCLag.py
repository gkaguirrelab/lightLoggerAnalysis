"""Derive one shared camera-response lag from raw light-logger recordings.

For each supplied ``GKA`` directory, this script loads the camera AGC product
and minispect illuminance, applies the empirical AGC kernel, and evaluates the
camera/minispect correlation across a common lag grid. The lag maximizing the
mean correlation across recordings is written to ``derived/cameraAGCLag.mat``.

Positive lag means that the camera response follows the minispect signal.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat, savemat


THIS_FILE: Path = Path(__file__).resolve()
PROJECT_ROOT: Path = THIS_FILE.parents[2]
LIGHT_LOGGER_ROOT: Path = PROJECT_ROOT.parent / "lightLogger"
PI_UTILITY_DIR: Path = LIGHT_LOGGER_ROOT / "raspberry_pi_firmware" / "utility"
SENSOR_UTILITY_DIR: Path = PROJECT_ROOT / "code" / "library" / "sensor_utility"
DEFAULT_KERNEL_PATH: Path = PROJECT_ROOT / "data" / "agc_empirical_kernels.mat"
DEFAULT_OUTPUT_PATH: Path = PROJECT_ROOT / "derived" / "cameraAGCLag.mat"
DEFAULT_LAG_RANGE_SECONDS: tuple[float, float] = (-10.0, 10.0)
DEFAULT_LAG_COUNT: int = 201


class InvalidRecordingError(ValueError):
    """A recording is readable but cannot contribute to lag estimation."""


@dataclass(frozen=True)
class RecordingSignals:
    """Camera and minispect signals needed by both calibration stages."""

    path: Path
    camera_time: np.ndarray
    camera_score: np.ndarray
    minispect_time: np.ndarray
    illuminance: np.ndarray


@dataclass(frozen=True)
class LagSearchResult:
    """Shared-lag result and the diagnostic curve supporting it."""

    output_path: Path
    shared_lag_seconds: float
    candidate_lags_seconds: np.ndarray
    mean_correlation_by_lag: np.ndarray
    recording_count: int
    skipped_recording_count: int


def start_matlab_engine() -> Any:
    """Start one MATLAB Engine for all minispect conversions in a run."""

    import matlab.engine

    # Project activation places the MATLAB calibration functions on the path.
    engine: Any = matlab.engine.start_matlab()
    engine.tbUseProject("lightLoggerAnalysis", nargout=0)
    return engine


def _load_raw_dependencies() -> tuple[Any, Any, Any]:
    """Import raw-recording utilities only when raw data are requested."""

    for import_path in (PI_UTILITY_DIR, SENSOR_UTILITY_DIR):
        import_path_string: str = os.fspath(import_path)
        if import_path_string not in sys.path:
            sys.path.append(import_path_string)

    from Pi_util import parse_chunks
    import world_util
    from ms_util import ms_counts_to_illuminance

    return parse_chunks, world_util, ms_counts_to_illuminance


def _flatten_minispect_chunks(
    chunks: Sequence[dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray]:
    """Flatten minispect AS values and timestamps from parsed chunks."""

    sample_chunks: list[np.ndarray] = []
    time_chunks: list[np.ndarray] = []
    for chunk in chunks:
        samples: np.ndarray = chunk["M"]["v"]["AS"]
        times: np.ndarray = chunk["M"]["t"]["AS"]
        if len(samples) == 0:
            continue
        if len(samples) != len(times):
            raise InvalidRecordingError(
                "Minispect sample and timestamp counts do not match."
            )
        sample_chunks.append(samples)
        time_chunks.append(times)

    if not sample_chunks:
        raise InvalidRecordingError("No minispect AS samples were found.")
    return (
        np.vstack(sample_chunks),
        np.concatenate(time_chunks).astype(float),
    )


def _extract_agc_product(metadata: Any) -> np.ndarray:
    """Return analog gain × digital gain × exposure from a known schema."""

    known_schemas: tuple[tuple[str, str, str], ...] = (
        ("Again", "Dgain", "exposure"),
        ("cameraAgain", "AGCDgain", "cameraExposure"),
    )
    selected_schema: tuple[str, str, str] | None = next(
        (
            schema
            for schema in known_schemas
            if all(column in metadata.columns for column in schema)
        ),
        None,
    )
    if selected_schema is None:
        raise InvalidRecordingError(
            "World metadata does not contain a recognized AGC schema."
        )

    analog_gain = metadata[selected_schema[0]].to_numpy(dtype=float)
    digital_gain = metadata[selected_schema[1]].to_numpy(dtype=float)
    exposure = metadata[selected_schema[2]].to_numpy(dtype=float)
    camera_score: np.ndarray = analog_gain * digital_gain * exposure
    camera_score[~np.isfinite(camera_score) | (camera_score <= 0)] = np.nan
    return camera_score


def load_recording_signals(
    recording_path: str | os.PathLike[str],
    *,
    matlab_engine: Any | None = None,
) -> RecordingSignals:
    """Load the camera score and calibrated minispect signal for one GKA."""

    path: Path = Path(recording_path).expanduser().resolve()
    if not path.is_dir() or path.name != "GKA":
        raise FileNotFoundError(
            f"Recording path must be an existing GKA directory:\n{path}"
        )

    parse_chunks, world_util, ms_counts_to_illuminance = (
        _load_raw_dependencies()
    )
    # Camera metadata provides one AGC product and timestamp per world frame.
    metadata: Any = world_util.world_metadata_from_chunks(
        os.fspath(path),
        convert_to_seconds=True,
        verbose=False,
    )
    camera_time: np.ndarray = metadata["timestamp"].to_numpy(dtype=float)
    camera_score: np.ndarray = _extract_agc_product(metadata)

    # Minispect samples live in M chunks and are converted to calibrated lux.
    chunks: list[dict[str, Any]] = parse_chunks(
        os.fspath(path),
        convert_time_units=True,
        convert_to_float=True,
        chunk_ranges={"M": (0, None)},
    )
    counts, minispect_time = _flatten_minispect_chunks(chunks)
    spectral_counts: np.ndarray = np.asarray(counts[:, :8], dtype=float)
    if not np.all(np.isfinite(spectral_counts) & (spectral_counts > 0)):
        raise InvalidRecordingError(
            "The first eight minispect channels must be finite and positive."
        )

    illuminance: np.ndarray = np.asarray(
        ms_counts_to_illuminance(counts, matlab_engine=matlab_engine),
        dtype=float,
    )
    if illuminance.shape != minispect_time.shape:
        raise InvalidRecordingError(
            "Illuminance and minispect timestamps do not match."
        )
    return RecordingSignals(
        path=path,
        camera_time=camera_time,
        camera_score=camera_score,
        minispect_time=minispect_time,
        illuminance=illuminance,
    )


def apply_empirical_kernel(
    time: np.ndarray,
    signal: np.ndarray,
    kernel_path: str | os.PathLike[str] = DEFAULT_KERNEL_PATH,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the empirical AGC kernel on its uniform timebase."""

    # Resample to the kernel's fixed interval before performing convolution.
    kernel_data: dict[str, Any] = loadmat(kernel_path)
    kernel_time: np.ndarray = np.squeeze(
        kernel_data["commonKernelTime"]
    ).astype(float)
    kernel: np.ndarray = np.squeeze(kernel_data["meanKernel"]).astype(float)
    time_step: float = float(np.median(np.diff(kernel_time)))
    uniform_time: np.ndarray = np.arange(time[0], time[-1], time_step)
    uniform_signal: np.ndarray = np.interp(uniform_time, time, signal)
    kernel_weights: np.ndarray = kernel / np.sum(kernel)
    # Edge padding avoids an artificial zero-valued transient at the start.
    padded_signal: np.ndarray = np.pad(
        uniform_signal,
        (len(kernel_weights) - 1, 0),
        mode="edge",
    )
    return uniform_time, np.convolve(
        padded_signal,
        kernel_weights,
        mode="valid",
    )


def lag_filtered_signal(
    time: np.ndarray,
    filtered_signal: np.ndarray,
    lag_seconds: float,
) -> np.ndarray:
    """Shift a signal; positive lag means the camera response occurs later."""

    return np.interp(
        time - lag_seconds,
        time,
        filtered_signal,
        left=np.nan,
        right=np.nan,
    )


def _zscore(values: np.ndarray) -> np.ndarray:
    """Return finite-sample z scores while preserving nonfinite positions."""

    values = np.asarray(values, dtype=float)
    finite: np.ndarray = np.isfinite(values)
    if np.sum(finite) < 2:
        raise InvalidRecordingError(
            "At least two finite samples are required for lag estimation."
        )
    standard_deviation: float = float(np.std(values[finite]))
    if standard_deviation == 0:
        raise InvalidRecordingError(
            "A constant signal cannot contribute to lag estimation."
        )
    result: np.ndarray = np.full(values.shape, np.nan, dtype=float)
    result[finite] = (
        values[finite] - float(np.mean(values[finite]))
    ) / standard_deviation
    return result


def _recording_lag_correlations(
    recording: RecordingSignals,
    candidate_lags_seconds: np.ndarray,
    kernel_path: str | os.PathLike[str],
) -> np.ndarray:
    """Calculate one recording's camera/minispect correlation by lag."""

    # Filtering models the AGC controller's temporal response to illumination.
    filtered_time, filtered_illuminance = apply_empirical_kernel(
        recording.minispect_time,
        recording.illuminance,
        kernel_path,
    )
    # Greater illumination requires less gain/exposure, hence the minus sign.
    camera_brightness: np.ndarray = -_zscore(recording.camera_score)
    filtered_brightness: np.ndarray = _zscore(filtered_illuminance)
    camera_on_filtered_time: np.ndarray = np.interp(
        filtered_time,
        recording.camera_time,
        camera_brightness,
        left=np.nan,
        right=np.nan,
    )

    correlations: np.ndarray = np.full(candidate_lags_seconds.shape, np.nan)
    for lag_index, lag_seconds in enumerate(candidate_lags_seconds):
        shifted_brightness: np.ndarray = lag_filtered_signal(
            filtered_time,
            filtered_brightness,
            float(lag_seconds),
        )
        valid: np.ndarray = (
            np.isfinite(camera_on_filtered_time)
            & np.isfinite(shifted_brightness)
        )
        if np.sum(valid) >= 3:
            correlations[lag_index] = np.corrcoef(
                camera_on_filtered_time[valid],
                shifted_brightness[valid],
            )[0, 1]
    return correlations


def _save_lag_result(result: LagSearchResult) -> None:
    """Write the compact MATLAB contract consumed by the next stage."""

    # Write-then-replace prevents a failed run from corrupting the prior MAT.
    result.output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path = result.output_path.with_suffix(".tmp.mat")
    savemat(
        temporary_path,
        {
            "cameraAGCLag": {
                "sharedLagSeconds": result.shared_lag_seconds,
                "candidateLagsSeconds": result.candidate_lags_seconds,
                "meanCorrelationByLag": result.mean_correlation_by_lag,
                "recordingCount": result.recording_count,
                "skippedRecordingCount": result.skipped_recording_count,
            }
        },
        do_compression=True,
        oned_as="column",
    )
    os.replace(temporary_path, result.output_path)


def derive_agc_lag(
    recording_paths: Sequence[str | os.PathLike[str]],
    *,
    output_path: str | os.PathLike[str] = DEFAULT_OUTPUT_PATH,
    kernel_path: str | os.PathLike[str] = DEFAULT_KERNEL_PATH,
    lag_range_seconds: tuple[float, float] = DEFAULT_LAG_RANGE_SECONDS,
    lag_count: int = DEFAULT_LAG_COUNT,
) -> LagSearchResult:
    """Estimate, save, and return one shared lag across raw recordings."""

    if not recording_paths:
        raise ValueError("At least one raw GKA recording path is required.")
    if (
        lag_count < 2
        or not np.all(np.isfinite(lag_range_seconds))
        or lag_range_seconds[0] >= lag_range_seconds[1]
    ):
        raise ValueError("The lag grid requires a finite range and at least 2 points.")

    resolved_recording_paths: list[Path] = [
        Path(path).expanduser().resolve() for path in recording_paths
    ]
    for recording_path in resolved_recording_paths:
        if not recording_path.is_dir() or recording_path.name != "GKA":
            raise FileNotFoundError(
                "Recording paths must be existing GKA directories:\n"
                f"{recording_path}"
            )

    # Every recording is evaluated on the same lag grid so their correlations
    # can be averaged without interpolation or recording-length weighting.
    candidate_lags: np.ndarray = np.linspace(
        lag_range_seconds[0],
        lag_range_seconds[1],
        lag_count,
    )
    correlation_rows: list[np.ndarray] = []
    skipped_recording_count: int = 0
    # Engine startup is expensive; reuse one instance for the complete run.
    matlab_engine: Any = start_matlab_engine()
    try:
        for recording_index, recording_path in enumerate(
            resolved_recording_paths,
            start=1,
        ):
            print(
                f"[{recording_index}/{len(resolved_recording_paths)}] "
                f"Loading {recording_path}"
            )
            try:
                recording: RecordingSignals = load_recording_signals(
                    recording_path,
                    matlab_engine=matlab_engine,
                )
                correlations: np.ndarray = _recording_lag_correlations(
                    recording,
                    candidate_lags,
                    kernel_path,
                )
            except InvalidRecordingError as error:
                skipped_recording_count += 1
                print(f"Skipping invalid recording: {error}", file=sys.stderr)
                continue
            correlation_rows.append(correlations)
    finally:
        matlab_engine.quit()

    if not correlation_rows:
        raise RuntimeError("No recording produced a valid lag-correlation curve.")
    correlations_by_recording: np.ndarray = np.vstack(correlation_rows)
    # Each valid recording contributes equally to the shared-lag decision.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_correlation: np.ndarray = np.nanmean(
            correlations_by_recording,
            axis=0,
        )
    if not np.any(np.isfinite(mean_correlation)):
        raise RuntimeError("No candidate lag produced a finite correlation.")

    best_index: int = int(np.nanargmax(mean_correlation))
    result = LagSearchResult(
        output_path=Path(output_path).expanduser().resolve(),
        shared_lag_seconds=float(candidate_lags[best_index]),
        candidate_lags_seconds=candidate_lags,
        mean_correlation_by_lag=mean_correlation,
        recording_count=len(correlation_rows),
        skipped_recording_count=skipped_recording_count,
    )
    _save_lag_result(result)
    return result


def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "recording_paths",
        nargs="+",
        help="Raw recording directories with the form .../activity/GKA.",
    )
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--kernel-path", type=Path, default=DEFAULT_KERNEL_PATH)
    parser.add_argument("--lag-min", type=float, default=-10.0)
    parser.add_argument("--lag-max", type=float, default=10.0)
    parser.add_argument("--lag-count", type=int, default=DEFAULT_LAG_COUNT)
    return parser


def main(argv: Sequence[str] | None = None) -> Path:
    """Run lag estimation and return the generated MATLAB path."""

    arguments = _build_argument_parser().parse_args(argv)
    result: LagSearchResult = derive_agc_lag(
        arguments.recording_paths,
        output_path=arguments.output_path,
        kernel_path=arguments.kernel_path,
        lag_range_seconds=(arguments.lag_min, arguments.lag_max),
        lag_count=arguments.lag_count,
    )
    print(
        f"Saved {result.shared_lag_seconds:g} s camera lag from "
        f"{result.recording_count} recording(s) to {result.output_path}"
    )
    return result.output_path


if __name__ == "__main__":
    main()
