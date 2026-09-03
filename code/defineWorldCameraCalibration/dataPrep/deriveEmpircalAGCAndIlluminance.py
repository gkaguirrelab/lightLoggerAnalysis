"""Export the empirical camera-AGC/minispect-illuminance point cloud.

The shared temporal lag is read from ``derived/MSIlluminanceToAGCLag.mat``. Each raw
``GKA`` recording is processed in memory, filtered for valid unsaturated
samples, and appended to the point cloud written to
``data/empircalAGCAndIlluminance.mat``.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat, savemat


THIS_FILE: Path = Path(__file__).resolve()
PROJECT_ROOT: Path = THIS_FILE.parents[3]
CALIBRATION_CODE_DIR: Path = THIS_FILE.parent.parent
VIDEO_IO_DIR: Path = (
    PROJECT_ROOT / "code" / "library" / "matlabIO" / "python_libraries"
)
DEFAULT_LAG_PATH: Path = (
    PROJECT_ROOT / "derived" / "MSIlluminanceToAGCLag.mat"
)
DEFAULT_OUTPUT_PATH: Path = (
    PROJECT_ROOT / "data" / "empircalAGCAndIlluminance.mat"
)
DEFAULT_MAXIMUM_SATURATION_PERCENT: float = 40.0
DEFAULT_INITIAL_SAMPLES_TO_EXCLUDE: int = 100


@dataclass(frozen=True)
class EmpiricalPointCloud:
    """Matched linear-scale calibration vectors selected for MATLAB."""

    camera_score: np.ndarray
    illuminance: np.ndarray


def _load_shared_lag_seconds(lag_path: Path) -> float:
    """Load and validate ``MSIlluminanceToAGCLag.sharedLagSeconds``."""

    if not lag_path.is_file():
        raise FileNotFoundError(
            "Run code/defineWorldCameraCalibration/"
            "defineMSIlluminanceToAGCLag.py first:\n"
            f"{lag_path}"
        )
    lag_file: dict[str, Any] = loadmat(lag_path, simplify_cells=True)
    lag_data: Any = lag_file.get("MSIlluminanceToAGCLag")
    if not isinstance(lag_data, dict) or "sharedLagSeconds" not in lag_data:
        raise KeyError(
            "The lag file must contain "
            "MSIlluminanceToAGCLag.sharedLagSeconds: "
            f"{lag_path}"
        )
    lag_seconds: float = float(lag_data["sharedLagSeconds"])
    if not np.isfinite(lag_seconds):
        raise ValueError(
            "MSIlluminanceToAGCLag.sharedLagSeconds must be finite."
        )
    return lag_seconds


def _processed_world_video_path(raw_recording_path: Path) -> Path:
    """Map a raw GKA path to its mirrored processed ``W.avi`` path."""

    path_parts: list[str] = list(raw_recording_path.parts)
    try:
        raw_root_index: int = path_parts.index("FLIC_raw")
    except ValueError as error:
        raise ValueError(
            "Raw recording paths must be beneath a directory named FLIC_raw."
        ) from error
    path_parts[raw_root_index] = "FLIC_processing"
    return Path(*path_parts) / "W.avi"


def _nearest_saturation_percent(
    sample_times: np.ndarray,
    frame_times: np.ndarray,
    frame_saturation: np.ndarray,
) -> np.ndarray:
    """Match every sample to the nearest decoded video frame."""

    if frame_times.size == 0:
        return np.full(sample_times.shape, np.nan, dtype=float)
    right_indices: np.ndarray = np.clip(
        np.searchsorted(frame_times, sample_times),
        0,
        frame_times.size - 1,
    )
    left_indices: np.ndarray = np.clip(
        right_indices - 1,
        0,
        frame_times.size - 1,
    )
    choose_right: np.ndarray = (
        np.abs(frame_times[right_indices] - sample_times)
        < np.abs(sample_times - frame_times[left_indices])
    )
    nearest_indices: np.ndarray = np.where(
        choose_right,
        right_indices,
        left_indices,
    )
    return 100.0 * frame_saturation[nearest_indices]


def _select_samples(
    camera_score: np.ndarray,
    illuminance: np.ndarray,
    saturation_percent: np.ndarray,
    maximum_saturation_percent: float,
    initial_samples_to_exclude: int,
) -> EmpiricalPointCloud | None:
    """Select valid calibration samples from one aligned recording."""

    if not (
        camera_score.shape == illuminance.shape == saturation_percent.shape
    ):
        raise ValueError("Aligned recording vectors must have matching shapes.")

    # The first samples are dropped after alignment, matching the historical
    # transient exclusion on each recording rather than on the pooled data.
    selected: np.ndarray = (
        np.isfinite(camera_score)
        & (camera_score > 0)
        & np.isfinite(illuminance)
        & (illuminance > 0)
        & np.isfinite(saturation_percent)
        & (saturation_percent <= maximum_saturation_percent)
    )
    selected[: min(initial_samples_to_exclude, selected.size)] = False
    if not np.any(selected):
        return None
    return EmpiricalPointCloud(
        camera_score=camera_score[selected],
        illuminance=illuminance[selected],
    )


def _derive_point_cloud(
    recording_paths: Sequence[Path],
    lag_seconds: float,
    maximum_saturation_percent: float,
    initial_samples_to_exclude: int,
) -> EmpiricalPointCloud:
    """Process raw recordings and build the selected point cloud in memory."""

    if not recording_paths:
        raise ValueError("At least one raw GKA recording path is required.")
    if not np.isfinite(maximum_saturation_percent) or not (
        0.0 <= maximum_saturation_percent <= 100.0
    ):
        raise ValueError("maximum_saturation_percent must be between 0 and 100.")
    if initial_samples_to_exclude < 0:
        raise ValueError("initial_samples_to_exclude must be nonnegative.")

    # Keep raw-data dependencies lazy so importing this module remains light.
    for import_path in (CALIBRATION_CODE_DIR, VIDEO_IO_DIR):
        import_path_string: str = os.fspath(import_path)
        if import_path_string not in sys.path:
            sys.path.insert(0, import_path_string)

    import defineMSIlluminanceToAGCLag
    import video_io

    selected_scores: list[np.ndarray] = []
    selected_illuminance: list[np.ndarray] = []
    processed_count: int = 0
    skipped_count: int = 0
    # Reuse one engine because starting MATLAB dominates short-recording cost.
    matlab_engine: Any = defineMSIlluminanceToAGCLag.start_matlab_engine()
    try:
        for recording_index, recording_path in enumerate(recording_paths, start=1):
            print(
                f"[{recording_index}/{len(recording_paths)}] "
                f"Processing {recording_path}"
            )
            try:
                recording: defineMSIlluminanceToAGCLag.RecordingSignals = (
                    defineMSIlluminanceToAGCLag.load_recording_signals(
                        recording_path,
                        matlab_engine=matlab_engine,
                    )
                )
            except defineMSIlluminanceToAGCLag.InvalidRecordingError as error:
                skipped_count += 1
                print(f"Skipping invalid recording: {error}", file=sys.stderr)
                continue

            # Apply the controller kernel and the single lag derived across all
            # recordings before putting both sensors on the minispect timebase.
            filtered_time, unshifted_illuminance = (
                defineMSIlluminanceToAGCLag.apply_empirical_kernel(
                    recording.minispect_time,
                    recording.illuminance,
                )
            )
            shifted_illuminance: np.ndarray = (
                defineMSIlluminanceToAGCLag.lag_filtered_signal(
                    filtered_time,
                    unshifted_illuminance,
                    lag_seconds,
                )
            )
            illuminance: np.ndarray = np.interp(
                recording.minispect_time,
                filtered_time,
                shifted_illuminance,
                left=np.nan,
                right=np.nan,
            )
            camera_score: np.ndarray = np.interp(
                recording.minispect_time,
                recording.camera_time,
                recording.camera_score,
                left=np.nan,
                right=np.nan,
            )
            paired: np.ndarray = np.isfinite(camera_score) & np.isfinite(illuminance)

            # Saturation is measured from the processed video, whose frames use
            # elapsed time relative to the first raw camera timestamp.
            video_path: Path = _processed_world_video_path(recording.path)
            if not video_path.is_file():
                raise FileNotFoundError(
                    f"Processed world video does not exist:\n{video_path}"
                )
            frame_saturation: np.ndarray = np.asarray(
                video_io.mesaure_video_saturation(
                    os.fspath(video_path),
                    verbose=True,
                ),
                dtype=float,
            )
            frame_rate: float = float(
                video_io.inspect_video_FPS(os.fspath(video_path))
            )
            if frame_rate <= 0:
                raise RuntimeError(f"Invalid video frame rate: {frame_rate}")
            frame_times: np.ndarray = (
                np.arange(frame_saturation.size, dtype=float) / frame_rate
            )
            sample_times: np.ndarray = (
                recording.minispect_time[paired] - recording.camera_time[0]
            )
            saturation_percent: np.ndarray = _nearest_saturation_percent(
                sample_times,
                frame_times,
                frame_saturation,
            )

            # Retain only selected vectors; no intermediate processing cache is
            # written or held after the final point cloud is concatenated.
            recording_points: EmpiricalPointCloud | None = _select_samples(
                camera_score[paired],
                illuminance[paired],
                saturation_percent,
                maximum_saturation_percent,
                initial_samples_to_exclude,
            )
            if recording_points is not None:
                selected_scores.append(recording_points.camera_score)
                selected_illuminance.append(recording_points.illuminance)
            processed_count += 1
    finally:
        matlab_engine.quit()

    if processed_count == 0:
        raise RuntimeError("No raw recording produced calibration data.")
    print(
        f"Processed {processed_count} recording(s); skipped {skipped_count} "
        f"invalid recording(s)."
    )

    if not selected_scores:
        raise RuntimeError("No AGC-to-illuminance samples passed selection.")
    return EmpiricalPointCloud(
        camera_score=np.concatenate(selected_scores),
        illuminance=np.concatenate(selected_illuminance),
    )


def derive_empircal_agc_and_illuminance(
    recording_paths: Sequence[str | os.PathLike[str]],
    *,
    lag_path: str | os.PathLike[str] = DEFAULT_LAG_PATH,
    output_path: str | os.PathLike[str] = DEFAULT_OUTPUT_PATH,
    maximum_saturation_percent: float = DEFAULT_MAXIMUM_SATURATION_PERCENT,
    initial_samples_to_exclude: int = DEFAULT_INITIAL_SAMPLES_TO_EXCLUDE,
) -> Path:
    """Process recordings, select the point cloud, and write MATLAB data."""

    resolved_lag_path: Path = Path(lag_path).expanduser().resolve()
    destination: Path = Path(output_path).expanduser().resolve()
    lag_seconds: float = _load_shared_lag_seconds(resolved_lag_path)

    raw_paths: list[Path] = [
        Path(path).expanduser().resolve() for path in recording_paths
    ]
    point_cloud: EmpiricalPointCloud = _derive_point_cloud(
        raw_paths,
        lag_seconds,
        maximum_saturation_percent,
        initial_samples_to_exclude,
    )
    # Atomic replacement keeps the checked-in calibration intact on failure.
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path = destination.with_suffix(".tmp.mat")
    savemat(
        temporary_path,
        {
            "empiralAGC": {
                "cameraScoreLinear": point_cloud.camera_score,
                "msIlluminance": point_cloud.illuminance,
                "sharedLagSeconds": lag_seconds,
            }
        },
        do_compression=True,
        oned_as="column",
    )
    os.replace(temporary_path, destination)
    return destination


def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "recording_paths",
        nargs="+",
        help="Raw GKA recording directories used to build the point cloud.",
    )
    parser.add_argument("--lag-path", type=Path, default=DEFAULT_LAG_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument(
        "--maximum-saturation-percent",
        type=float,
        default=DEFAULT_MAXIMUM_SATURATION_PERCENT,
    )
    parser.add_argument(
        "--initial-samples-to-exclude",
        type=int,
        default=DEFAULT_INITIAL_SAMPLES_TO_EXCLUDE,
    )
    return parser


def main(argv: Sequence[str] | None = None) -> Path:
    """Run the export and return the generated MATLAB path."""

    arguments = _build_argument_parser().parse_args(argv)
    output_path: Path = derive_empircal_agc_and_illuminance(
        arguments.recording_paths,
        lag_path=arguments.lag_path,
        output_path=arguments.output_path,
        maximum_saturation_percent=arguments.maximum_saturation_percent,
        initial_samples_to_exclude=arguments.initial_samples_to_exclude,
    )
    print(f"Saved empirical AGC-to-illuminance data to {output_path}")
    return output_path


if __name__ == "__main__":
    main()
