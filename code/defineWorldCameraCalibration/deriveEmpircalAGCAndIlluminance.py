"""Derive the empirical camera-AGC/minispect-illuminance calibration data.

Running this file without arguments rebuilds ``data/empircalAGC.mat`` from
the checked-in processed-data cache. Optional raw ``GKA`` recording paths may
be supplied to rebuild that cache before the MATLAB file is exported.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from scipy.io import savemat


THIS_FILE: Path = Path(__file__).resolve()
LIGHT_LOGGER_ANALYSIS_ROOT: Path = THIS_FILE.parents[2]
FIT_AGC_TO_ILLUMINANCE_DIR: Path = (
    LIGHT_LOGGER_ANALYSIS_ROOT / "code" / "fitAGCtoIlluminance"
)
DEFAULT_CACHE_PATH: Path = (
    FIT_AGC_TO_ILLUMINANCE_DIR
    / "cached_processing_data"
    / "frame_saturation_calibration_data.npz"
)
DEFAULT_OUTPUT_PATH: Path = (
    LIGHT_LOGGER_ANALYSIS_ROOT / "data" / "empircalAGC.mat"
)

def _rebuild_processed_cache(
    recording_paths: Sequence[str],
    cache_path: Path,
    illuminance_diagnostics: bool,
) -> None:
    """Run the raw-recording workflow only when paths were supplied."""

    fit_module_path: str = os.fspath(FIT_AGC_TO_ILLUMINANCE_DIR)
    if fit_module_path not in sys.path:
        sys.path.insert(0, fit_module_path)

    # The raw workflow has additional recording/video dependencies. Keeping
    # this import local lets the normal cache-only derivation run with just
    # NumPy and SciPy.
    from fit_agc_to_illuminance_util import fit_agc_to_illuminance

    fit_agc_to_illuminance(
        list(recording_paths),
        illuminance_diagnostics=illuminance_diagnostics,
        frame_saturation_data_output_path=cache_path,
    )


def _selected_linear_point_cloud(
    cache_path: Path,
    maximum_saturation_percent: float | None,
    initial_samples_to_exclude: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Load and select the matched linear-scale vectors from the cache."""

    if (
        maximum_saturation_percent is not None
        and (
            not np.isfinite(maximum_saturation_percent)
            or not 0.0 <= maximum_saturation_percent <= 100.0
        )
    ):
        raise ValueError("maximum_saturation_percent must be between 0 and 100.")
    if initial_samples_to_exclude < 0:
        raise ValueError("initial_samples_to_exclude must be nonnegative.")

    selected_camera_scores: list[np.ndarray] = []
    selected_ms_illuminance: list[np.ndarray] = []

    with np.load(cache_path, allow_pickle=False) as archive:
        format_version: int = int(archive["format_version"][0])
        if format_version not in (1, 2):
            raise ValueError(
                f"Unsupported frame-saturation cache version: {format_version}."
            )

        recording_count: int = int(archive["recording_count"][0])
        for recording_index in range(recording_count):
            key_suffix: str = str(recording_index)
            cached_camera_score: np.ndarray = np.asarray(
                archive[f"camera_score_{key_suffix}"],
                dtype=float,
            )
            if format_version == 1:
                camera_score: np.ndarray = np.full(
                    cached_camera_score.shape,
                    np.nan,
                    dtype=float,
                )
                valid_reciprocal: np.ndarray = (
                    np.isfinite(cached_camera_score) & (cached_camera_score > 0)
                )
                camera_score[valid_reciprocal] = (
                    1.0 / cached_camera_score[valid_reciprocal]
                )
            else:
                camera_score = cached_camera_score

            ms_illuminance: np.ndarray = np.asarray(
                archive[f"filtered_illuminance_{key_suffix}"],
                dtype=float,
            )
            saturation_percent: np.ndarray = np.asarray(
                archive[f"frame_saturation_percent_{key_suffix}"],
                dtype=float,
            )
            if not (
                camera_score.shape
                == ms_illuminance.shape
                == saturation_percent.shape
            ):
                raise ValueError(
                    "Cached camera-score, illuminance, and saturation vectors "
                    f"do not match for recording index {recording_index}."
                )

            camera_score_derivative: np.ndarray
            if camera_score.size == 0:
                camera_score_derivative = np.array([], dtype=float)
            elif camera_score.size == 1:
                camera_score_derivative = np.zeros(1, dtype=float)
            else:
                camera_score_derivative = np.gradient(camera_score)

            selected: np.ndarray = (
                np.isfinite(camera_score)
                & (camera_score > 0)
                & np.isfinite(ms_illuminance)
                & (ms_illuminance > 0)
                & np.isfinite(saturation_percent)
                & np.isfinite(camera_score_derivative)
            )
            if maximum_saturation_percent is not None:
                selected &= saturation_percent <= maximum_saturation_percent

            excluded_count: int = min(
                initial_samples_to_exclude,
                selected.size,
            )
            selected[:excluded_count] = False
            if np.any(selected):
                selected_camera_scores.append(camera_score[selected])
                selected_ms_illuminance.append(ms_illuminance[selected])

    if not selected_camera_scores:
        raise RuntimeError("No AGC-to-illuminance samples passed the selection.")

    return (
        np.concatenate(selected_camera_scores),
        np.concatenate(selected_ms_illuminance),
    )


def derive_empircal_agc_and_illuminance(
    recording_paths: Sequence[str | os.PathLike[str]] = (),
    *,
    processed_data_path: str | os.PathLike[str] = DEFAULT_CACHE_PATH,
    output_path: str | os.PathLike[str] = DEFAULT_OUTPUT_PATH,
    maximum_saturation_percent: float | None = 40.0,
    initial_samples_to_exclude: int = 100,
    illuminance_diagnostics: bool = False,
) -> Path:
    """Create ``empircalAGC.mat`` from raw recordings or the existing cache.

    The MATLAB file contains one scalar struct named ``empiralAGC`` with the
    fields ``cameraScoreLinear`` and ``msIlluminance``. Both fields contain
    the selected linear-scale point cloud in matched column vectors.

    Args:
        recording_paths: Optional raw ``GKA`` recording directories. When
            supplied, the processed-data cache is rebuilt before export.
        processed_data_path: Cached point-cloud data used for the final point
            selection and MATLAB export.
        output_path: Destination MATLAB file.
        maximum_saturation_percent: Maximum retained frame saturation.
        initial_samples_to_exclude: Samples removed from the beginning of
            every recording before export.
        illuminance_diagnostics: Whether raw processing should generate the
            optional eight-versus-nine-channel diagnostic products.

    Returns:
        Absolute path to the generated MATLAB file.
    """

    cache_path: Path = Path(processed_data_path).expanduser().resolve()
    destination: Path = Path(output_path).expanduser().resolve()
    resolved_recording_paths: list[str] = [
        os.fspath(Path(path).expanduser().resolve()) for path in recording_paths
    ]

    if resolved_recording_paths:
        _rebuild_processed_cache(
            resolved_recording_paths,
            cache_path,
            illuminance_diagnostics,
        )
    elif not cache_path.is_file():
        raise FileNotFoundError(
            "The processed AGC-to-illuminance cache does not exist. Supply "
            f"one or more raw GKA recording paths or restore the cache:\n{cache_path}"
        )

    camera_score_linear, ms_illuminance = _selected_linear_point_cloud(
        cache_path,
        maximum_saturation_percent,
        initial_samples_to_exclude,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    savemat(
        destination,
        {
            "empiralAGC": {
                "cameraScoreLinear": camera_score_linear,
                "msIlluminance": ms_illuminance,
            }
        },
        do_compression=True,
        oned_as="column",
    )

    if not destination.is_file():
        raise RuntimeError(f"The expected MATLAB file was not created: {destination}")

    return destination


def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for the derivation script."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "recording_paths",
        nargs="*",
        help="Optional raw GKA recording directories used to rebuild the cache.",
    )
    parser.add_argument(
        "--processed-data-path",
        type=Path,
        default=DEFAULT_CACHE_PATH,
        help="Processed-data cache used for the final export.",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="MATLAB output path (default: data/empircalAGC.mat).",
    )
    parser.add_argument(
        "--maximum-saturation-percent",
        type=float,
        default=40.0,
        help="Maximum retained frame saturation percentage.",
    )
    parser.add_argument(
        "--initial-samples-to-exclude",
        type=int,
        default=100,
        help="Samples excluded from the start of each recording.",
    )
    parser.add_argument(
        "--illuminance-diagnostics",
        action="store_true",
        help="Generate optional channel-count diagnostics during raw processing.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> Path:
    """Run the derivation and return the generated MATLAB file path."""

    arguments = _build_argument_parser().parse_args(argv)
    output_path: Path = derive_empircal_agc_and_illuminance(
        arguments.recording_paths,
        processed_data_path=arguments.processed_data_path,
        output_path=arguments.output_path,
        maximum_saturation_percent=arguments.maximum_saturation_percent,
        initial_samples_to_exclude=arguments.initial_samples_to_exclude,
        illuminance_diagnostics=arguments.illuminance_diagnostics,
    )
    print(f"Saved empirical AGC-to-illuminance data to {output_path}")
    return output_path


if __name__ == "__main__":
    main()
