"""Shared helpers for locating synchronized samples in raw recordings."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import chunk_io
import numpy as np
import world_util


# -----------------------------------------------------------------------------
# Raw recording frame counts
# Count only arrays physically stored in naturally ordered chunk files. This
# intentionally excludes any synthetic rows inserted later for missing frames.
# -----------------------------------------------------------------------------
def count_num_frames(
    recording_path: str,
) -> dict[Literal["W", "M"], int]:
    """Count physically stored frames for each supported recording sensor.

    Parameters
    ----------
    recording_path
        Directory containing the raw world-camera and/or minispect chunk files.

    Returns
    -------
    dict
        Frame counts keyed by ``"W"`` for the world camera and ``"M"`` for
        the minispectrometer. A sensor that has no chunks in the recording is
        included with a count of zero.
    """
    chunks_by_sensor: dict[str, list[tuple[str, str]]] = (
        chunk_io.group_sensors_files(recording_path)
    )

    frame_count_by_sensor: dict[Literal["W", "M"], int] = {
        sensor: sum(
            int(np.load(frame_path, mmap_mode="r").shape[0])
            for _, frame_path in chunks_by_sensor.get(sensor, [])
        )
        for sensor in ("W", "M")
    }
    return frame_count_by_sensor


def find_target_frame(
    raw_chunks: str | Path,
    global_frame_index: int,
) -> dict[str, Any]:
    """Return one raw world frame and its synchronized sensor context.

    Parameters
    ----------
    raw_chunks
        Directory containing the original world-camera and minispect chunks.
    global_frame_index
        Zero-based index across the physically stored world frames in natural
        chunk order. Missing frames are not synthesized or counted.
    Returns
    -------
    dict
        MATLAB-ready context containing ``worldTimestampSeconds``,
        ``worldFrame``, normalized ``AGCSettings``,
        ``minispectTimestampSeconds``, and ``minispectValue``. The caller is
        responsible for adding provenance fields and saving any output file.
    """
    # Resolve the recording path without changing or writing its contents.
    raw_chunks: Path = Path(raw_chunks).expanduser().resolve()

    # Reject invalid paths and indices before loading recording metadata.
    if not raw_chunks.is_dir():
        raise FileNotFoundError(
            f"Raw recording chunk directory does not exist: {raw_chunks}"
        )
    if (
        isinstance(global_frame_index, bool)
        or not isinstance(global_frame_index, (int, np.integer))
        or global_frame_index < 0
    ):
        raise ValueError("global_frame_index must be a nonnegative integer.")

    # Load the naturally ordered world metadata with timestamps in seconds.
    world_metadata = world_util.world_metadata_from_chunks(
        str(raw_chunks), convert_to_seconds=True, verbose=False
    )
    # Synthetic gap rows contain NaN camera settings. Remove every row that
    # contains a NaN so iloc indexes only frames physically present on disk.
    captured_metadata = world_metadata.dropna().reset_index(drop=True)
    if global_frame_index >= len(captured_metadata):
        raise IndexError(
            f"Global frame index {global_frame_index} exceeds "
            f"{len(captured_metadata)} physically stored raw frames."
        )
    # Translate the physical global index to its exact shared-clock timestamp.
    timestamp_seconds: float = float(
        captured_metadata.iloc[int(global_frame_index)]["timestamp"]
    )

    # Retrieve that world frame and the temporally nearest minispect sample.
    context: dict[str, dict[str, Any]] = chunk_io.find_nearest_neighbor(
        str(raw_chunks), timestamp_seconds, sensors=("W", "M")
    )

    # Normalize legacy and modern camera metadata to the three fields consumed
    # by the MATLAB reconstruction pipeline (AGain, DGain, exposure)
    agc_settings: dict[str, Any] = context["W"]["AGCSettings"]
    if all(name in agc_settings for name in ("Again", "Dgain", "exposure")):
        source_names: tuple[str, str, str] = ("Again", "Dgain", "exposure")
    elif all(
        name in agc_settings
        for name in ("cameraAgain", "AGCDgain", "cameraExposure")
    ):
        source_names = ("cameraAgain", "AGCDgain", "cameraExposure")
    else:
        raise ValueError(
            "World metadata does not contain a recognized AGC-setting schema."
        )

    normalized_agc_settings: dict[str, float] = {
        target_name: float(agc_settings[source_name])
        for target_name, source_name in zip(
            ("Again", "Dgain", "exposure"), source_names, strict=True
        )
    }

    # Return only synchronized frame context; output policy belongs to callers.
    return {
        "worldTimestampSeconds": context["W"]["timestamp"],
        "worldFrame": context["W"]["value"],
        "AGCSettings": normalized_agc_settings,
        "minispectTimestampSeconds": context["M"]["timestamp"],
        "minispectValue": context["M"]["value"],
    }
