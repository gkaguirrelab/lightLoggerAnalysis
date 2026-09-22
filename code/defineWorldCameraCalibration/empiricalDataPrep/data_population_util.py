"""Shared helpers for materializing synchronized world-camera MAT files."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import chunk_io
import numpy as np
import pandas as pd
import world_util
from scipy.io import savemat


DEFAULT_MAT_README: str = (
    "Synchronized context for one selected raw world-camera frame. "
    "sourceImageFilename is the frame's descriptive label; "
    "globalWorldFrameIndex is its zero-based index across naturally ordered "
    "raw world chunks, without synthetic gap-filled frames; "
    "worldTimestampSeconds is on the shared logger clock; worldFrame is the "
    "unprocessed Bayer frame; AGCSettings contains the Again, Dgain, and "
    "exposure values used by the reconstruction pipeline; "
    "minispectTimestampSeconds and minispectValue contain the nearest "
    "minispect packet and its parsed AS, TS, LS, and TEMP values."
)

LEGACY_AGC_SETTING_NAMES: tuple[str, str, str] = (
    "Again",
    "Dgain",
    "exposure",
)
MODERN_AGC_SETTING_NAMES: tuple[str, str, str] = (
    "cameraAgain",
    "AGCDgain",
    "cameraExposure",
)


def _world_timestamp_at_global_frame_index(
    raw_chunks: Path, global_frame_index: int
) -> float:
    """Return the timestamp of one physically stored raw world frame."""
    # Use the established world-camera metadata loader so natural chunk
    # ordering, legacy schemas, and nanosecond-to-second conversion remain
    # centralized in world_util.
    world_metadata: pd.DataFrame = world_util.world_metadata_from_chunks(
        str(raw_chunks), convert_to_seconds=True, verbose=False
    )

    setting_columns: list[str] = [
        column for column in world_metadata.columns if column != "timestamp"
    ]
    if not setting_columns:
        raise ValueError("World metadata contains no camera-setting columns.")

    # world_metadata_from_chunks represents missing frames by inserting rows
    # with a valid interpolated timestamp and NaN in every camera-setting field.
    # Remove exactly those rows so iloc once again addresses captured frames
    # only—the same zero-based index convention used in the data READMEs.
    gap_filled_rows: pd.Series = world_metadata[setting_columns].isna().all(axis=1)
    captured_metadata: pd.DataFrame = world_metadata.loc[
        ~gap_filled_rows
    ].reset_index(drop=True)

    if global_frame_index >= len(captured_metadata):
        raise IndexError(
            f"Global frame index {global_frame_index} exceeds "
            f"{len(captured_metadata)} physically stored raw frames."
        )

    timestamp_seconds: float = float(
        captured_metadata.iloc[global_frame_index]["timestamp"]
    )
    if not np.isfinite(timestamp_seconds):
        raise ValueError(
            f"Captured frame {global_frame_index} has a non-finite timestamp."
        )
    return timestamp_seconds


def _normalize_agc_settings(
    agc_settings: Mapping[str, float],
) -> dict[str, float]:
    """Normalize legacy or modern metadata to Again/Dgain/exposure."""
    # Older recordings already use the three MATLAB-facing names. Preserve
    # those values exactly when that schema is available.
    if all(name in agc_settings for name in LEGACY_AGC_SETTING_NAMES):
        source_names: tuple[str, str, str] = LEGACY_AGC_SETTING_NAMES

    # Newer recordings distinguish camera-applied settings from the values
    # requested by the AGC algorithm. The applied values produced this frame,
    # so those are the correct values to serialize for reconstruction.
    elif all(name in agc_settings for name in MODERN_AGC_SETTING_NAMES):
        source_names = MODERN_AGC_SETTING_NAMES
    else:
        raise ValueError(
            "World metadata carries neither the legacy settings "
            f"{LEGACY_AGC_SETTING_NAMES} nor the modern settings "
            f"{MODERN_AGC_SETTING_NAMES}; found {sorted(agc_settings)}."
        )

    # Always expose one stable schema to MATLAB regardless of recording age.
    return {
        target_name: float(agc_settings[source_name])
        for target_name, source_name in zip(
            LEGACY_AGC_SETTING_NAMES, source_names, strict=True
        )
    }


def find_target_frame(
    raw_chunks: str | Path,
    global_frame_index: int,
    output_path: str | Path,
    *,
    source_image_filename: str,
    readme: str = DEFAULT_MAT_README,
    overwrite: bool = False,
) -> Path:
    """Write one synchronized raw world-frame measurement to a MAT file.

    Parameters
    ----------
    raw_chunks
        Directory containing the original world-camera and minispect chunks.
    global_frame_index
        Zero-based index across the physically stored world frames in natural
        chunk order. Missing frames are not synthesized or counted.
    output_path
        Destination ``.mat`` path.
    source_image_filename
        Descriptive historical frame label stored in the MAT file.
    readme
        Human-readable description stored in the MAT file's ``README`` field.
    overwrite
        Replace an existing output only when explicitly enabled.

    Returns
    -------
    pathlib.Path
        The resolved output path, whether newly written or preserved.
    """
    # Normalize path-like inputs once so validation, creation, and the returned
    # value all refer to the same absolute locations.
    raw_chunks = Path(raw_chunks).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()

    # Reject malformed inputs before opening any potentially large chunk files.
    if not raw_chunks.is_dir():
        raise FileNotFoundError(
            f"Raw world-camera chunk directory does not exist: {raw_chunks}"
        )
    if (
        isinstance(global_frame_index, bool)
        or not isinstance(global_frame_index, (int, np.integer))
        or global_frame_index < 0
    ):
        raise ValueError("global_frame_index must be a nonnegative integer.")
    if output_path.suffix.lower() != ".mat":
        raise ValueError(f"Output path must end in .mat: {output_path}")
    if not source_image_filename.strip():
        raise ValueError("source_image_filename must not be empty.")

    # Preservation is the default. Callers must opt in explicitly before an
    # existing scientific data product can be replaced.
    if output_path.exists() and not overwrite:
        return output_path

    # Translate the physical global index to the selected frame's exact clock
    # time. This lookup never constructs a video or inserts missing frames.
    timestamp_seconds: float = _world_timestamp_at_global_frame_index(
        raw_chunks, int(global_frame_index)
    )

    # The exact world timestamp lets chunk_io retrieve that raw frame together
    # with its AGC metadata and the temporally nearest minispect packet.
    context: dict[str, dict[str, Any]] = chunk_io.find_nearest_neighbor(
        str(raw_chunks), timestamp_seconds, sensors=("W", "M")
    )

    # Guard against accidentally pairing the requested index with an adjacent
    # frame if timestamp conversion or source metadata ever changes.
    if not np.isclose(
        float(context["W"]["timestamp"]), timestamp_seconds, rtol=0.0, atol=1e-9
    ):
        raise ValueError(
            "The indexed world timestamp did not resolve to the same raw frame."
        )

    # Create only the requested destination directory; source data remains
    # read-only throughout this operation.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Match the established exampleWorldCameraImages schema exactly so the
    # same MATLAB readers can consume either the example or Macbeth outputs.
    savemat(
        output_path,
        {
            "README": readme,
            "sourceImageFilename": source_image_filename,
            "globalWorldFrameIndex": np.int64(global_frame_index),
            "worldTimestampSeconds": context["W"]["timestamp"],
            "worldFrame": context["W"]["value"],
            "AGCSettings": _normalize_agc_settings(
                context["W"]["AGCSettings"]
            ),
            "minispectTimestampSeconds": context["M"]["timestamp"],
            "minispectValue": context["M"]["value"],
        },
        do_compression=True,
        long_field_names=True,
    )
    return output_path
