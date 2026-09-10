"""Utilities that operate on complete light-logger recordings."""

import os
import sys
from typing import Literal
from collections.abc import Callable, Iterable

import hdf5storage
import numpy as np
from natsort import natsorted
from tqdm.auto import tqdm

sys.path.append(os.path.dirname(__file__))
import ms_util
import world_util


def group_sensors_files(recording_path: str) -> dict[str, list[tuple[str, str]]]:
    """Return naturally sorted metadata/data chunk pairs for each sensor.

    Args:
        recording_path: Directory containing the raw sensor chunk files.

    Returns:
        A dictionary keyed by ``"W"``, ``"P"``, and ``"M"``. Each value is a
        list of ``(metadata_path, data_path)`` tuples in natural chunk order.
    """

    def group_sensor_files(sensor_name: str) -> list[tuple[str, str]]:
        """Return metadata/data chunk pairs for one sensor.

        Args:
            sensor_name: Sensor name used to identify matching filenames.

        Returns:
            A naturally sorted list of ``(metadata_path, data_path)`` tuples.
        """
        # Find each sensor's metadata file and pair it with the data file that
        # has the same name without the "_metadata" suffix.
        return [
            (
                os.path.join(recording_path, filename),
                os.path.join(recording_path, filename.replace("_metadata", "")),
            )
            for filename in natsorted(os.listdir(recording_path))
            if sensor_name in filename and "metadata" in filename
        ]

    # Group all three sensor types by their single-letter sensor identifier.
    return {
        sensor[0].upper(): group_sensor_files(sensor)
        for sensor in ("world", "pupil", "ms")
    }


def _process_raw_world_helper(
    files: list[tuple[str, str]],
    output_path: str,
    chunk_range: tuple[int | None]=(0, None), 
    overwrite_existing: bool = False,
    verbose: bool = False,
) -> None:
    """Process world-camera chunks and save them as MATLAB files.

    Args:
        files: Metadata/data path pairs for the world-camera chunks.
        output_path: Directory in which to save the processed chunks.
        overwrite_existing: Whether to replace existing output files.
        verbose: Whether to display a processing progress bar.

    Returns:
        None. Processed chunks are written to ``output_path``.

    Raises:
        AssertionError: If a world metadata array does not use the supported
            legacy or modern column layout.
    """
    # Get the start and end chunk we wish to process 
    # while replacing the None sentiel if it exists
    start_chunk, end_chunk = chunk_range
    end_chunk = len(files) if end_chunk is None else end_chunk
    desired_chunk_range: object = range(start_chunk, end_chunk)

    # Ensure the world-camera output directory exists.
    os.makedirs(output_path, exist_ok=True)

    # Add a progress bar to the chunk iterator when verbose output is enabled.
    file_iterator: Iterable[tuple[int, tuple[str, str]]] = enumerate(files)
    if verbose:
        file_iterator = tqdm(
            file_iterator,
            total=len(files),
            desc="Processing world chunks",
        )

    # Iterate over the paired metadata and frame-buffer files.
    for file_num, (metadata_path, data_path) in file_iterator:
        # Skip files not in the desired range 
        if(file_num not in desired_chunk_range):
            continue

        # Form the output path and skip an existing file unless overwrite was
        # requested.
        output_filepath: str = os.path.join(
            output_path,
            f"world_chunk{file_num}.mat",
        )
        if os.path.exists(output_filepath) and not overwrite_existing:
            continue

        # Load the metadata and its associated frame buffer.
        metadata_buffer: np.ndarray = np.load(metadata_path)
        frame_buffer: np.ndarray = np.load(data_path)

        # The final chunk of a recording may be empty.
        if len(frame_buffer) == 0:
            continue

        # The metadata is either legacy shaped or modern shaped
        # Legacy shape is timestamp, Again, DGain, Exposure
        # Modern shape is timestamp, "cameraAgain", "AGCDgain", "cameraExposure", "AGCAgain", "AGCExposure"
        Again_idx: int
        Dgain_idx: int
        exposure_idx: int

        # First check to see if the metadata is properly shaped
        assert metadata_buffer.shape[1] in (4, 6), f"Metadata buffer must have cols: (timestamp, Again, DGain, Exposure) or timestamp, cameraAgain, AGCDgain, cameraExposure, AGCAgain, AGCExposure"

        if(metadata_buffer.shape[1] == 4):
            Again_idx = 1
            Dgain_idx = 2
            exposure_idx = 3
        else:
            # Add the timestamp column to the modern AGC column names so the
            # resulting indices match the complete metadata buffer.
            modern_metadata_columns: tuple[str, ...] = (
                "timestamp",
                *world_util.WORLD_AGC_METADATA_COLS,
            )
            Again_idx = modern_metadata_columns.index("cameraAgain")
            Dgain_idx = modern_metadata_columns.index("AGCDgain")
            exposure_idx = modern_metadata_columns.index("cameraExposure")

        # Repackage the AGC settings specifically into the shape required for the pipeline
        agc_settings: dict[str, np.ndarray] = {
            "Again": metadata_buffer[:, Again_idx],
            "Dgain": metadata_buffer[:, Dgain_idx],
            "exposure": metadata_buffer[:, exposure_idx],
        }

        # Transform the raw frames and combine them with their named metadata.
        data_dict: dict[str, object] = {
            "data": world_util.world_transformation_pipeline(
                frame_buffer,
                agc_settings,
            ),
            "metadata": agc_settings | {"timestamps": metadata_buffer[:, 0]},
        }
        # Save as MATLAB v7.3 so large world-camera arrays are not limited by
        # the 2 GB matrix limit of the older MATLAB v5 format.
        hdf5storage.savemat(
            output_filepath,
            data_dict,
            fmt="7.3",
            store_python_metadata=False,
            truncate_existing=True,
        )


def _process_raw_ms_helper(
    files: list[tuple[str, str]],
    output_path: str,
    chunk_range: tuple[int | None]=(0, None),
    overwrite_existing: bool = False,
    verbose: bool = False,
) -> None:
    """Process minispectrometer chunks and save them as MATLAB files.

    Args:
        files: Metadata/data path pairs for the minispectrometer chunks.
        output_path: Directory in which to save the processed chunks.
        overwrite_existing: Whether to replace existing output files.
        verbose: Whether to display a processing progress bar.

    Returns:
        None. Processed chunks are written to ``output_path``.
    """

    # Get the start and end chunk we wish to process 
    # while replacing the None sentiel if it exists
    start_chunk, end_chunk = chunk_range
    end_chunk = len(files) if end_chunk is None else end_chunk
    desired_chunk_range: object = range(start_chunk, end_chunk)

    # Ensure the minispectrometer output directory exists.
    os.makedirs(output_path, exist_ok=True)

    # Add a progress bar to the chunk iterator when verbose output is enabled.
    file_iterator: Iterable[tuple[int, tuple[str, str]]] = enumerate(files)
    if verbose:
        file_iterator = tqdm(
            file_iterator,
            total=len(files),
            desc="Processing MS chunks",
        )

    # Iterate over the paired metadata and reading-buffer files.
    for file_num, (metadata_path, data_path) in file_iterator:
        # Skip files not in the desired range 
        if(file_num not in desired_chunk_range):
            continue

        # Form the output path and skip an existing file unless overwrite was
        # requested.
        output_filepath: str = os.path.join(
            output_path,
            f"ms_chunk{file_num}.mat",
        )
        if os.path.exists(output_filepath) and not overwrite_existing:
            continue

        # Load the timestamps and their associated minispectrometer readings.
        timestamps: np.ndarray = np.load(metadata_path).flatten()
        frame_buffer: np.ndarray = np.load(data_path)

        # The final chunk of a recording may be empty.
        if len(frame_buffer) == 0:
            continue

        # Estimate the radiance spectrum for each minispectrometer reading.
        radiance: np.ndarray
        sampling: np.ndarray
        fit_value: float | np.ndarray
        fit_errors: np.ndarray
        radiance, sampling, fit_value, fit_errors = (
            ms_util.estimate_radiance_spectrum_form_ms(
                frame_buffer,
                visualize_results=False,
            )
        )
        # Combine the reconstructed data and timestamps into named structures.
        data_dict: dict[str, object] = {
            "data": {
                "radiance": radiance,
                "S": sampling,
                "fVal": fit_value,
                "fitErrors": fit_errors,
            },
            "metadata": {"timestamps": timestamps},
        }
        # Save as MATLAB v7.3 to use the same output format as the world chunks.
        hdf5storage.savemat(
            output_filepath,
            data_dict,
            fmt="7.3",
            store_python_metadata=False,
            truncate_existing=True,
        )


def process_raw_recording(
    path_to_raw: str,
    output_path: str,
    overwrite_existing: bool = False,
    verbose: bool = False,
    chunk_ranges: dict[Literal["W", "M"], tuple[int | None]] = {sensor_name: (0, None) for sensor_name in "WM"}
) -> None:
    """Process all world-camera and minispectrometer chunks in a recording.

    Args:
        path_to_raw: Existing, nonempty directory containing raw sensor chunks.
        output_path: Destination directory for the processed sensor folders.
        overwrite_existing: Whether to replace existing processed chunks.
        verbose: Whether to display progress bars during processing.

    Returns:
        None. World and minispectrometer results are written beneath
        ``output_path`` in ``W`` and ``M`` subdirectories.

    Raises:
        AssertionError: If ``path_to_raw`` is not an existing, nonempty
            directory.
        AssertionError: If world metadata has an unsupported column layout.
    """
    # Ensure the raw recording is an existing, non-empty directory.
    assert (
        os.path.isdir(path_to_raw) and len(os.listdir(path_to_raw)) > 0
    ), f"{path_to_raw} must be an existing, non-empty directory"

    # Find the raw chunk files for each sensor.
    files_by_sensor: dict[str, list[tuple[str, str]]] = group_sensors_files(
        path_to_raw
    )

    # Map each supported sensor to its processing function.
    helper_map: dict[
        str,
        Callable[[list[tuple[str, str]], str, bool, bool], None],
    ] = {
        "W": _process_raw_world_helper,
        "M": _process_raw_ms_helper,
    }

    # Create the output directory if it does not already exist.
    os.makedirs(output_path, exist_ok=True)

    # Process each supported sensor into its own output subdirectory.
    for sensor_name, helper in helper_map.items():
        # Skip this because it was from a legacy implementation. 
        # This does not exist in any of our recordings now
        if(sensor_name == "P"):
            continue

        helper(
            files_by_sensor[sensor_name],
            os.path.join(output_path, sensor_name),
            chunk_ranges[sensor_name], 
            overwrite_existing,
            verbose,
        )


def main() -> None:
    """Run the module as a script.

    Returns:
        None.
    """
    pass


if __name__ == "__main__":
    main()
