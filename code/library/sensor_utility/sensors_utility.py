"""Utilities that operate on complete light-logger recordings."""

import os
import shutil
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from collections.abc import Iterable

import dill
import scipy.io
import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D
from matplotlib.axes import Axes
from matplotlib.axes._secondary_axes import SecondaryAxis
from matplotlib.collections import QuadMesh
from matplotlib.figure import Figure
from PIL import Image
from natsort import natsorted
from tqdm.auto import tqdm

sys.path.append(os.path.dirname(__file__))
import ms_util
import world_util

# Resolve the video helpers independently of the caller's working directory.
video_library: Path = Path(__file__).resolve().parents[1] / "matlabIO" / "python_libraries"
if str(video_library) not in sys.path:
    sys.path.append(str(video_library))

import video_io

# Raw world timestamps use nanoseconds; MS timestamps use seconds.
WORLD_NANOSECONDS_TO_SECONDS: float = 1e-9


def group_sensors_files(recording_path: str | Path) -> dict[str, list[tuple[str, str]]]:
    """Return naturally sorted raw metadata/data path pairs for each sensor.

    Args:
        recording_path: String or Path. Directory to scan for world, pupil, and ms filenames.

    Returns:
        Dictionary with W, P, and M keys, each containing a list of
        (metadata_path, data_path) pairs. Missing sensors have empty lists.
        Data names are derived by removing "_metadata" from metadata names;
        this function does not check that the corresponding data files exist.

    Raises:
        OSError: If the recording directory cannot be listed.
    """

    def group_sensor_files(sensor_name: str) -> list[tuple[str, str]]:
        """Find metadata files for one sensor in the enclosing recording directory.

        Args:
            sensor_name: Filename substring to match, such as "world" or "ms".

        Returns:
            Naturally sorted (metadata_path, data_path) pairs. Matching names
            must contain both sensor_name and "metadata".
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
    recording_metadata: dict,
    files: list[tuple[str | Path, str | Path]],
    output_path: str | Path,
    age: int | None = None,
    gaze_angles_xy: np.ndarray | None = None,
    chunk_range: tuple[int | None]=(0, None),
    overwrite_existing: bool = False,
    verbose: bool = False,
    n_workers: int = world_util.WORLD_IMPUTATION_WORKERS,
    demosaic_n_workers: int = world_util.WORLD_DEMOSAIC_WORKERS,
) -> None:
    """Transform selected world chunks into one compressed native HDF5 file.

    Args:
        recording_metadata: Configuration from config.pkl; requires a finite positive world FPS
            at ["sensors"]["W"]["sensor_mode"]["fps"].
        files: Pairs of string or Path filenames. Ordered (metadata .npy path, frame-buffer .npy path) pairs.
        output_path: String or Path. HDF5 filename, not a directory. Parent directories are created.
        age: Optional participant age forwarded to the transformation, currently unused.
        gaze_angles_xy: Camera pixel X/Y pairs aligned to all input frames, or a single
            pair shared by every frame. Coordinates are one-based, as in MATLAB.
        chunk_range: Zero-based (start, stop) chunk indices; stop is exclusive.
            None means the beginning/end. Defaults to all chunks. Empty buffers
            are skipped, and selected frames are concatenated in input order.
        overwrite_existing: Replace the entire output file when True. Otherwise,
            skip a completed file and raise for an incomplete file; no resuming.
        verbose: Show progress over input chunk pairs.
        n_workers: Imputation worker count, default WORLD_IMPUTATION_WORKERS (6).
        demosaic_n_workers: Independent demosaicing worker count, default
            WORLD_DEMOSAIC_WORKERS (16). Workers are reused across chunks.
            Set both worker counts to 1 for fully serial processing.

    Returns:
        None. Writes /data/frames as float64 (frames, azimuth, elevation,
        post-receptoral channel), preserving signed fractions and NaNs; /data/average_frames retains floating-point means,
        and /agcSettings/{Again,Dgain,exposure,timestamps} as aligned per-frame
        arrays. Timestamps retain raw nanosecond units. These datasets use
        lossless gzip level 4 with shuffle and automatic HDF5 storage chunks.
        /data/average_frames stores one spatial RGB mean per frame (N, 3),
        aligned with /agcSettings/timestamps and using the same compression.
        /metadata receives the first nonempty chunk's processing metadata,
        excluding AGC settings, which vary by frame and are stored separately.
        /data/cone_weights stores per-frame cone weights (N, 3).
        /data/scene_illuminance stores aligned scene illuminance (N,) in lux.
        /metadata/fps stores the world camera rate from config.pkl.
        /data/num_dark_and_saturated_pixels stores raw counts (N, 2): pixels
        <= the calibrated dark signal, then pixels >= 255, before transformation.

    Raises:
        AssertionError: World metadata has neither 4 legacy nor 6 modern columns.
        KeyError: Required world FPS configuration fields are missing.
        ValueError: FPS is invalid, an existing file is incomplete, frame counts disagree, or
            no nonempty chunks are selected.
        OSError: Input arrays cannot be read or the HDF5 file cannot be written.

    Notes:
        The root completed attribute becomes True only after all selected chunks
        are written. Exceptions close the file and pool; partial files remain
        marked incomplete. Each chunk is flushed before processing the next.
        Script callers using worker processes need an if __name__ == "__main__"
        guard. Worker benchmarks are recorded in tmp/benchmark_imputation_workers.json
        and tmp/benchmark_demosaic_large_buffers.json.
    """

    # Get the FPS the world captured at. We get this from the metadata file saved with the recording
    world_fps: float = float(recording_metadata["sensors"]["W"]["sensor_mode"]["fps"])
    if not np.isfinite(world_fps) or world_fps <= 0:
        raise ValueError("World camera FPS in config.pkl must be finite and positive.")

    # Load the floor value for pixel values that we calculated (the DarkSignal)
    dark_threshold: float = float(scipy.io.loadmat(
        Path(__file__).resolve().parents[3] / "derived" / "darkSignal.mat"
    )["darkSignal"].item())
    assert np.isfinite(dark_threshold), f"Dark threshold must be finite; got {dark_threshold}."

    def _write_hdf5_metadata(group: h5py.Group, values: dict[str, object]) -> None:
        """Recursively store this run's calibration metadata beneath an HDF5 group.

        Args:
            group: Open destination group; entry names must not already exist.
            values: Dictionary of scalars, arrays, nested dictionaries, or None.
                MATLAB structured/object arrays from scipy.loadmat are supported.

        Returns:
            None. Writes directly into group. Dictionaries and MATLAB structs
            become subgroups; object arrays become numbered entries with a shape
            attribute. None becomes a group with is_none=True. Strings are UTF-8
            datasets; ordinary numeric values keep their NumPy shape and dtype.
            These one-time metadata datasets are not compressed.

        Raises:
            TypeError: A value cannot be represented by the supported HDF5 types.
            ValueError: A name conflicts with an existing group or dataset.
            OSError: HDF5 cannot write to the output file.
        """

        # Examples below assume group is the file's /metadata group.
        # Handle one named dictionary entry at a time.
        for name, value in values.items():
            # Nested dictionaries need a container rather than an array dataset.
            # {"calibration": {"gain": 2.0}} -> /metadata/calibration/gain = 2.0
            if isinstance(value, dict):
                # Create that container and repeat these same steps for its entries.
                _write_hdf5_metadata(group.create_group(name), value)

            # HDF5 has no direct equivalent of Python's None value.
            # {"optional": None} -> empty /metadata/optional group, is_none=True attribute.
            elif value is None:
                # Mark an empty group so None is distinguishable from missing metadata.
                group.create_group(name).attrs["is_none"] = True
            else:
                # Normalize scalars/lists/arrays so we can inspect their NumPy dtype.
                array: np.ndarray = np.asarray(value)

                # Named fields identify a structured array, often a MATLAB struct.
                # {"fit": np.array([(2.0,)], dtype=[("gain", "f8")])}
                # -> /metadata/fit/gain is the one-element dataset [2.0].
                if array.dtype.names:
                    # Turn its fields into dictionary entries and write them recursively.
                    _write_hdf5_metadata(group.create_group(name), {
                        field: array[field] for field in array.dtype.names
                    })

                # Object arrays may contain cells with different shapes and types.
                # A two-cell object array named "cells", containing [1, 2] and [3],
                # -> /metadata/cells/0 = [1, 2], /metadata/cells/1 = [3];
                # the /metadata/cells group has shape=(2,) as an attribute.
                elif array.dtype.kind == "O":
                    # Give this collection its own subgroup.
                    child: h5py.Group = group.create_group(name)

                    # Preserve the original dimensions before flattening its elements.
                    child.attrs["shape"] = array.shape

                    # Number the cells in NumPy row-major order; recurse on each value.
                    _write_hdf5_metadata(child, {
                        str(index): item for index, item in enumerate(array.flat)
                    })

                # NumPy Unicode strings need an explicit HDF5 string encoding.
                # {"units": "W/m²/sr"} -> /metadata/units is a scalar UTF-8 dataset.
                elif array.dtype.kind == "U":
                    # Use UTF-8 to preserve characters such as the squared sign in units.
                    group.create_dataset(name, data=array.astype(object), dtype=h5py.string_dtype("utf-8"))
                else:
                    # Numeric scalars/arrays and byte strings can be written directly.
                    # {"bit_depth": 8} -> /metadata/bit_depth has scalar shape ().
                    # {"weights": [1.0, 2.0]} -> /metadata/weights has shape (2,).
                    # {"map": np.ones((2, 3))} -> /metadata/map has shape (2, 3).
                    # {"tag": b"RGB"} -> /metadata/tag is a scalar byte-string dataset.
                    group.create_dataset(name, data=array)

    # Get the start and end chunk we wish to process
    # while replacing the None sentiel if it exists
    start_chunk, end_chunk = chunk_range
    start_chunk = 0 if start_chunk is None else start_chunk
    end_chunk = len(files) if end_chunk is None else min(end_chunk, len(files))
    desired_chunk_range: object = range(start_chunk, end_chunk)

    # Decide whether to skip before creating the file; never check inside the loop.
    if os.path.exists(output_path) and not overwrite_existing:

        # We check if the previous existing output was fully completed here
        with h5py.File(output_path, "r") as existing:
            if not existing.attrs.get("completed", False):
                raise ValueError(f"Incomplete HDF5 output; retry with overwrite_existing=True: {output_path}")

            # Older opponent maps used per-frame neutral balancing. They must be
            # regenerated before the new analysis can reuse this completed file.
            if "metadata/balance_neutral" not in existing or bool(existing["metadata/balance_neutral"][()]):
                raise ValueError(
                    "World output does not record balance_neutral=False; regenerate "
                    f"with overwrite_existing=True: {output_path}"
                )
        return

    # Normalize gaze coordinates before opening or replacing the output file.
    # None lets the opponent-map helper use the image center. A single X/Y pair
    # is reused for every frame; otherwise, each input world frame needs its own pair.
    gaze: np.ndarray | None = (
        None if gaze_angles_xy is None else np.asarray(gaze_angles_xy, dtype=float)
    )

    if gaze is not None and gaze.shape != (2,):
        # Gaze rows refer to the entire input recording, including chunks outside
        # chunk_range. Memory mapping reads each buffer's length without loading it all.
        total_frames: int = sum(len(np.load(path, mmap_mode="r")) for _, path in files)

        if gaze.shape != (total_frames, 2):
            raise ValueError(
                "Gaze must have one X/Y pixel pair per input world frame; "
                f"expected {(total_frames, 2)}, got {gaze.shape}."
            )

    # Make the path to the output file if it does not already exist
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    # Add a progress bar to the chunk iterator when verbose output is enabled.
    file_iterator: Iterable[tuple[int, tuple[str, str]]] = enumerate(files)
    if verbose:
        file_iterator = tqdm(
            file_iterator,
            total=len(files),
            desc="Processing world chunks",
        )

    # Tune demosaicing separately from imputation and reuse its pool across chunks.
    # Worker counts are independent; context exit always releases the pool.
    processing_metadata: dict | None = None

    # "x" creates a new file and raises FileExistsError if it already exists;
    # "w" creates a file or erases/replaces an existing one when overwrite is requested.
    with h5py.File(output_path, "w" if overwrite_existing else "x") as output, \
        world_util.demosaic_worker_pool(n_workers=demosaic_n_workers) as demosaic_pool:
        # We mark the output as not completed to start with as a flag there
        # is still procesing to do with it
        output.attrs["completed"] = False

        # Create groups which we will populate with datasets later
        data_group: h5py.Group = output.create_group("data")
        agc_group: h5py.Group = output.create_group("agcSettings")
        metadata_group: h5py.Group = output.create_group("metadata")

        # Keep track of how many frames we have written
        frames_written: int = 0

        # Track this chunk's first gaze row in the full recording, not the output.
        gaze_offset: int = 0

        # Iterate over the paired metadata and frame-buffer files.
        for file_num, (metadata_path, data_path) in file_iterator:
            frame_count: int = len(np.load(data_path, mmap_mode="r"))
            chunk_gaze: np.ndarray | None

            # Preserve the image-center default or the shared X/Y pair unchanged.
            # For per-frame gaze, select the rows corresponding to this chunk.
            if gaze is None or gaze.shape == (2,):
                chunk_gaze = gaze
            else:
                chunk_gaze = gaze[gaze_offset:gaze_offset + frame_count]

            # Advance even when this chunk is skipped. Otherwise later chunks would
            # receive gaze rows from earlier frames and be centered incorrectly.
            gaze_offset += frame_count

            # Skip files not in the desired range
            if(file_num not in desired_chunk_range):
                continue

            # Load the metadata and its associated frame buffer.
            metadata_buffer: np.ndarray = np.load(metadata_path)
            frame_buffer: np.ndarray = np.load(data_path)

            # The final chunk of a recording may be empty.
            if len(frame_buffer) == 0:
                continue

            if len(metadata_buffer) != len(frame_buffer):
                raise ValueError(f"Frame/metadata length mismatch: {data_path}")

            # The metadata is either legacy shaped or modern shaped
            # Legacy shape is timestamp, Again, DGain, Exposure
            # Modern shape is timestamp, "cameraAgain", "AGCDgain", "cameraExposure", "AGCAgain", "AGCExposure"
            Again_idx: int
            Dgain_idx: int
            exposure_idx: int

            # First check to see if the metadata is properly shaped
            assert metadata_buffer.shape[1] in (4, 6), f"Metadata buffer must have cols: (timestamp, Again, DGain, Exposure) or timestamp, cameraAgain, AGCDgain, cameraExposure, AGCAgain, AGCExposure"

            # Check to see if the metadata buffer is of the old format
            if(metadata_buffer.shape[1] == 4):
                Again_idx = 1
                Dgain_idx = 2
                exposure_idx = 3

            # or the NEW format
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

            # Count raw Bayer pixels before the pipeline can modify their values.
            assert frame_buffer.ndim == 3, (
                f"Raw Bayer frames must have shape (N, height, width); got {frame_buffer.shape}."
            )

            # Let's now count the the number of dark and saturated pixels
            dark_and_saturated_counts: np.ndarray = np.zeros((len(frame_buffer), 2), dtype=np.int64)
            for frame_index, frame in enumerate(frame_buffer):
                # One coordinate entry per matching raw pixel: dark first, saturated second.
                dark_and_saturated_counts[frame_index, 0] = len(np.where(frame <= dark_threshold)[0])
                dark_and_saturated_counts[frame_index, 1] = len(np.where(frame >= 255)[0])

            # Pass the aligned gaze rows with this frame buffer so the opponent maps
            # are centered on the corresponding gaze position for each world frame.
            data_and_processing_metadata: dict[str, object] = world_util.world_transformation_pipeline(
                frame_buffer, agc_settings, chunk_gaze, age, n_workers=n_workers, demosaic_pool=demosaic_pool,
            )
            transformed_data: np.ndarray = data_and_processing_metadata["data"]

            # Save calibration/settings from the first processed chunk only because it's the
            # same for all processed chunks
            # Exclude AGC here because its per-frame values are appended separately.
            if processing_metadata is None:
                processing_metadata = {
                    key: value for key, value in data_and_processing_metadata["metadata"].items()
                    if key != "agc_settings"
                }

                # Add the FPS to the metadata that we will store
                processing_metadata["fps"] = world_fps

                # Write all the metadata to the HDF5
                _write_hdf5_metadata(metadata_group, processing_metadata)

            # If the shapes are not matching anymore, that means something has gone wayyy out of wack
            if len(transformed_data) != len(frame_buffer):
                raise ValueError("Transformation must return one output per input frame")

            # If we have not yet created the frames dataset, create it now
            # HDF5 compresses each storage chunk independently; shuffle improves
            # compression of numeric arrays without changing their values.
            if "frames" not in data_group:
                data_group.create_dataset(
                    "frames", shape=(0, *transformed_data.shape[1:]),
                    maxshape=(None, *transformed_data.shape[1:]),
                    chunks=True, dtype=np.float64,
                    compression="gzip", compression_opts=4, shuffle=True,
                )

            # Append this data buffer
            data_dataset: h5py.Dataset = data_group["frames"]
            end_frame: int = frames_written + len(transformed_data)

            # Reset the dataset to accomodate the new frames we are going to put in
            data_dataset.resize(end_frame, axis=0)

            # Preserve fractional, signed radiance and NaNs outside camera coverage.
            data_dataset[frames_written:end_frame] = transformed_data

            # Scene balancing varies per frame, so save these weights alongside maps.
            if "cone_weights" not in data_group:
                data_group.create_dataset("cone_weights", shape=(0, 3), maxshape=(None, 3),
                                          chunks=True, dtype=np.float64)
            weights_dataset = data_group["cone_weights"]
            weights_dataset.resize(end_frame, axis=0)
            weights_dataset[frames_written:end_frame] = data_and_processing_metadata["cone_weights"]

            # Save full-frame camera illuminance alongside the same frame timestamps.
            if "scene_illuminance" not in data_group:
                data_group.create_dataset(
                    "scene_illuminance", shape=(0,), maxshape=(None,),
                    chunks=True, dtype=np.float64,
                    compression="gzip", compression_opts=4, shuffle=True,
                ).attrs["units"] = "lux"
            illuminance_dataset: h5py.Dataset = data_group["scene_illuminance"]
            illuminance_dataset.resize(end_frame, axis=0)
            illuminance_dataset[frames_written:end_frame] = data_and_processing_metadata["scene_illuminance"]

            # MS fitting still requires camera RGB radiance, not opponent channels.
            mean_rgb: np.ndarray = data_and_processing_metadata["average_rgb"]
            if "average_frames" not in data_group:
                data_group.create_dataset(
                    "average_frames", shape=(0, 3), maxshape=(None, 3),
                    chunks=True, dtype=mean_rgb.dtype,
                    compression="gzip", compression_opts=4, shuffle=True,
                )
            averages_dataset: h5py.Dataset = data_group["average_frames"]
            averages_dataset.resize(end_frame, axis=0)
            averages_dataset[frames_written:end_frame] = mean_rgb

            # Append dark then saturated counts in the same frame order as the maps.
            if "num_dark_and_saturated_pixels" not in data_group:
                counts_dataset: h5py.Dataset = data_group.create_dataset(
                    "num_dark_and_saturated_pixels", shape=(0, 2), maxshape=(None, 2),
                    chunks=True, dtype=np.int64,
                    compression="gzip", compression_opts=4, shuffle=True,
                )
                counts_dataset.attrs["columns"] = ["dark", "saturated"]
                counts_dataset.attrs["dark_threshold"] = dark_threshold
                counts_dataset.attrs["saturation_threshold"] = 255
            counts_dataset = data_group["num_dark_and_saturated_pixels"]
            counts_dataset.resize(end_frame, axis=0)
            counts_dataset[frames_written:end_frame] = dark_and_saturated_counts

            # Append this AGC buffer
            for name, values in {**agc_settings, "timestamps": metadata_buffer[:, 0]}.items():
                # If this dataset does not already exist, create it
                if name not in agc_group:
                    agc_group.create_dataset(
                        name, shape=(0,), maxshape=(None,), chunks=True,
                        dtype=values.dtype,
                        compression="gzip", compression_opts=4, shuffle=True,
                    )

                # Get this dataset, resize, and write to it
                dataset: h5py.Dataset = agc_group[name]
                dataset.resize(end_frame, axis=0)
                dataset[frames_written:end_frame] = values

            # Update the number of frames written
            frames_written = end_frame

            # Flush HDF5 buffers to the file without closing it; the OS may still cache writes.
            output.flush()

            # Drop every reference to this large result before processing the next buffer.
            del transformed_data, data_and_processing_metadata, frame_buffer, metadata_buffer, agc_settings

        # Ensure we at least wrote some frames in the course of the video
        if "frames" not in data_group:
            raise ValueError("No nonempty world chunks selected; no transformed data written")

        # Mark this dataset as completed
        # so that it's clear it was completed without error
        output.attrs["completed"] = True

        return



def _process_raw_ms_helper(
    recording_metadata: dict,
    ms_files: list[tuple[str | Path, str | Path]],
    world_files: list[tuple[str | Path, str | Path]],
    output_path: str | Path,
    chunk_range: tuple[int | None, int | None] = (0, None),
    overwrite_existing: bool = False,
    verbose: bool = False,
    compression_level: int = 4,
    base_alpha: float = ms_util.MS_JOINT_BASE_ALPHA,
    reference_intensity: float = ms_util.MS_JOINT_REFERENCE_INTENSITY,
    beta: float = ms_util.MS_JOINT_BETA,
) -> None:
    """Fit joint MS/camera spectra and append chunks to one native HDF5 file.

    Args:
        recording_metadata: Configuration from config.pkl, accepted for parity with world
            processing and future metadata support; currently unused.
        ms_files: Pairs of string or Path filenames. Ordered (timestamp .npy path, reading-buffer .npy path) pairs.
            Each timestamp array flattens to (N,) in seconds on the world
            recording clock. Each reading buffer has shape (N, channels), with
            the first seven columns containing AS7341 F1-F7 measurements.
        output_path: String or Path. MS destination directory, typically GKA/<attempt>/M.
            Processed world data are read from its sibling W directory.
        chunk_range: Zero-based (start, stop) chunk indices, stop exclusive;
            None means beginning/end. Defaults to all chunks.
        overwrite_existing: Replace the entire MS HDF5 file when True; otherwise skip a
            completed file and raise for an incomplete file.
        verbose: Show progress over input chunk pairs.
        compression_level: Gzip compression level (0-9), default 4.
        base_alpha: Joint-fit smoothness strength; defaults to the MATLAB model setting.
        reference_intensity: Reference intensity used to scale smoothness.
        beta: Joint-fit shrinkage strength; defaults to the MATLAB model setting.
        world_files: Pairs of string or Path filenames. Raw world metadata/data pairs required for camera alignment.
            Requires the sibling W/transformed.hdf5 with /data/average_frames.

    Returns:
        None. Writes transformed.hdf5 with these /data fields (NumPy shapes):
            counts/readings: (N, channels) complete raw MS reading buffers.
            illuminance/lux: (N,) photopic lux derived from the joint spectrum.
            illuminance/melanopic_edi: (N,) melanopic equivalent daylight illuminance in lux.
            radiance: (N, L) spectral radiance in W/m²/sr/nm, L wavelengths.
            S: (3,) sampling descriptor [start_nm, step_nm, wavelength_count].
            fVal: (N,) Euclidean norms of the ten measurement residuals; these
                are fit errors, not percentages or the regularized objective.
            fitErrors: (N, 10) signed observed-minus-predicted residuals for
                seven MS channels followed by camera R, G, B.
        /timestamps stores the MS times in seconds. /metadata stores window-mean
        RGB (N, 3), world_frame_index_bounds (N, 2) indexing /data/average_frames with zero-based
        start-inclusive/stop-exclusive indices, world frame counts
        (N,), MS reading time bounds in seconds (N, 2), and fit settings.
        Covered windows require finite averaged RGB values. Empty buffers are
        skipped; uncovered windows retain NaN spectra and fit diagnostics.

    Raises:
        AssertionError: A covered MS window has nonfinite averaged RGB.
        ValueError: World input list is empty, timestamp/reading counts differ,
            or camera-window validation fails.
        RuntimeError: The joint spectrum estimator fails to converge.
        OSError: A required input cannot be read or an output cannot be written.

    Notes:
        Window boundaries use neighboring readings even across chunk boundaries.
        Neighbor metadata is loaded on demand, including outside a selected chunk
        range, so selecting a chunk does not change its alignment windows.
        For MS times [5, 15, 25], midpoint boundaries are 10 and 20: world
        times 1..30 form windows 1..10, 11..20, and 21..30. The first and last
        windows extend to available recording edges; no time offset is added.
        Spatial RGB means are averaged with equal weight per camera frame,
        including frames spanning world chunk boundaries. World timestamps span
        the recording; MS timestamps and readings are loaded chunkwise.
        Windows without world frames retain their rows with NaN spectra and fit
        diagnostics. Invalid measurements in covered windows still cause errors.
        Per-reading arrays use extendible chunks, the requested gzip level, and shuffle.
        S and constant fit metadata are written once. completed=True is set only
        after all selected chunks have been appended; failed ms_files require overwrite.
    """

    # 1. Resolve the requested half-open chunk range: start included, stop excluded.
    start_chunk: int | None
    end_chunk: int | None
    start_chunk, end_chunk = chunk_range
    start_chunk = 0 if start_chunk is None else start_chunk
    end_chunk = len(ms_files) if end_chunk is None else min(end_chunk, len(ms_files))
    desired_chunk_range: range = range(start_chunk, end_chunk)


    # See if the output path exists, is completed, and if we should overwrite or not
    output_filepath: str = os.path.join(output_path, "transformed.hdf5")
    if os.path.exists(output_filepath) and not overwrite_existing:
        with h5py.File(output_filepath, "r") as existing:
            if not existing.attrs.get("completed", False):
                raise ValueError(f"Incomplete MS HDF5 output; retry with overwrite_existing=True: {output_filepath}")
            for required_path in ("data/counts/readings", "data/illuminance/lux", "data/illuminance/melanopic_edi"):
                if required_path not in existing:
                    raise ValueError(
                        f"MS output lacks {required_path}; regenerate with "
                        f"overwrite_existing=True: {output_filepath}"
                    )
        return

    # Ensure the minispectrometer output directory exists.
    os.makedirs(output_path, exist_ok=True)

    # Add a progress bar to the chunk iterator when verbose output is enabled.
    file_iterator: Iterable[tuple[int, tuple[str, str]]] = enumerate(ms_files)
    if verbose:
        file_iterator = tqdm(
            file_iterator,
            total=len(ms_files),
            desc="Processing MS chunks",
        )

    # World output lives beside the MS output: .../M -> .../W.
    world_output_path: str = os.path.join(
        os.path.dirname(os.path.abspath(output_path)), "W"
    )
    if not world_files:
        raise ValueError("Joint spectral reconstruction requires world files")

    world_filepath: str = os.path.join(world_output_path, "transformed.hdf5")

    # Create the file before processing buffers; failures leave completed=False.
    # "x" creates a new file and raises FileExistsError if it already exists;
    # "w" creates a file or erases/replaces an existing one when overwrite is requested.
    with (
        h5py.File(world_filepath, "r") as world,
        h5py.File(output_filepath, "w" if overwrite_existing else "x") as output,
    ):
        # Initialize the output for this HDF5 file, we will have a data group
        # and a metadata group
        output.attrs["completed"] = False
        data_group: h5py.Group = output.create_group("data")
        data_group.create_group("counts")
        illuminance_group: h5py.Group = data_group.create_group("illuminance")
        illuminance_group.attrs["units"] = "lux"
        metadata_group: h5py.Group = output.create_group("metadata")

        # Track the next output row as successive MS chunks are appended.
        output_row_start: int = 0

        # The world processing must be done at this point to be considered
        # ready to process the MS
        if not world.attrs.get("completed", False):
            raise ValueError(f"Incomplete world output: {world_filepath}")

        # Load in the world timestsamps. We will need this to calculate what
        # world timestamps overlap with each MS reading
        world_timestamps: np.ndarray = world["agcSettings/timestamps"][:] * WORLD_NANOSECONDS_TO_SECONDS

        # Ensure nothing weird is happening with the world tiemstamps
        if not np.all(np.isfinite(world_timestamps)) or np.any(np.diff(world_timestamps) <= 0):
            raise ValueError("World timestamps must be finite and strictly increasing")

        # We will only load in the precomputed averages so processing is faster
        average_frames: h5py.Dataset = world["data/average_frames"]

        # Ensure the sahpe is NFrames, color1, color2, color3
        if average_frames.shape != (len(world_timestamps), 3):
            raise ValueError("World averages must have one RGB row per timestamp")

        # Carry the last processed timestamp across readings and chunk boundaries.
        # None means there is no earlier reading; zero is a valid timestamp.
        previous_timestamp: float | None = None

        # When starting partway through a recording, use the last timestamp
        # from the immediately preceding chunk.
        if start_chunk != 0:
            previous_chunk_timestamps: np.ndarray = np.load(
                ms_files[start_chunk - 1][0]
            ).reshape(-1)
            if len(previous_chunk_timestamps) == 0:
                raise ValueError("The chunk preceding start_chunk has no MS timestamps")

            if not np.all(np.isfinite(previous_chunk_timestamps)) or np.any(np.diff(previous_chunk_timestamps) <= 0):
                raise ValueError(f"MS timestamps must be finite and increasing: {metadata_path}")

            previous_timestamp = float(previous_chunk_timestamps[-1])
            del previous_chunk_timestamps

        # 3. Load, fit, and save one MS chunk at a time to keep memory use bounded.
        for file_num, (metadata_path, data_path) in file_iterator:

            # Skip ms_files not in the desired range
            if file_num not in desired_chunk_range:
                continue

            # Load the timestamps and their associated minispectrometer readings.
            timestamps: np.ndarray = np.load(metadata_path).reshape(-1)
            reading_buffer: np.ndarray = np.load(data_path)

            # Every reading must have a timestamp, including in empty chunks.
            if len(timestamps) != len(reading_buffer):
                raise ValueError(f"MS timestamp/reading count mismatch: {data_path}")

            # Skip an empty chunk
            if len(timestamps) == 0:
                continue

            # Ensure nothing weird is going on with the timestamps of the MS
            if not np.all(np.isfinite(timestamps)) or np.any(np.diff(timestamps) <= 0):
                raise ValueError(f"MS timestamps must be finite and increasing: {metadata_path}")

            # Preserve every raw channel for fitting and later inspection.
            assert reading_buffer.ndim == 2, (
                f"MS readings must be two-dimensional; got {reading_buffer.shape}: {data_path}"
            )
            assert reading_buffer.shape[1] >= ms_util.MS_AS_CHANNELS, (
                f"MS readings require at least {ms_util.MS_AS_CHANNELS} AS channels; "
                f"got {reading_buffer.shape[1]}: {data_path}"
            )

            # Joint-spectrum lux and melanopic EDI are populated after fitting below.
            # Missing camera coverage leaves these estimates NaN, like the spectrum.
            illuminance_lux: np.ndarray = np.full(len(timestamps), np.nan)
            melanopic_edi: np.ndarray = np.full(len(timestamps), np.nan)

            # 3a. For each reading in this buffer, find its first world frame and the index just
            # after its last frame. world_frame_index_bounds is (N, 2), using global world row indices.
            world_frame_index_bounds: np.ndarray = np.empty((len(timestamps), 2), dtype=np.int64)
            ms_reading_bounds: np.ndarray = np.empty((len(timestamps), 2), dtype=float)

            # Iterate over the reading num and current time for this reading
            for reading_index, current_time in enumerate(timestamps):

                next_time: float | None = None

                # At the last reading, we have to read into the next chunk to get the timestamp
                if reading_index == len(timestamps) - 1:
                    # Read only the immediately following chunk, if one exists.
                    if file_num + 1 < len(ms_files):
                        next_chunk_timestamps: np.ndarray = np.load(
                            ms_files[file_num + 1][0]
                        ).reshape(-1)

                        # Make sure nothing weird has gone on with them
                        if len(next_chunk_timestamps) == 0:
                            raise ValueError("The next MS chunk has no timestamps")

                        if not np.all(np.isfinite(next_chunk_timestamps)) or np.any(np.diff(next_chunk_timestamps) <= 0):
                            raise ValueError(f"MS timestamps must be finite and increasing: {metadata_path}")

                        next_time = float(next_chunk_timestamps[0])

                        del next_chunk_timestamps

                # Otherwise, inside the chunk we just look at the next timestamp
                else:
                    # Inside the chunk, the next timestamp is already in memory.
                    next_time = float(timestamps[reading_index + 1])

                # Define the camera averaging window associated with this MS reading.
                # The camera samples more frequently than the MS sensor, so several
                # world frames contribute to the RGB average paired with one MS row.
                # We assign frames using the halfway points between neighboring MS
                # timestamps: these separate frames nearest to this reading from
                # frames nearest to the previous or next reading.
                #
                # Example: for MS timestamps [5, 15, 25] seconds, the reading at
                # 15 seconds gets a window from 10 to 20 seconds. Our searches below
                # use (10, 20]: exclude a world frame at 10, include one at 20.
                # This assigns a frame exactly at a shared midpoint to the earlier
                # MS reading, avoiding double-counting it in adjacent windows.
                #
                # Each reading gets its own window; there is no fixed window width.
                # For irregular MS times [5, 15, 35], the reading at 15 instead
                # uses (10, 25]. Neighbors can come from adjacent chunk files;
                # splitting a recording into files must not change these windows.
                #
                # At the recording edges, there may be no previous or next reading.
                # We use -inf or +inf for that missing boundary so the search includes
                # the available world frames at that edge; it creates no extra frames.
                # This is a nearest-timestamp averaging convention, not a measurement
                # of the MS sensor's actual integration/exposure interval.
                #
                # First calculate these time bounds in seconds, save them in
                # ms_reading_bounds, then search the sorted world timestamps to get
                # world_frame_index_bounds. Those indices select the saved RGB rows
                # that will be averaged before fitting this reading's spectrum.

                # The lower bound is halfway between the previous and current MS times.
                # Example: previous=5 s and current=15 s gives a lower bound of 10 s.
                left_time: float
                if previous_timestamp is None:
                    # This is the first MS reading: impose no lower time limit.
                    # The world timestamp search will therefore start at frame index 0.
                    left_time = -np.inf
                else:
                    # A previous reading must have a finite, earlier timestamp.
                    if not np.isfinite(previous_timestamp) or previous_timestamp >= current_time:
                        raise ValueError("The previous MS timestamp must be finite and earlier than the current timestamp")

                    # Move halfway forward from the previous reading toward this one.
                    previous_gap: float = float(current_time - previous_timestamp)
                    left_time = previous_timestamp + previous_gap / 2

                # The upper bound is halfway between the current and next MS times.
                # Example: current=15 s and next=25 s gives an upper bound of 20 s.
                right_time: float
                if next_time is None:
                    # This is the last MS reading: impose no upper time limit.
                    # The world timestamp search will include all remaining frames.
                    right_time = np.inf
                else:
                    # A next reading must have a finite, later timestamp.
                    if not np.isfinite(next_time) or next_time <= current_time:
                        raise ValueError("The next MS timestamp must be finite and later than the current timestamp")

                    # Move halfway forward from this reading toward the next one.
                    next_gap: float = float(next_time - current_time)
                    right_time = current_time + next_gap / 2

                # Save this MS reading's averaging window in seconds, e.g. [10, 20].
                # These are midpoint boundaries, not the sensor's exposure start/end.
                ms_reading_bounds[reading_index] = (left_time, right_time)

                # side="right" preserves (left, right] time windows: a frame
                # exactly at a midpoint belongs to the earlier MS reading.
                world_frame_index_bounds[reading_index, 0] = np.searchsorted(
                    world_timestamps, left_time, side="right"
                )
                world_frame_index_bounds[reading_index, 1] = np.searchsorted(
                    world_timestamps, right_time, side="right"
                )

                # The next iteration reuses this value, even if it starts a new chunk.
                previous_timestamp = float(current_time)

            # 3b. Average each range of saved RGB rows. For world_frame_index_bounds [0, 10],
            # average_frames[0:10] includes rows 0 through 9 (ten camera frames).
            frame_counts: np.ndarray = world_frame_index_bounds[:, 1] - world_frame_index_bounds[:, 0]
            has_world_frames: np.ndarray = frame_counts > 0

            # Preserve missing windows as NaN without averaging an empty slice.
            rgb: np.ndarray = np.full((len(timestamps), 3), np.nan, dtype=float)
            for reading_index in np.flatnonzero(has_world_frames):
                first, stop = world_frame_index_bounds[reading_index]
                rgb[reading_index] = average_frames[first:stop].mean(axis=0)

            assert np.all(np.isfinite(rgb[has_world_frames])), "Nonfinite camera RGB averages in MS windows"

            # 3c. Fit only covered readings, retaining every original output row.
            # Initialize explicitly so even an entirely uncovered chunk can be saved.
            sampling: np.ndarray = ms_util.MS_SPECTRAL_SAMPLING.copy()
            radiance: np.ndarray = np.full(
                (len(timestamps), ms_util.MS_SPECTRAL_WAVELENGTHS.size), np.nan,
            )
            fit_value: np.ndarray = np.full(len(timestamps), np.nan)
            fit_errors: np.ndarray = np.full(
                (len(timestamps), ms_util.MS_JOINT_CHANNEL_COUNT + 3), np.nan,
            )
            if np.any(has_world_frames):
                (
                    fitted_radiance, sampling, fitted_lux, fitted_melanopic_edi,
                    fitted_value, fitted_errors,
                ) = ms_util.estimate_joint_radiance_spectrum(
                    reading_buffer[has_world_frames], rgb[has_world_frames],
                    base_alpha=base_alpha, reference_intensity=reference_intensity, beta=beta,
                )
                radiance[has_world_frames] = fitted_radiance
                illuminance_lux[has_world_frames] = fitted_lux
                melanopic_edi[has_world_frames] = fitted_melanopic_edi
                fit_value[has_world_frames] = fitted_value
                fit_errors[has_world_frames] = fitted_errors

            # 5. Keep fixed run settings separate from arrays that grow with each chunk.
            # These settings describe the actual fit above; calibration coefficients
            # come directly from ms_util, rather than being copied by hand here.
            processing_settings: dict[str, object] = {
                "world_frame_index_bounds_convention": "zero-based world row indices, start inclusive, stop exclusive",
                "estimator": "estimateJointRadianceSpectrum",
                "illuminance_method": "joint spectrum projected through photoreceptor sensitivities",
                "photoreceptor_calibration": "sensorToPhotoreceptorMatrix.mat",
                "window_convention": "left-open, right-closed; outer time bounds extend to recording edges",
                "window_time_units": "seconds",
                "world_timestamp_scale": WORLD_NANOSECONDS_TO_SECONDS,
                "radiance_units": "W/m^2/sr/nm",
                "base_alpha": base_alpha,
                "reference_intensity": reference_intensity,
                "beta": beta,
                "solver_tolerance": ms_util.MS_JOINT_SOLVER_TOLERANCE,
                "solver_max_iterations": ms_util.MS_JOINT_SOLVER_MAX_ITERATIONS,
                "minispect_radiance_weights": ms_util.MS_RADIANCE_WEIGHTS[:ms_util.MS_JOINT_CHANNEL_COUNT],
            }

            # Append each per-reading array along its first axis, preserving float precision.
            chunk_arrays_by_path: dict[str, np.ndarray] = {
                "timestamps": timestamps,
                "data/counts/readings": reading_buffer,
                "data/illuminance/lux": illuminance_lux,
                "data/illuminance/melanopic_edi": melanopic_edi,
                "data/radiance": radiance,
                "data/fVal": fit_value,
                "data/fitErrors": fit_errors,
                "metadata/world_frame_index_bounds": world_frame_index_bounds,
                "metadata/mean_camera_rgb": rgb,
                "metadata/world_frame_counts": frame_counts,
                "metadata/ms_reading_bounds": ms_reading_bounds,
            }

            # Reserve the same output rows in every dataset. If 100 readings were
            # already saved and this chunk has 20, write rows [100:120].
            output_row_stop: int = output_row_start + len(timestamps)

            for dataset_path, chunk_values in chunk_arrays_by_path.items():

                # The first chunk creates an empty dataset with this array's trailing
                # dimensions. None allows only the first (reading) axis to grow.
                if dataset_path not in output:
                    output.create_dataset(
                        dataset_path, shape=(0, *chunk_values.shape[1:]),
                        maxshape=(None, *chunk_values.shape[1:]), dtype=chunk_values.dtype,
                        chunks=True, compression="gzip", compression_opts=compression_level, shuffle=True,
                    )

                # Extend the existing dataset and write this chunk into its new rows.
                dataset: h5py.Dataset = output[dataset_path]
                dataset.resize(output_row_stop, axis=0)
                dataset[output_row_start:output_row_stop] = chunk_values

            # The wavelength support and fixed estimator settings do not vary by chunk.
            # The first chunk writes S and fixed settings once. Later chunks only
            # check that their wavelength support still matches the saved spectra.
            if "S" not in data_group:
                data_group.create_dataset(
                    "S", data=sampling, chunks=True,
                    compression="gzip", compression_opts=compression_level, shuffle=True,
                )

                # Strings need UTF-8 storage; numeric settings retain their dtype.
                for setting_name, setting_value in processing_settings.items():
                    if isinstance(setting_value, str):
                        metadata_group.create_dataset(
                            setting_name, data=setting_value, dtype=h5py.string_dtype("utf-8")
                        )
                    else:
                        setting_array: np.ndarray = np.asarray(setting_value)

                        # HDF5 scalar datasets do not support chunk compression.
                        dataset_options: dict[str, object] = {}
                        if setting_array.ndim > 0:
                            dataset_options = dict(
                                chunks=True, compression="gzip",
                                compression_opts=compression_level, shuffle=True,
                            )
                        metadata_group.create_dataset(setting_name, data=setting_array, **dataset_options)

            elif not np.array_equal(data_group["S"][:], sampling):
                raise ValueError("MS wavelength support changed between chunks")

            # Advance the cursor only after all arrays in this chunk are written.
            output_row_start = output_row_stop

            # Push HDF5 buffers to the file while keeping it open for the next chunk.
            output.flush()

        if output_row_start == 0:
            raise ValueError("No nonempty MS chunks selected")

        # Mark success only after every selected chunk has been saved.
        output.attrs["completed"] = True

    return



def process_raw_recording(
    path_to_raw: str | Path,
    output_path: str | Path,
    age: int,
    gaze_angles_xy: np.ndarray,
    overwrite_existing: bool = False,
    verbose: bool = False,
    chunk_ranges: dict[Literal["W", "M"], tuple[int | None]] = {sensor_name: (0, None) for sensor_name in "WM"},
    n_workers: int = world_util.WORLD_IMPUTATION_WORKERS,
    demosaic_n_workers: int = world_util.WORLD_DEMOSAIC_WORKERS,
) -> None:
    """Process world-camera data first, then joint minispectrometer spectra.

    Args:
        path_to_raw: String or Path. Existing nonempty directory of raw world/MS chunks and config.pkl.
            The configuration must be a dictionary with world FPS at
            ["sensors"]["W"]["sensor_mode"]["fps"]; it is passed to both processors.
        output_path: String or Path. Parent output directory. World data go to W/transformed.hdf5;
            MS results go to M/transformed.hdf5.
        overwrite_existing: Replace both recording HDF5 files when True. Otherwise skip
            completed files independently and raise for incomplete files.
        verbose: Display chunk progress bars.
        chunk_ranges: W and M keys mapped to zero-based (start, stop) chunk ranges,
            stop exclusive. None means beginning/end; defaults to all chunks.
        n_workers: World imputation workers, default WORLD_IMPUTATION_WORKERS (6).
        demosaic_n_workers: Independent world demosaicing workers, default
            WORLD_DEMOSAIC_WORKERS (16). Set both counts to 1 for serial execution.

        age: Required participant age forwarded to world reconstruction.
        gaze_angles_xy: Camera pixel X/Y pairs aligned to all world frames (one-based).

    Returns:
        None. Writes native HDF5 world and MS output files, including world
        FPS in W/transformed.hdf5 at /metadata/fps.
        Pupil files may be discovered but are not processed by this function.

    Raises:
        AssertionError: The input directory is missing/empty or world metadata
            has an unsupported column count.
        ValueError: World output is incomplete, selected world chunks are empty,
            or data validation/alignment fails in either processing helper.
        TypeError: config.pkl does not contain a dictionary.
        KeyError: Required configuration fields or chunk-range keys are missing.
        OSError: Input or output files cannot be accessed.

    Notes:
        Selected world chunks must cover the requested MS windows. MS uses saved
        per-frame RGB averages when available. Script callers using worker
        processes need an if __name__ == "__main__" guard.
    """

    # Ensure the raw recording is an existing, non-empty directory.
    assert (
        os.path.isdir(path_to_raw) and len(os.listdir(path_to_raw)) > 0
    ), f"{path_to_raw} must be an existing, non-empty directory"

    # Both sensor processors receive the same recording-level configuration.
    with open(os.path.join(path_to_raw, "config.pkl"), "rb") as config_file:
        recording_metadata: dict = dill.load(config_file)
    if not isinstance(recording_metadata, dict):
        raise TypeError("config.pkl must contain a recording metadata dictionary.")

    # Find the raw chunk files for each sensor.
    files_by_sensor: dict[str, list[tuple[str, str]]] = group_sensors_files(
        path_to_raw
    )

    if not files_by_sensor["W"]:
        raise ValueError(f"No world-camera chunks found in: {path_to_raw}")
    if not files_by_sensor["M"]:
        raise ValueError(f"No minispectrometer chunks found in: {path_to_raw}")

    # Complete and close the world file before MS opens it to read the RGB averages.
    world_output_path: str = os.path.join(output_path, "W")
    _process_raw_world_helper(
        recording_metadata,
        files_by_sensor["W"], os.path.join(world_output_path, "transformed.hdf5"),
        age=age, gaze_angles_xy=gaze_angles_xy, chunk_range=chunk_ranges["W"],
        overwrite_existing=overwrite_existing, verbose=verbose,
        n_workers=n_workers, demosaic_n_workers=demosaic_n_workers,
    )

    # This still runs when world processing skips an already completed file.
    _process_raw_ms_helper(
        recording_metadata,
        ms_files=files_by_sensor["M"], world_files=files_by_sensor["W"],
        output_path=os.path.join(output_path, "M"),
        chunk_range=chunk_ranges["M"], overwrite_existing=overwrite_existing,
        verbose=verbose,
    )


def _plot_raw_world_helper(
    path_to_raw: str | Path,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Plot raw world channels and camera settings in a 2-by-1 layout.

    Args:
        path_to_raw: String or Path to raw world chunks and config.pkl.
        axes: Two primary axes in panel order, or None to create a 2-by-1 figure.
        verbose: Show loading and plotting progress bars.


    Returns:
        Figure and primary axes array of shape (2,). Tile 0 shows raw Bayer
        R, G, B and whole-frame means. Tile 1 shows analog/digital gains on
        the left (0-11) and exposure on the right (0-10000 microseconds), with
        a combined legend. All world tiles use time in seconds and frame numbers
        above. The exposure twin axis is not included in the returned array.
        Does not show or save the figure.

    Raises:
        ValueError: The axes count, frame data, or metadata are incompatible.
        KeyError: Required recording configuration fields are missing.
        OSError: Raw chunks or config.pkl cannot be read.
    """

    # The world metadata reader expects a string rather than a Path object.
    path_to_raw = os.fspath(path_to_raw)

    # Validate supplied axes before setting up the figure.
    if axes is not None:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        assert axes.size == 2, f"Raw world plotting requires two axes; got {axes.size}." f" Recording: {path_to_raw}"

    # If we didn't receive axes to plot on, we will create our own figure now
    figure: Figure
    if axes is None:
        figure, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)

    # Otherwise, use the figure belonging to the supplied axes.
    else:
        figure = axes[0].figure

    # First, we get the world metadata from the chunks. This would include things
    # such as the timestamp and AGC settings. We do not fill missing frames here
    # with NaN, as it's unneeded
    metadata: pd.DataFrame = world_util.world_metadata_from_chunks(
        path_to_raw, fill_missing_frames=False, verbose=verbose,
    )

    # Then, we plot the raw RGB channels of the camera sensor
    world_util.plot_world_camera_channels(path_to_raw, ax=axes[0], verbose=verbose)


    # The settings plot handles both value scales and their combined legend.
    world_util.plot_world_camera_settings(
        metadata, ax=axes[1], time_origin=0, verbose=verbose,
    )
    return figure, axes


def _plot_raw_ms_helper(
    path_to_raw: str | Path,
    ax: Axes | None = None,
    verbose: bool = False,
) -> tuple[Figure, Axes]:
    """Plot raw AS spectral channel counts in a single tile.

    Args:
        path_to_raw: String or Path to paired raw MS timestamp and reading chunks.
        ax: Destination axis, or None to create a 1-by-1 figure.
        verbose: Show chunk loading and plotting progress bars.

    Returns:
        Figure and primary Axes. The single tile shows each AS channel's counts
        against its recorded time in seconds, with a channel legend. TS, LS,
        and temperature are omitted. Does not show or save the figure.

    Raises:
        ValueError: No nonempty chunks exist or readings and timestamps disagree.
        OSError: The directory or chunk files cannot be read.
    """

    # If no ax provided, make a figure to plot on
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)

    # Create arrays to hold the timestamp buffer and the readings buffer as
    # we read them in
    timestamp_buffer: list[np.ndarray] = []
    reading_buffer: list[np.ndarray] = []

    # Use the MS timestamps directly; they are already in seconds.
    for metadata_path, data_path in tqdm(group_sensors_files(path_to_raw)["M"], desc="Loading raw MS chunks", disable=not verbose):
        # Read in the 1D timestamp buffer
        timestamps: np.ndarray = np.load(metadata_path).reshape(-1)

        # Read in the 2D (nReads, bufferSize) readings array
        readings: np.ndarray = np.load(data_path)

        # Readings buffer needs to be 2D and the number of readings and timestamps must agree
        if readings.ndim != 2 or len(readings) != len(timestamps):
            raise ValueError(f"MS readings and timestamps disagree: {data_path}")

        # Possible to have an empty chunk, in that case just skip
        if len(readings) == 0:
            continue

        # Otherwise, append this chunk's readings to the buffer
        timestamp_buffer.append(timestamps)
        reading_buffer.append(readings)

    # Make sure we have some chunk at all
    if not reading_buffer:
        raise ValueError(f"No nonempty MS chunks found: {path_to_raw}")

    # Join nonempty chunks in recording order, keeping counts aligned with time.
    timestamps = np.concatenate(timestamp_buffer)
    readings = np.concatenate(reading_buffer)

    # Only provide the AS panel so the shared plotter skips TS, LS, and TEMP.
    ms_util.plot_chunk_readings(
        timestamps, readings, axes={"AS_ax": ax}, verbose=verbose,
    )

    return ax.figure, ax


def _plot_raw_eye_helper(
    recording_path: str | Path,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Plot raw eye summaries in a 6-by-1 layout.

    Args:
        recording_path: GKA/number recording directory. Raw and processed
            paths follow the FLIC_raw/FLIC_processing and sibling Neon/number layout.
        axes: 6 primary axes, or None to create a figure.
        verbose: Reserved for consistency with the other plotting helpers.

    Returns:
        Figure and primary axes array of shape (6,). Panels show pupil
        diameter, mean eyelid aperture, left-minus-right aperture, and gaze X/Y.
        Raw plots also show blink events. Gaze comes from the raw gaze.csv;
        all signals retain their native samples, outliers, and missing values.
        Times are Neon seconds relative to the first raw gaze sample.
        Does not show or save the figure.

    Raises:
        ValueError: Axes count, export selection, or timestamps are invalid.
        KeyError: A required CSV column is missing.
        OSError: A required recording file cannot be read.
    """
    # Check supplied axes before loading the recording.
    if axes is not None:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if axes.size != 6:
            raise ValueError(f"Raw eye plotting requires 6 axes; got {axes.size}." f" Recording: {recording_path}")

    # Locate the numbered recording's Neon export and, for cleaned data, its mapper output.
    recording_path = Path(recording_path)
    raw_recording_path: Path = Path(
        str(recording_path).replace("/FLIC_processing/", "/FLIC_raw/")
    )
    raw_neon_path: Path = (
        raw_recording_path.parent.parent / "Neon" / raw_recording_path.name
    )

    # Find one raw eye export for this recording.
    eye_state_paths: list[Path] = list(raw_neon_path.rglob("3d_eye_states.csv"))
    if len(eye_state_paths) != 1:
        raise ValueError(f"Expected one eye-state export in {raw_neon_path}; found {len(eye_state_paths)}.")

    eye_table: pd.DataFrame = pd.read_csv(eye_state_paths[0])
    eye_timestamps_ns: np.ndarray = eye_table["timestamp [ns]"].to_numpy()
    gaze_table: pd.DataFrame = pd.read_csv(eye_state_paths[0].parent / "gaze.csv")
    gaze_timestamps_ns: np.ndarray = gaze_table["timestamp [ns]"].to_numpy()

    # Validate the timelines before aligning their clocks.
    for name, timestamps in (
        ("Eye state", eye_timestamps_ns),
        ("Raw gaze", gaze_timestamps_ns),
    ):
        if (not len(timestamps) or not np.all(np.isfinite(timestamps))
                or np.any(timestamps[1:] <= timestamps[:-1])):
            raise ValueError(f"{name} timestamps must be finite and strictly increasing." f" Recording: {recording_path}")

    # Raw Neon signals share their own clock, independent of the GKA sensors.
    # Use the first gaze sample as zero for gaze, eye states, and blink events.
    neon_origin_ns: int | float = gaze_timestamps_ns[0]

    # Subtract the large epoch values before converting to seconds.
    eye_time_s: np.ndarray = (eye_timestamps_ns - neon_origin_ns).astype(float) * 1e-9
    gaze_time_s: np.ndarray = (gaze_timestamps_ns - neon_origin_ns).astype(float) * 1e-9
    panel_times_s: list[np.ndarray] = [eye_time_s] * 3 + [gaze_time_s] * 2
    panel_values: list[pd.Series] = [
        eye_table[["pupil diameter left [mm]", "pupil diameter right [mm]"]].mean(axis=1),
        eye_table[["eyelid aperture left [mm]", "eyelid aperture right [mm]"]].mean(axis=1),
        eye_table["eyelid aperture left [mm]"] - eye_table["eyelid aperture right [mm]"],
        gaze_table["gaze x [px]"], gaze_table["gaze y [px]"],
    ]

    blinks: pd.DataFrame = pd.read_csv(eye_state_paths[0].parent / "blinks.csv")
    blink_start_timestamps_ns: np.ndarray = blinks["start timestamp [ns]"].to_numpy()
    blink_end_timestamps_ns: np.ndarray = blinks["end timestamp [ns]"].to_numpy()
    if (not np.all(np.isfinite(blink_start_timestamps_ns)) or not np.all(np.isfinite(blink_end_timestamps_ns))
            or np.any(blink_end_timestamps_ns < blink_start_timestamps_ns)):
        raise ValueError("Blink events must have finite, ordered start/end times." f" Recording: {recording_path}")

    # Create the figure only after the input data have been validated.
    figure: Figure
    if axes is None:
        figure, axes = plt.subplots(6, 1, figsize=(14, 24), constrained_layout=True)
    else:
        figure = axes[0].figure

    titles: list[str] = ["Pupil Diameter", "Mean Eyelid Aperture (L/R Combined)",
              "Eyelid Aperture Difference (Left - Right)", "Gaze X", "Gaze Y"]
    labels: list[str] = ["Diameter [mm]", "Aperture [mm]", "Difference [mm]",
              "Position [px]", "Position [px]"]

    for axis, time, value, title, label in zip(axes, panel_times_s, panel_values, titles, labels):
        axis.plot(time, value, "k.", markersize=3, linewidth=1)
        axis.set_title(title, fontsize=12, fontweight="bold", pad=20)
        axis.set_ylabel(label)

    # Include every event edge, including blinks between gaze samples.
    blink_timestamps_ns: np.ndarray = np.unique(np.concatenate((gaze_timestamps_ns, blink_start_timestamps_ns, blink_end_timestamps_ns)))
    blink_active: np.ndarray = np.zeros(len(blink_timestamps_ns), dtype=float)
    for start, end in zip(blink_start_timestamps_ns, blink_end_timestamps_ns):
        blink_active[(blink_timestamps_ns >= start) & (blink_timestamps_ns < end)] = 1

    blink_time_s: np.ndarray = (blink_timestamps_ns - neon_origin_ns).astype(float) * 1e-9
    axes[5].step(blink_time_s, blink_active, where="post", color="black")
    axes[5].set_title("Raw Blink Events", fontsize=12, fontweight="bold", pad=20)
    axes[5].set_ylabel("Blink Active")
    axes[5].set_ylim(-0.1, 1.1)

    for axis in axes:
        axis.set_xlabel("Raw Neon time from first gaze [s]")
        axis.tick_params(axis="x", labelbottom=True)
        axis.grid(True, alpha=0.3)
        axis.margins(x=0)

    return figure, axes


def _plot_processed_eye_helper(
    recording_path: str | Path,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray, np.ndarray, np.ndarray]:
    """Plot cleaned eye summaries and build the shared world-to-Neon clock.

    Args:
        recording_path: GKA/number recording directory. Raw and processed
            paths follow the FLIC_raw/FLIC_processing and sibling Neon/number layout.
        axes: 5 primary axes, or None to create a figure.
        verbose: Reserved for consistency with the other plotting helpers.

    Returns:
        Figure, primary axes array of shape (5,), original world timestamps
        in seconds, and their paired Neon elapsed seconds. The caller uses
        these clock pairs to align world and MS plots. Panels show pupil
        diameter, mean eyelid aperture, left-minus-right aperture, and gaze X/Y.
        Only saved cleaned signals are plotted; no blink channel is inferred.
        Times are Neon seconds relative to the first mapped world frame.
        Does not show or save the figure.

    Raises:
        ValueError: Axes count, export selection, or timestamps are invalid.
        KeyError: A required CSV column is missing.
        OSError: A required recording file cannot be read.
    """
    # Check supplied axes before loading the recording.
    if axes is not None:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if axes.size != 5:
            raise ValueError(f"Cleaned eye plotting requires 5 axes; got {axes.size}." f" Recording: {recording_path}")

    recording_path = Path(recording_path)

    # Processed GKA and Neon attempts are siblings under the same activity.
    # Their raw inputs follow the same layout on the FLIC_raw volume.
    mapper_path: Path = (
        recording_path.parent.parent / "Neon" / recording_path.name / "egocentric_mapper_results"
    )
    raw_recording_path: Path = Path(
        str(recording_path).replace("/FLIC_processing/", "/FLIC_raw/")
    )
    raw_neon_path: Path = (
        raw_recording_path.parent.parent / "Neon" / raw_recording_path.name
    )

    # Each mapper row pairs a world-video frame with its Neon timestamp.
    mapped_world_timestamp_table: pd.DataFrame = pd.read_csv(mapper_path / "alternative_camera_timestamps.csv")
    
    # Subtract the shared epoch before converting to float seconds, preserving
    # precision. Every Neon timeline below uses this same world-frame origin.
    if mapped_world_timestamp_table.empty:
        raise ValueError(f"Mapper video timestamps are empty. Recording: {recording_path}")
    
    # Let's get the relative time of the recording frome the Neon and then convert it to seconds
    mapped_world_time_s: np.ndarray = (
        mapped_world_timestamp_table["timestamp [ns]"].to_numpy() - mapped_world_timestamp_table["timestamp [ns]"].iloc[0]
    ).astype(float) * 1e-9

    # Include missing-frame rows so the metadata follows the mapper video layout.
    # Inserted rows have timestamps but NaN gain settings.
    raw_world_metadata: pd.DataFrame = world_util.world_metadata_from_chunks(
        os.fspath(raw_recording_path),
        convert_to_seconds=True,
        verbose=verbose,
        fill_missing_frames=True,
    )

    # A one-to-one frame pairing is required to translate world/MS times into
    # Neon time. Keep inserted frames in both arrays to prevent clock shifts.
    assert len(raw_world_metadata) == len(mapped_world_time_s), (
        f"Filled world metadata does not match the mapper video frame layout "
        f"for recording: {recording_path}"
    )
    raw_world_time_s: np.ndarray = raw_world_metadata["timestamp"].to_numpy()

    # Interpolation requires finite, strictly increasing source and target clocks.
    for clock_name, timestamps in (
        ("World", raw_world_time_s),
        ("Mapper video", mapped_world_time_s),
    ):
        if (len(timestamps) < 2 or not np.all(np.isfinite(timestamps))
                or np.any(np.diff(timestamps) <= 0)):
            raise ValueError(
                f"{clock_name} clock needs at least two finite, increasing timestamps. "
                f"Recording: {recording_path}"
            )

    # Keep every timestamp pair, including inserted frames. Zero is the first
    # mapped world frame in Neon time, shared by eyes, world, and MS.
    neon_clock_s: np.ndarray = mapped_world_time_s

    cleaned_eye_table: pd.DataFrame = pd.read_csv(mapper_path / "cleaned_eye_data.csv")

    # The current cleaner writes time [s] from the first overlapping gaze
    # sample, without absolute timestamps. Recover its samples from the
    # source overlap and express them relative to the shared world origin.
    eye_state_export_paths: list[Path] = list(raw_neon_path.rglob("3d_eye_states.csv"))
    if len(eye_state_export_paths) != 1:
        raise ValueError(f"Expected one eye-state export in {raw_neon_path}; found {len(eye_state_export_paths)}.")

    eye_state_table: pd.DataFrame = pd.read_csv(eye_state_export_paths[0])
    eye_state_time_s: np.ndarray = (
        eye_state_table["timestamp [ns]"].to_numpy() - first_world_frame_neon_ns
    ).astype(float) * 1e-9
    mapped_gaze_table: pd.DataFrame = pd.read_csv(mapper_path / "alternative_camera_gaze.csv")
    mapped_gaze_time_s: np.ndarray = (
        mapped_gaze_table["timestamp [ns]"].to_numpy() - first_world_frame_neon_ns
    ).astype(float) * 1e-9

    # Validate the timelines before aligning their clocks.
    for name, timestamps in (("Eye state", eye_state_time_s), ("Mapped gaze", mapped_gaze_time_s)):
        if (not len(timestamps) or not np.all(np.isfinite(timestamps))
                or np.any(timestamps[1:] <= timestamps[:-1])):
            raise ValueError(f"{name} timestamps must be finite and strictly increasing." f" Recording: {recording_path}")

    # Reproduce the cleaner's sample selection to verify that this CSV belongs
    # to these source exports before recovering its absolute position in time.
    overlapping_gaze_time_s: np.ndarray = mapped_gaze_time_s[
        (mapped_gaze_time_s >= max(mapped_gaze_time_s[0], eye_state_time_s[0]))
        & (mapped_gaze_time_s <= min(mapped_gaze_time_s[-1], eye_state_time_s[-1]))
    ]
    # The cleaner converts epoch nanoseconds to float before subtraction;
    # allow sub-microsecond rounding when comparing its elapsed seconds.
    if (not len(overlapping_gaze_time_s) or len(overlapping_gaze_time_s) != len(cleaned_eye_table)
            or not np.allclose(cleaned_eye_table["time [s]"],
                               (overlapping_gaze_time_s - overlapping_gaze_time_s[0]),
                               rtol=0, atol=1e-6)):
        raise ValueError("Cleaned eye time does not match the source gaze/eye overlap." f" Recording: {recording_path}")

    # Plot every saved cleaned row. The overlap above identifies the CSV's
    # time origin; it does not trim the saved signals or zero each sensor
    # independently. Preserve the delay relative to the first world frame.
    # Do not clip to the world-video interval: negative times and eye samples
    # after the last world frame remain visible on the shared time axis.
    cleaned_time_s: np.ndarray = (
        cleaned_eye_table["time [s]"].to_numpy() + overlapping_gaze_time_s[0]
    )
    if not np.all(np.isfinite(cleaned_time_s)) or np.any(np.diff(cleaned_time_s) <= 0):
        raise ValueError(f"Cleaned eye times must be finite and strictly increasing. Recording: {recording_path}")
    panel_times_s: list[np.ndarray] = [cleaned_time_s] * 5
    panel_values: list[np.ndarray] = [cleaned_eye_table[column].to_numpy() for column in (
        "pupil diameter [mm]", "eyelid aperture [mm]",
        "eyelid aperture difference [mm]", "gaze x [px]", "gaze y [px]",
    )]

    # Create the figure only after the input data have been validated.
    figure: Figure
    if axes is None:
        figure, axes = plt.subplots(5, 1, figsize=(14, 20), constrained_layout=True)
    else:
        figure = axes[0].figure

    titles: list[str] = ["Pupil Diameter", "Mean Eyelid Aperture (L/R Combined)",
              "Eyelid Aperture Difference (Left - Right)", "Gaze X", "Gaze Y"]
    labels: list[str] = ["Diameter [mm]", "Aperture [mm]", "Difference [mm]",
              "Position [px]", "Position [px]"]

    for axis, time, value, title, label in zip(axes, panel_times_s, panel_values, titles, labels):
        axis.plot(time, value, "r-", markersize=3, linewidth=1)
        axis.set_title(title, fontsize=12, fontweight="bold", pad=20)
        axis.set_ylabel(label)

    for axis in axes:
        axis.set_xlabel("Neon time from first mapped world frame [s]")
        axis.tick_params(axis="x", labelbottom=True)
        axis.grid(True, alpha=0.3)
        axis.margins(x=0)

    return figure, axes, raw_world_time_s, neon_clock_s


def plot_raw_recording(
    path_to_raw: str | Path,
    verbose: bool = False,
    plot_title: str | None = None,
    time_range: tuple[int | float, int | float] | None = None,
) -> tuple[Figure, np.ndarray]:
    """Plot raw sensors and uncleaned eye summaries in a 9-by-1 figure.

    Args:
        path_to_raw: String or Path to raw world/MS chunks and config.pkl.
        verbose: Show loading and plotting progress bars.
        plot_title: Figure title; None uses the recording directory basename.
        time_range: (start, end) in each panel's displayed seconds, including start
            and excluding end. None plots all saved samples, including cleaned
            eye samples outside the world-video interval. Frame numbers
            retain their original positions in the recording. Raw eyes use seconds
            since the first raw gaze sample, independently of the GKA clock.

    Returns:
        Figure and primary axes array of shape (9,). Top: raw Bayer R, G, B
        and whole-frame means. Middle: analog/digital gains on the left and
        exposure in microseconds on the right, with a combined legend.
        Third: AS spectral channel counts with a channel legend. All tiles
        use time in seconds; world tiles also show frame numbers above.
        The exposure twin axis is excluded from the array. The figure is
        neither shown nor saved, and is closed if plotting fails.

    Raises:
        ValueError: Data are empty, malformed, or misaligned.
        KeyError: Required recording configuration fields are missing.
        OSError: The recording directory, chunks, or config.pkl cannot be read.
    """

    # Time bounds use the displayed clock described above.
    if time_range is not None:
        assert len(time_range) == 2, f"time_range must contain (start, end); got {time_range}." f" Recording: {path_to_raw}"
        assert all(isinstance(value, (int, float, np.integer, np.floating)) for value in time_range), (
            f"time_range bounds must be numbers; got {time_range}." f" Recording: {path_to_raw}"
        )
        assert np.all(np.isfinite(time_range)), f"time_range bounds must be finite; got {time_range}." f" Recording: {path_to_raw}"
        assert time_range[0] < time_range[1], f"time_range requires start < end; got {time_range}." f" Recording: {path_to_raw}"

    path_to_raw = os.path.abspath(os.path.expanduser(path_to_raw))

    # Reserve two world panels, one MS panel, and six eye panels.
    # Keep independent horizontal axes because raw sensor clocks are not aligned.
    figure: Figure
    axes: np.ndarray
    figure, axes = plt.subplots(
        9, 1, figsize=(14, 36), constrained_layout=True,
    )

    # As in the processed plotter, load eyes, world, and MS through their helpers.
    # Raw eye data need no world-to-Neon clock conversion.
    try:
        _plot_raw_eye_helper(path_to_raw, axes[3:], verbose=verbose)
        _plot_raw_world_helper(path_to_raw, axes[:2], verbose=verbose)
        _plot_raw_ms_helper(path_to_raw, axes[2], verbose=verbose)

        for axis in axes:
            axis.set_xlabel("GKA time [s]" if axis in axes[:3] else "Raw Neon time from first gaze [s]")
            axis.tick_params(axis="x", labelbottom=True)
            axis.set_title(axis.get_title(), fontsize=12, fontweight="bold", pad=20)
            for child in axis.child_axes:
                if isinstance(child, SecondaryAxis):
                    child.set_xlabel(child.get_xlabel(), fontsize=9, labelpad=3)

    # If an error happens, close the figure
    except Exception:
        plt.close(figure)
        raise

    # Use the recording name unless the caller provides a custom title.
    if plot_title is None:
        plot_title = os.path.basename(path_to_raw)
    figure.suptitle(plot_title)

    # Apply the selected time display window after all panels and labels exist.
    # This changes the plotted data only; the raw sensor files remain unchanged.
    if time_range is not None:
        # Include samples at the start time, but exclude samples at the end time.
        start_time: int | float = time_range[0]
        end_time: int | float = time_range[1]

        # `axes` contains the primary time-series panels. `figure.axes` also
        # contains the exposure twin axis created by the camera-settings helper.
        # Include primary panels and any twin axes sharing their time coordinates
        # so exposure is masked along with the gain and raw sensor traces.
        time_axes: list[Axes] = [
            axis for axis in figure.axes
            if any(
                axis is primary or axis.get_shared_x_axes().joined(axis, primary)
                for primary in axes.flat
            )
        ]

        axis: Axes
        line: Line2D
        for axis in time_axes:
            # A panel may contain multiple lines, each with its own sensor timestamps.
            # For example, world-camera frames and MS readings have different sampling rates.
            for line in axis.lines:
                line_time: np.ndarray = np.asarray(line.get_xdata(), dtype=float)

                # Copy into floating-point storage so NaNs can be assigned safely,
                # even when the original plotted values were integer counts.
                line_values: np.ndarray = np.asarray(line.get_ydata(), dtype=float).copy()
                outside_range: np.ndarray = (line_time < start_time) | (line_time >= end_time)

                # Hide values outside the window without deleting their timestamps.
                # Keeping the time coordinates preserves the existing frame-number mapping.
                line_values[outside_range] = np.nan
                line.set_ydata(line_values)

            # Recalculate data limits from the modified lines. Autoscaled value axes
            # then fit the remaining finite samples; explicitly fixed scales stay fixed.
            axis.relim()
            axis.autoscale_view()

            # Set the requested horizontal window explicitly, including intervals
            # near the boundaries that contain no recorded samples.
            axis.set_xlim(start_time, end_time)

    return figure, axes


def _interpolate(values: np.ndarray, source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Interpolate coordinates and extend endpoint slopes into plot margins.

    Args:
        values: Coordinates to transform.
        source: Increasing source coordinates with at least two samples.
        target: Corresponding target coordinates, the same length as source.

    Returns:
        Array with the shape of values, linearly interpolated inside source
        bounds and linearly extrapolated using the endpoint slopes outside.
    """
    result: np.ndarray = np.interp(values, source, target)
    result = np.where(
        values < source[0],
        target[0] + (values - source[0]) * (target[1] - target[0]) / (source[1] - source[0]),
        result,
    )
    return np.where(
        values > source[-1],
        target[-1] + (values - source[-1]) * (target[-1] - target[-2]) / (source[-1] - source[-2]),
        result,
    )


def _plot_processed_world_helper(
    world_path: str | Path,
    axes: np.ndarray | None = None,
    verbose: bool = False,
    *,
    GKA_timestamps_raw_s: np.ndarray | None = None,
    neon_clock_s: np.ndarray | None = None,
) -> tuple[Figure, np.ndarray]:
    """Plot processed world data in a 4-by-1 layout.

    Args:
        world_path: String or Path to W, containing transformed.hdf5.
        axes: Four primary axes in panel order, or None to create a 4-by-1 figure.
        verbose: Show loading and averaging progress bars.


        GKA_timestamps_raw_s: Original world timestamps paired with neon_clock_s.
        neon_clock_s: Corresponding Neon elapsed seconds. Both arrays must be
            supplied together. Linear interpolation preserves individual sensor
            samples; endpoint slopes extend the clock beyond camera coverage.

    Returns:
        Figure and four primary axes: L+M, L-M, S-(L+M), and saved scene
        illuminance in lux. Each opponent channel has its own scale. Channel
        means are computed together in one pass and stored in an (N, 3)
        matrix before plotting. Twin axes are excluded. Does not show or save.

    Raises:
        AssertionError: Timestamp or frame shapes do not match the
            expected dimensions or frame counts.
        ValueError: Axes count or timestamp values are invalid.
        OSError: The HDF5 file or camera calibration cannot be read.
        KeyError: A required dataset is missing.
    """

    # A clock correspondence must contain matching, increasing timestamp pairs.
    if (GKA_timestamps_raw_s is None) != (neon_clock_s is None):
        raise ValueError("Supply both GKA_timestamps_raw_s and neon_clock_s." f" Recording: {Path(world_path).parent}")

    if GKA_timestamps_raw_s is not None:
        if (len(GKA_timestamps_raw_s) < 2 or len(GKA_timestamps_raw_s) != len(neon_clock_s)
                or not np.all(np.isfinite(GKA_timestamps_raw_s))
                or not np.all(np.isfinite(neon_clock_s))
                or np.any(np.diff(GKA_timestamps_raw_s) <= 0)
                or np.any(np.diff(neon_clock_s) <= 0)):
            raise ValueError("Clock pairs must be finite, increasing, and equal in length." f" Recording: {Path(world_path).parent}")

    # Validate supplied axes before setting up the figure.
    if axes is not None:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        assert len(axes) == 4, f"World plotting requires four axes; got {len(axes)}." f" Recording: {Path(world_path).parent}"

    # Create a local figure only when the caller has not supplied axes.
    figure: Figure
    if axes is None:
        figure, axes = plt.subplots(
            4, 1, figsize=(12, 14), sharex=True, constrained_layout=True
        )
    else:
        figure = axes[0].figure

    # First, let's open the hdf5 file
    with h5py.File(os.path.join(world_path, "transformed.hdf5"), "r") as hdf5_filepointer:
        # Let's first get the sensor timestamps, which are stored as nanoseconds.
        # These are all read into memory
        time: np.ndarray = hdf5_filepointer["agcSettings/timestamps"][:] * WORLD_NANOSECONDS_TO_SECONDS

        # Convert the camera clock before creating traces and frame-number axes.
        if GKA_timestamps_raw_s is not None:
            time = _interpolate(time, GKA_timestamps_raw_s, neon_clock_s)

        # Scene illuminance is stored in lux, with one value per world frame.
        scene_illuminance: np.ndarray = hdf5_filepointer["data/scene_illuminance"][:]
        assert scene_illuminance.shape == time.shape, (
            "data/scene_illuminance must have one value per world timestamp; "
            f"expected {time.shape}, got {scene_illuminance.shape}."
            f" Recording: {Path(world_path).parent}"
        )

        # Then, let's get the frames datagroup. These are not yet read into memory
        frames: h5py.Dataset = hdf5_filepointer["data/frames"]

        # Validate dimensions before indexing shapes or comparing frame counts.
        assert time.ndim == 1, (
            "agcSettings/timestamps must have shape (nFrames,); "
            f"got shape {time.shape} with {time.ndim} dimensions."
            f" Recording: {Path(world_path).parent}"
        )
        assert frames.ndim == 4, (
            "data/frames must have shape (nFrames, azimuth, elevation, channel); "
            f"got shape {frames.shape} with {frames.ndim} dimensions."
            f" Recording: {Path(world_path).parent}"
        )

        # Each image must have one timestamp and the three opponent channels.
        assert frames.shape[0] == len(time), (
            "data/frames and agcSettings/timestamps must have matching frame counts; "
            f"got {frames.shape[0]} frames and {len(time)} timestamps."
            f" Recording: {Path(world_path).parent}"
        )
        assert frames.shape[-1] == 3, (
            "data/frames must contain 3 channels (L+M, L-M, S-(L+M)); "
            f"got {frames.shape[-1]} channels in shape {frames.shape}."
            f" Recording: {Path(world_path).parent}"
        )

        # Define the size of chunks of frames we will read into memory
        chunk_size: int = 1000

        # Calculate the number of chunks the dataset contains
        num_chunks: int = (len(time) + chunk_size - 1) // chunk_size


        opponent_channel_means: np.ndarray = np.full((len(time), 3), np.nan)

        # Iterate over the chunks in the dataset
        for chunk_index in tqdm(
            range(num_chunks), desc="Loading world chunks", disable=not verbose,
        ):
            # Get the start stop of this chunk in the dataset, then load it into a
            # buffer in memory
            start: int = chunk_index * chunk_size
            stop: int = min(start + chunk_size, len(time))
            buffer: np.ndarray = frames[start:stop]

            # We are interested in storing the average opponent channel and plotting these
            # so average the finite locations in this frame and save it into the buffer
            for frame_index in tqdm(
                range(len(buffer)), desc="Averaging opponent frames",
                disable=not verbose, leave=False,
            ):
                frame: np.ndarray = buffer[frame_index]
                output_frame_index: int = start + frame_index

                # Average L+M, L-M, and S-(L+M) separately over their valid pixels.
                for channel_index in range(frame.shape[-1]):
                    channel: np.ndarray = frame[:, :, channel_index]
                    valid_pixels: np.ndarray = channel[np.isfinite(channel)]

                    # Leave the initialized NaN when this channel has no valid pixels.
                    if valid_pixels.size == 0:
                        continue

                    opponent_channel_means[output_frame_index, channel_index] = valid_pixels.mean()

    # First three panels: separate spatial means of each opponent channel.
    channel_labels: tuple[str, ...] = ("L+M", "L-M", "S-(L+M)")
    channel_colors: tuple[str, ...] = ("black", "darkorange", "purple")
    for channel_index, (label, color) in enumerate(zip(channel_labels, channel_colors)):
        axes[channel_index].plot(
            time, opponent_channel_means[:, channel_index], color=color,
        )
        axes[channel_index].set_title(f"Average {label} radiance")
        axes[channel_index].set_ylabel(f"Mean {label}")

    # Fourth panel: scene illuminance computed from the reconstructed RGB radiance.
    axes[3].plot(time, scene_illuminance, color="gray", label="Camera photopic illuminance")
    axes[3].set_title("Camera illuminance")
    axes[3].set_ylabel("Illuminance [lux]")

    # Add stored-frame coordinates above the world panels.
    for axis in axes:
        if not axis.get_visible():
            continue
        axis.set_xlabel("time [s]")
        axis.tick_params(axis="x", labelbottom=True)
        axis.set_title(axis.get_title(), fontsize=12, fontweight="bold", pad=20)
        axis.grid(True, alpha=0.25)
        frame_timestamps: np.ndarray = np.asarray(time, dtype=float)
        if (frame_timestamps.ndim != 1 or frame_timestamps.size == 0
                or not np.all(np.isfinite(frame_timestamps)) or np.any(np.diff(frame_timestamps) <= 0)):
            raise ValueError("Frame timestamps must be finite, nonempty, and increasing." f" Recording: {Path(world_path).parent}")

        # Interpolate frame positions so dropped frames do not distort the time scale.
        frame_indices: np.ndarray = np.arange(len(frame_timestamps), dtype=float)
        if len(frame_timestamps) == 1:
            frame_timestamps = np.array([frame_timestamps[0], frame_timestamps[0] + 1])
            frame_indices = np.array([0., 1.])

        frame_axis: SecondaryAxis = axis.secondary_xaxis(
            "top", functions=(
                lambda times: _interpolate(times, frame_timestamps, frame_indices),
                lambda frames: _interpolate(frames, frame_indices, frame_timestamps),
            ),
        )
        frame_axis.set_xlabel("Frame number", fontsize=9, labelpad=3)
        frame_axis.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))

    return figure, axes


def _plot_processed_ms_helper(
    ms_path: str | Path,
    axes: np.ndarray | None = None,
    verbose: bool = False,
    *,
    GKA_timestamps_raw_s: np.ndarray | None = None,
    neon_clock_s: np.ndarray | None = None,
) -> tuple[Figure, np.ndarray]:
    """Plot processed MS spectra and fit values in a 2-by-1 layout.

    Args:
        ms_path: String or Path to M, containing transformed.hdf5.
        axes: Two primary axes in panel order, or None to create a 2-by-1 figure.
        verbose: Reserved for consistency with the other plotting helpers.

        GKA_timestamps_raw_s: Original world timestamps paired with neon_clock_s.
        neon_clock_s: Corresponding Neon elapsed seconds. Both arrays must be
            supplied together. Linear interpolation preserves individual sensor
            samples; endpoint slopes extend the clock beyond camera coverage.

    Returns:
        Figure and primary axes array of shape (2,). Tile 0 is a spectral
        radiance heatmap with time horizontally, wavelength vertically
        (400-900 nm), and a radiance colorbar. Tile 1 shows the overall fit
        value divided by the per-reading mean calibrated MS measurement over
        time. The MS calibration offset is applied in seconds before the optional
        conversion to Neon time. The colorbar axis is excluded
        from the returned array. Does not show or save the figure.

    Raises:
        AssertionError: Axes count, timestamps, wavelength sampling, or dataset
            shapes are incompatible.
        OSError: The HDF5 file cannot be opened.
        KeyError: A required dataset is missing.
    """

    # A clock correspondence must contain matching, increasing timestamp pairs.
    if (GKA_timestamps_raw_s is None) != (neon_clock_s is None):
        raise ValueError("Supply both GKA_timestamps_raw_s and neon_clock_s." f" Recording: {Path(ms_path).parent}")

    if GKA_timestamps_raw_s is not None:
        if (len(GKA_timestamps_raw_s) < 2 or len(GKA_timestamps_raw_s) != len(neon_clock_s)
                or not np.all(np.isfinite(GKA_timestamps_raw_s))
                or not np.all(np.isfinite(neon_clock_s))
                or np.any(np.diff(GKA_timestamps_raw_s) <= 0)
                or np.any(np.diff(neon_clock_s) <= 0)):
            raise ValueError("Clock pairs must be finite, increasing, and equal in length." f" Recording: {Path(ms_path).parent}")

    # Validate supplied axes before setting up the figure.
    if axes is not None:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        assert len(axes) == 2, (
            "MS plotting requires two axes; "
            f"got {len(axes)} axes after flattening."
            f" Recording: {Path(ms_path).parent}"
        )

    # Create a local figure only when the caller has not supplied axes.
    figure: Figure
    if axes is None:
        figure, axes = plt.subplots(
            2, 1, figsize=(12, 7), sharex=True, constrained_layout=True
        )
    else:
        figure = axes[0].figure

    # Open the processed MS file and load one spectrum and fit value per reading.
    with h5py.File(
        os.path.join(ms_path, "transformed.hdf5"), "r"
    ) as hdf5_file:
        # MS timestamps are already stored in seconds.
        timestamps: np.ndarray = hdf5_file["timestamps"][:]

        # Each radiance row contains the wavelength samples for one timestamp.
        spectral_radiance: np.ndarray = hdf5_file["data/radiance"][:]

        # S defines the wavelength grid; fVal stores the fit residual norms.
        wavelength_sampling: np.ndarray = hdf5_file["data/S"][:].reshape(-1)
        fit_residual_norms: np.ndarray = hdf5_file["data/fVal"][:]

        # Use the original counts and calibration saved with this fit.
        raw_counts: np.ndarray = hdf5_file["data/counts/readings"][:]
        radiance_weights: np.ndarray = hdf5_file["metadata/minispect_radiance_weights"][:]

    # MS timestamps are seconds; the calibration constant is milliseconds.
    # The transformed writer stores uncorrected times, so apply this once here.
    timestamps = timestamps + ms_util.MS_TIME_OFFSET / 1000.0

    if GKA_timestamps_raw_s is not None:
        timestamps = _interpolate(timestamps, GKA_timestamps_raw_s, neon_clock_s)

    # Validate timestamps before computing intervals or checking aligned shapes.
    assert timestamps.ndim == 1, (
        "MS timestamps must have shape (nReadings,); "
        f"got shape {timestamps.shape} with {timestamps.ndim} dimensions."
        f" Recording: {Path(ms_path).parent}"
    )
    assert timestamps.size > 0, (
        f"MS timestamps must be nonempty; got shape {timestamps.shape}." f" Recording: {Path(ms_path).parent}"
    )
    assert np.all(np.isfinite(timestamps)), (
        "MS timestamps must be finite; "
        f"nonfinite values at indices {np.flatnonzero(~np.isfinite(timestamps))[:10]}."
        f" Recording: {Path(ms_path).parent}"
    )

    timestamp_intervals: np.ndarray = np.diff(timestamps)
    assert np.all(timestamp_intervals > 0), (
        "MS timestamps must be strictly increasing; "
        f"nonpositive intervals start at indices "
        f"{np.flatnonzero(timestamp_intervals <= 0)[:10]}, with differences "
        f"{timestamp_intervals[timestamp_intervals <= 0][:10]}."
        f" Recording: {Path(ms_path).parent}"
    )

    # Validate the wavelength grid before indexing it or converting its count.
    assert wavelength_sampling.size == 3, (
        "MS wavelength sampling must contain [start, step, count]; "
        f"got {wavelength_sampling.size} values: {wavelength_sampling}."
        f" Recording: {Path(ms_path).parent}"
    )
    assert np.all(np.isfinite(wavelength_sampling)), (
        "MS wavelength sampling must be finite; "
        f"got {wavelength_sampling}."
        f" Recording: {Path(ms_path).parent}"
    )
    assert wavelength_sampling[1] > 0, (
        "MS wavelength step must be positive; "
        f"got {wavelength_sampling[1]} nm."
        f" Recording: {Path(ms_path).parent}"
    )
    assert wavelength_sampling[2] > 0, (
        "MS wavelength sample count must be positive; "
        f"got {wavelength_sampling[2]}."
        f" Recording: {Path(ms_path).parent}"
    )
    assert wavelength_sampling[2] == int(wavelength_sampling[2]), (
        "MS wavelength sample count must be an integer; "
        f"got {wavelength_sampling[2]}."
        f" Recording: {Path(ms_path).parent}"
    )

    # Require one spectrum and one fit residual norm per timestamp.
    assert spectral_radiance.shape == (len(timestamps), int(wavelength_sampling[2])), (
        "MS spectral radiance must have shape (nReadings, nWavelengths); "
        f"expected {(len(timestamps), int(wavelength_sampling[2]))}, "
        f"got {spectral_radiance.shape}."
        f" Recording: {Path(ms_path).parent}"
    )
    assert fit_residual_norms.shape == (len(timestamps),), (
        "MS fit residual norms must have shape (nReadings,); "
        f"expected {(len(timestamps),)}, got {fit_residual_norms.shape}."
        f" Recording: {Path(ms_path).parent}"
    )

    # The joint fit uses seven MS channels, calibrated into radiance units.
    assert raw_counts.ndim == 2, (
        f"MS counts must have shape (nReadings, channels); got {raw_counts.shape}." f" Recording: {Path(ms_path).parent}"
    )
    assert raw_counts.shape[0] == len(timestamps), (
        f"MS counts must have {len(timestamps)} rows; got {raw_counts.shape[0]}." f" Recording: {Path(ms_path).parent}"
    )
    assert radiance_weights.shape == (7,), (
        f"Expected seven saved MS calibration weights; got {radiance_weights.shape}." f" Recording: {Path(ms_path).parent}"
    )
    assert raw_counts.shape[1] >= len(radiance_weights), (
        f"MS counts need at least {len(radiance_weights)} channels; got {raw_counts.shape[1]}." f" Recording: {Path(ms_path).parent}"
    )
    assert np.all(np.isfinite(radiance_weights)), (
        f"MS calibration weights must be finite; got {radiance_weights}." f" Recording: {Path(ms_path).parent}"
    )

    # Match the count-to-radiance conversion used by the joint estimator.
    calibrated_ms_measurements: np.ndarray = (
        raw_counts[:, :len(radiance_weights)] / 10 ** radiance_weights
    )
    mean_ms_measurements: np.ndarray = calibrated_ms_measurements.mean(axis=1)

    # Normalize each residual independently; undefined ratios remain missing.
    valid_readings: np.ndarray = (
        np.isfinite(mean_ms_measurements)
        & (mean_ms_measurements > 0)
        & np.isfinite(fit_residual_norms)
    )
    normalized_fit_values: np.ndarray = np.full(fit_residual_norms.shape, np.nan)
    normalized_fit_values[valid_readings] = (
        fit_residual_norms[valid_readings] / mean_ms_measurements[valid_readings]
    )

    # S stores the starting wavelength, step size, and sample count.
    wavelengths: np.ndarray = (
        wavelength_sampling[0]
        + wavelength_sampling[1] * np.arange(int(wavelength_sampling[2]))
    )

    # Top panel: spectral radiance across wavelength and time.
    # Transpose spectra so rows represent wavelength and columns represent time.
    radiance_mesh: QuadMesh = axes[0].pcolormesh(
        timestamps,
        wavelengths,
        np.ma.masked_invalid(spectral_radiance.T),
        shading="auto",
        cmap="viridis",
    )
    figure.colorbar(radiance_mesh, ax=axes[0], label="Spectral radiance")
    axes[0].set_title("Spectral radiance over time")
    axes[0].set_ylabel("Wavelength [nm]")
    axes[0].set_ylim(400, 900)

    # Bottom panel: fit residual relative to mean calibrated MS signal over time.
    axes[1].plot(timestamps, normalized_fit_values, color="tab:blue")
    axes[1].set_title("Normalized fit value")
    axes[1].set_ylabel("Residual norm / mean MS measurement")

    # Label both time axes without reporting progress for routine formatting.
    axis: Axes
    for axis in axes:
        axis.set_xlabel("time [s]")
        axis.tick_params(axis="x", labelbottom=True)
        axis.set_title(axis.get_title(), fontsize=12, fontweight="bold", pad=20)

    # Add a grid to the fit panel; the heatmap already fills the spectral panel.
    axes[1].grid(True, alpha=0.25)

    return figure, axes


def plot_processed_recording(
    recording_directory: str | Path,
    verbose: bool = False,
    *,
    plot_title: str | None = None,
    time_range: tuple[int | float, int | float] | None = None,
) -> tuple[Figure, np.ndarray]:
    """Plot processed world, MS, and cleaned eye time series.

    Args:
        recording_directory: String or Path to the parent directory containing
            W/transformed.hdf5 and M/transformed.hdf5.
        verbose: Show loading and plotting progress bars.
        plot_title: Figure title; None uses the recording directory basename.
        time_range: (start, end) in Neon seconds since the first mapped world
            frame, including start
            and excluding end. None plots all saved samples, including cleaned
            eye samples outside the world-video interval. Frame numbers
            retain their original positions in the recording.


    Returns:
        Figure and primary axes array of shape (6, 2).
        The final five panels contain cleaned eye summaries; the last tile is hidden.
        First three rows, left column: L+M,
        L-M, and S-(L+M) on the left; MS spectral radiance, normalized fit
        value, and photopic/melanopic illuminance on the right. Every subplot
        displays its own time labels. The colorbar axis is excluded.
        Does not show or save; closes the figure on failure.

    Raises:
        OSError: A sensor file or calibration cannot be read.
        KeyError: A required dataset is missing.
        AssertionError: Dataset shapes or timestamps are incompatible.
    """

    # Time bounds use the displayed clock described above.
    if time_range is not None:
        assert len(time_range) == 2, f"time_range must contain (start, end); got {time_range}." f" Recording: {recording_directory}"
        assert all(isinstance(value, (int, float, np.integer, np.floating)) for value in time_range), (
            f"time_range bounds must be numbers; got {time_range}." f" Recording: {recording_directory}"
        )
        assert np.all(np.isfinite(time_range)), f"time_range bounds must be finite; got {time_range}." f" Recording: {recording_directory}"
        assert time_range[0] < time_range[1], f"time_range requires start < end; got {time_range}." f" Recording: {recording_directory}"

    # First, we will get the recording directory that was input
    recording_directory = os.path.abspath(os.path.expanduser(recording_directory))

    GKA_timestamps_raw_s: np.ndarray
    neon_clock_s: np.ndarray

    # Place world and MS panels side by side on one shared horizontal axis.
    figure: Figure
    axes: np.ndarray
    figure, axes = plt.subplots(
        6, 2, figsize=(18, 28.2),
        sharex=True, constrained_layout=True,
    )

    try:
        # We first load and plot the eye data. This is so that
        # we have a common clock for all of the other sensors and do not need
        # to update it later
        _, _, GKA_timestamps_raw_s, neon_clock_s = _plot_processed_eye_helper(
            recording_directory, axes[3:].flat[:5], verbose=verbose,
        )
        axes[-1, -1].set_visible(False)

        # First, let's plot the world plots along the first column
        _plot_processed_world_helper(
            os.path.join(recording_directory, "W"),
            np.array([*axes[:3, 0], axes[2, 1]], dtype=object),
            verbose=verbose, GKA_timestamps_raw_s=GKA_timestamps_raw_s, neon_clock_s=neon_clock_s,
        )

        # Next, let's plot the MS plots along the second column
        _plot_processed_ms_helper(
            os.path.join(recording_directory, "M"), axes[:2, 1],
            verbose=verbose, GKA_timestamps_raw_s=GKA_timestamps_raw_s, neon_clock_s=neon_clock_s,
        )

        # Compare saved illuminance values using each sensor's own timestamps.
        with h5py.File(
            os.path.join(recording_directory, "M", "transformed.hdf5"), "r"
        ) as hdf5_file:
            ms_timestamps: np.ndarray = hdf5_file["timestamps"][:]
            ms_illuminance: np.ndarray = hdf5_file["data/illuminance/lux"][:]
            ms_melanopic_edi: np.ndarray = hdf5_file["data/illuminance/melanopic_edi"][:]

        # Use exactly the same correction and clock mapping as the MS panels.
        ms_timestamps = ms_timestamps + ms_util.MS_TIME_OFFSET / 1000.0
        ms_timestamps = _interpolate(ms_timestamps, GKA_timestamps_raw_s, neon_clock_s)

        axes[2, 1].plot(
            ms_timestamps, ms_illuminance,
            color="black", label="MS photopic illuminance",
        )
        axes[2, 1].plot(
            ms_timestamps, ms_melanopic_edi,
            color="blue", label="MS melanopic EDI",
        )
        axes[2, 1].set_title("Photopic illuminance and melanopic EDI", fontsize=12, fontweight="bold", pad=20)
        axes[2, 1].legend(loc="upper right")

        for axis in axes.flat:
            if axis.get_visible():
                axis.set_xlabel("Neon time from first mapped world frame [s]")

    # If any error occurs, we need to close the figure
    except Exception:
        # Avoid leaving a partially populated figure open if loading fails.
        plt.close(figure)
        raise

    # Use the recording name unless the caller provides a custom title.
    if plot_title is None:
        plot_title = os.path.basename(recording_directory)
    figure.suptitle(plot_title)

    # Apply the selected time display window after all panels and labels exist.
    # This changes the plotted data only; the saved sensor datasets remain unchanged.
    if time_range is not None:
        # Include samples at the start time, but exclude samples at the end time.
        start_time: int | float = time_range[0]
        end_time: int | float = time_range[1]

        # Select the time-series panels and exclude the spectral colorbar,
        # whose coordinates represent radiance rather than recording time.
        time_axes: list[Axes] = [
            axis for axis in figure.axes
            if any(
                axis is primary or axis.get_shared_x_axes().joined(axis, primary)
                for primary in axes.flat
            )
        ]

        axis: Axes
        line: Line2D
        for axis in time_axes:
            # A panel may contain multiple lines, each with its own sensor timestamps.
            # The world and MS panels have different sampling rates.
            for line in axis.lines:
                line_time: np.ndarray = np.asarray(line.get_xdata(), dtype=float)

                # Copy into floating-point storage so NaNs can be assigned safely,
                # even when the original plotted values were integer counts.
                line_values: np.ndarray = np.asarray(line.get_ydata(), dtype=float).copy()
                outside_range: np.ndarray = (line_time < start_time) | (line_time >= end_time)

                # Hide values outside the window without deleting their timestamps.
                # Keeping the time coordinates preserves the existing frame-number mapping.
                line_values[outside_range] = np.nan
                line.set_ydata(line_values)

            # Recalculate data limits from the modified lines. Autoscaled value axes
            # then fit the remaining finite samples; explicitly fixed scales stay fixed.
            axis.relim()
            axis.autoscale_view()

            # Set the requested horizontal window explicitly, including intervals
            # near the boundaries that contain no recorded samples.
            axis.set_xlim(start_time, end_time)

        # The spectral heatmap is a QuadMesh, not a line, so the loop above did not
        # mask its values. Each heatmap column corresponds to one MS reading.
        # The fit-value line retains the same complete MS timestamp array.
        ms_time: np.ndarray = np.asarray(axes[1, 1].lines[0].get_xdata(), dtype=float)
        outside_ms_range: np.ndarray = (ms_time < start_time) | (ms_time >= end_time)
        spectral_mesh: QuadMesh = axes[0, 1].collections[0]

        # Preserve existing missing spectra as NaNs and work on a separate copy.
        # Restore wavelength-by-time shape if Matplotlib exposes a flattened array.
        spectral_values: np.ndarray = np.asarray(spectral_mesh.get_array().filled(np.nan)).copy()
        spectral_values = spectral_values.reshape(-1, len(ms_time))

        # Hide every wavelength for excluded readings, using their timestamps rather
        # than the inferred left/right edges of the heatmap cells. The color scale stays unchanged.
        spectral_values[:, outside_ms_range] = np.nan
        spectral_mesh.set_array(np.ma.masked_invalid(spectral_values))

    return figure, axes


def animate_processed_recording(
    recording_directory: str | Path,
    verbose: bool = False,
    *,
    plot_title: str | None = None,
    time_range: tuple[int | float, int | float] | None = None,
    output_path: str | Path,
    gif_fps: float = 5.0,
    dpi: int = 60,
    overwrite_existing: bool = False,
) -> Path:
    """Save an L+M camera video above the processed plots as a lossless AVI.

    Eye data are always included. All sensor paths are derived from the recording directory.
    Reuse plot_processed_recording to draw the selected interval, then move red
    time bars across its panels alongside the corresponding camera frames.

    Args:
        recording_directory: String or Path to the parent directory containing
            W/transformed.hdf5 and M/transformed.hdf5.
        verbose: Show loading, plotting, and animation-rendering progress bars.
        plot_title: Figure title; None uses the recording directory basename.
        time_range: (start, end) on the same clock as plot_processed_recording,
            including start
            and excluding end. None animates the complete camera recording.
        output_path: Required destination string or Path ending in .avi.
            Parent directories are created if needed.
        gif_fps: Positive output sampling and playback rate in frames per second.
            Defaults to 5.0; playback preserves recording time.
        dpi: Positive integer resolution for rendering the processed plot.
            Defaults to 60 dots per inch.
        overwrite_existing: Allow replacing an existing output video.
            Defaults to False.

    Returns:
        Absolute Path to the saved AVI video.

    Raises:
        AssertionError: Animation settings, time bounds, or camera data are invalid.
        FileExistsError: The destination exists and overwrite_existing is False.
        OSError: Input data cannot be read or output files cannot be written.
        KeyError: A required sensor dataset is missing.
        RuntimeError: The video encoder cannot open its output file.

    Notes:
        The plots describe the selected interval while the red time bars move.
        Camera timestamp gaps appear black instead
        of holding an outdated frame. L+M display scaling uses a one-second
        source window, including neighboring frames outside the selected interval.

        Rendered PNGs are staged in a temporary directory and encoded through
        video_io.dir_to_video. Memory use does not grow with animation duration,
        but temporary disk usage does. The figure and temporary frames are cleaned
        up on completion or failure. The finished video is moved to the requested
        destination after encoding.
    """
    assert np.isfinite(gif_fps), f"Animation gif_fps must be finite; got {gif_fps}." f" Recording: {recording_directory}"
    assert gif_fps > 0, f"Animation gif_fps must be positive; got {gif_fps}." f" Recording: {recording_directory}"
    assert isinstance(dpi, int), f"dpi must be an integer; got {dpi}." f" Recording: {recording_directory}"
    assert dpi > 0, f"dpi must be positive; got {dpi}." f" Recording: {recording_directory}"

    # Resolve paths once so later file operations do not depend on the working directory.
    recording_directory = Path(recording_directory).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()

    # Ensure that the output path is correctly a .avi path and that if it exists,
    # we will not overwrite it unless we intended to
    assert output_path.suffix.lower() == ".avi", f"Output must be an AVI file; got {output_path}."
    if output_path.exists() and not overwrite_existing:
        raise FileExistsError(f"Animation already exists: {output_path}")

    # Build the figure once without the video component for this time range
    figure: Figure
    axes: np.ndarray
    figure, axes = plot_processed_recording(
        recording_directory, verbose=verbose, plot_title=plot_title, time_range=time_range,
    )

    # Next, we will generate the animation on top of this plotted basis
    try:

        # Next, we will set the DPI of the figure, this lowers quality but increases speed
        figure.set_dpi(dpi)

        # Render directly to pixels without opening an interactive plotting window.
        canvas: FigureCanvasAgg = FigureCanvasAgg(figure)

        # Camera plot timestamps correspond directly to saved opponent-map frames.
        camera_time: np.ndarray = np.asarray(axes[2, 1].lines[0].get_xdata(), dtype=float)
        assert camera_time.size > 0, "Animation requires at least one camera frame." f" Recording: {recording_directory}"
        assert np.all(np.isfinite(camera_time)), "Camera timestamps contain nonfinite values." f" Recording: {recording_directory}"
        assert np.all(np.diff(camera_time) > 0), "Camera timestamps must be strictly increasing." f" Recording: {recording_directory}"

        # Open up the saved processed world camera information
        with h5py.File(recording_directory / "W" / "transformed.hdf5", "r") as saved:
            # Get the frames dataset and the saved FPS
            frames: h5py.Dataset = saved["data/frames"]
            camera_fps: float = float(saved["metadata/fps"][()])
            assert np.isfinite(camera_fps), f"Saved camera FPS must be finite; got {camera_fps}." f" Recording: {recording_directory}"
            assert camera_fps > 0, f"Saved camera FPS must be positive; got {camera_fps}." f" Recording: {recording_directory}"
            assert len(frames) == len(camera_time), "Camera timestamps and saved frames must have equal lengths." f" Recording: {recording_directory}"

            # Include the final stored frame for its nominal duration by default.
            start_time: float
            end_time: float

            # Define the range of the video to illustrate. The above plot already
            # took care of this, we just also need to do this for the video part
            # If if no specific time range, just take the first adn last timestamp
            if time_range is None:
                start_time = float(camera_time[0])

                # The end is exlusive, so need to add 1/camera_fps at the end
                end_time = float(camera_time[-1] + 1 / camera_fps)

            # Otherwise, use the specified time range
            else:
                start_time = float(time_range[0])
                end_time = float(time_range[1])

            # Each timestamp below becomes one output video frame. For example,
            # gif_fps=5 samples the recording every 0.2 seconds. Encoding those
            # images at 5 frames per second keeps playback at the recorded speed.
            # Stop before end_time to match the plotter's [start, end) selection.
            animation_times: np.ndarray = np.arange(start_time, end_time, 1 / gif_fps)

            # Add one red vertical line to each visible time-series panel.
            # Keep the line objects so we can move them instead of rebuilding plots.
            # The returned axes array excludes the colorbar.
            cursors: list[Line2D] = [
                axis.axvline(start_time, color="red", linewidth=1.5)
                for axis in axes.flat if axis.get_visible()
            ]

            # Render the figure before reading its pixel dimensions. These dimensions
            # determine the width of the video and the space reserved for the plot.
            canvas.draw()

            plot_width: int
            plot_height: int
            plot_width, plot_height = canvas.get_width_height()

            # Reserve a camera area above the plot, half as tall as the plot is wide.
            # The camera image will be resized and centered within this area.
            camera_height: int = plot_width // 2

            # Stack the camera area above the plot to get the total video height.
            # Adding dimension % 2 adds one pixel only when that dimension is odd.
            # Even dimensions avoid encoder cropping, and all frames use this same size.
            video_width: int = plot_width + plot_width % 2
            video_height: int = plot_height + camera_height
            video_height += video_height % 2

            # The display converter uses neighboring source frames to calculate
            # the mean and standard deviation used to scale L+M into display values.
            # Use one second at camera_fps, not gif_fps: these are camera frames,
            # not output video frames. A shorter recording uses all available frames.
            # The actual window position is chosen separately for each frame below.
            window_frames: int = min(len(frames), max(1, int(round(camera_fps))))
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Save rendered images on disk so we do not keep the entire video in RAM.
            # This directory and its PNGs are removed automatically when the block exits.
            with TemporaryDirectory(prefix="processed_animation_") as temporary_directory:
                image_directory: Path = Path(temporary_directory) / "frames"
                image_directory.mkdir()

                frame_number: int
                timestamp: np.floating
                cursor: Line2D

                # frame_number names the output PNG; timestamp selects the camera
                # image and positions every red time bar for that same instant.
                for frame_number, timestamp in enumerate(tqdm(
                    animation_times, desc="Rendering processed animation", disable=not verbose,
                )):
                    # Move both endpoints of each vertical line to the current time.
                    # The sensor curves still describe the full selected
                    # interval; only the red markers change from one frame to the next.
                    for cursor in cursors:
                        cursor.set_xdata([timestamp, timestamp])

                    canvas.draw()

                    # Capture the newly rendered plot as an image. Copy the pixel
                    # buffer because Matplotlib reuses it on the next canvas.draw().
                    # Convert RGBA to RGB so it matches the camera image and video.
                    plot_image: Image.Image = Image.fromarray(
                        np.asarray(canvas.buffer_rgba()).copy(),
                    ).convert("RGB")

                    # searchsorted(..., side="right") finds the first camera time
                    # strictly after this output time. Subtracting one selects the
                    # most recent camera frame at or before it, including exact matches.
                    # A result of -1 means the output time precedes the first camera frame.
                    source_index: int = int(np.searchsorted(camera_time, timestamp, side="right") - 1)
                    camera_image: Image.Image = Image.new("RGB", (plot_width, camera_height), "black")

                    # Start with a black camera area and replace it only if a recent
                    # source frame exists. One nominal camera frame lasts 1/camera_fps
                    # seconds; the 1.5 multiplier allows modest timestamp jitter.
                    # If the latest frame is older than that, leave the area black
                    # instead of holding an old image across a missing-frame gap.
                    if source_index >= 0 and timestamp - camera_time[source_index] < 1.5 / camera_fps:
                        # Center the normalization window on the selected source frame.
                        # max(..., 0) prevents starting before the recording; the outer
                        # min keeps the full window inside the end of the dataset.
                        # Near either end, shift the window inward without shortening it.
                        # This window may include frames outside the requested time range
                        # because those neighbors are used for display scaling only.
                        window_start: int = min(
                            max(source_index - window_frames // 2, 0),
                            len(frames) - window_frames,
                        )
                        opponent_frames: np.ndarray = frames[window_start:window_start + window_frames]
                        display_frames: np.ndarray = world_util.opponentmap_to_videoframe(
                            opponent_frames, fps=camera_fps, window_seconds=1.0,
                        )

                        # source_index refers to the full recording; subtract window_start
                        # to locate that frame within the small buffer we just converted.
                        # Select channel zero (L+M) and represent its grayscale values
                        # as RGB. The converter already handles the image orientation.
                        local_index: int = source_index - window_start
                        lm_image: Image.Image = Image.fromarray(
                            display_frames[local_index, ..., 0],
                        ).convert("RGB")

                        # Calculate the scale needed to fit each dimension, then choose
                        # the smaller scale so neither width nor height exceeds the area.
                        # Apply that same scale to both dimensions to preserve proportions.
                        # LANCZOS smooths the resized image; max(1, ...) avoids zero pixels.
                        image_scale: float = min(plot_width / lm_image.width, camera_height / lm_image.height)
                        lm_image = lm_image.resize(
                            (max(1, int(lm_image.width * image_scale)),
                             max(1, int(lm_image.height * image_scale))),
                            Image.Resampling.LANCZOS,
                        )

                        # Split the unused width and height equally on both sides.
                        # These offsets center the image within the black camera area.
                        camera_position: tuple[int, int] = (
                            (plot_width - lm_image.width) // 2,
                            (camera_height - lm_image.height) // 2,
                        )
                        camera_image.paste(lm_image, camera_position)

                    # Compose the final frame: camera at the top, plot directly below.
                    # PIL positions are (left, top), so the plot starts at camera_height.
                    # The white background also fills any one-pixel padding added above.
                    combined_image: Image.Image = Image.new(
                        "RGB", (video_width, video_height), "white",
                    )
                    combined_image.paste(camera_image, (0, 0))
                    combined_image.paste(plot_image, (0, camera_height))

                    # Save one PNG per output timestamp. Eight-digit numbering gives
                    # each frame a unique name in playback order (00000000.png, etc.).
                    # dir_to_video reads these images after the rendering loop finishes.
                    combined_image.save(image_directory / f"{frame_number:08d}.png")

                # Let's turn this directory of frames into a video, output at a temporary location for now
                temporary_video: Path = Path(temporary_directory) / "animation.avi"
                video_io.dir_to_video(str(image_directory), str(temporary_video), fps=gif_fps)

                # Move the finished video to the requested output location.
                shutil.move(str(temporary_video), str(output_path))

        return output_path

    finally:
        # Release the plot even if reading, rendering, or video encoding fails.
        plt.close(figure)


def main() -> None:
    """Provide a placeholder script entry point without processing data.

    Args:
        None.

    Returns:
        None.
    """
    pass


if __name__ == "__main__":
    main()
