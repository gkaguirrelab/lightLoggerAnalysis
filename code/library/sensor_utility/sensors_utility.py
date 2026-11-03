"""Utilities that operate on complete light-logger recordings."""

import os
import sys
from typing import Literal
from collections.abc import Iterable

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from natsort import natsorted
from tqdm.auto import tqdm

sys.path.append(os.path.dirname(__file__))
import ms_util
import world_util

# Raw world timestamps use nanoseconds; MS timestamps use seconds.
WORLD_NANOSECONDS_TO_SECONDS: float = 1e-9


def group_sensors_files(recording_path: str) -> dict[str, list[tuple[str, str]]]:
    """Return naturally sorted raw metadata/data path pairs for each sensor.

    Args:
        recording_path: Directory to scan for world, pupil, and ms filenames.

    Returns:
        Dictionary with W, P, and M keys, each containing a list of
        (metadata_path, data_path) pairs. Missing sensors have empty lists.
        Data names are derived by removing "_metadata" from metadata names;
        this function does not check that the corresponding data files exist.

    Raises:
        OSError: If the recording directory cannot be listed."""

    def group_sensor_files(sensor_name: str) -> list[tuple[str, str]]:
        """Find metadata files for one sensor in the enclosing recording directory.

        Args:
            sensor_name: Filename substring to match, such as "world" or "ms".

        Returns:
            Naturally sorted (metadata_path, data_path) pairs. Matching names
            must contain both sensor_name and "metadata"."""
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
        files: Ordered (metadata .npy path, frame-buffer .npy path) pairs.
        output_path: HDF5 filename, not a directory. Parent directories are created.
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

    Raises:
        AssertionError: World metadata has neither 4 legacy nor 6 modern columns.
        ValueError: An existing file is incomplete, frame counts disagree, or
            no nonempty chunks are selected.
        OSError: Input arrays cannot be read or the HDF5 file cannot be written.

    Notes:
        The root completed attribute becomes True only after all selected chunks
        are written. Exceptions close the file and pool; partial files remain
        marked incomplete. Each chunk is flushed before processing the next.
        Script callers using worker processes need an if __name__ == "__main__"
        guard. Worker benchmarks are recorded in tmp/benchmark_imputation_workers.json
        and tmp/benchmark_demosaic_large_buffers.json."""

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
        return

    # Validate gaze before opening (and potentially replacing) the output file.
    gaze = None if gaze_angles_xy is None else np.asarray(gaze_angles_xy, dtype=float)
    if gaze is not None and gaze.shape != (2,):
        total_frames = sum(len(np.load(path, mmap_mode="r")) for _, path in files)
        if gaze.shape != (total_frames, 2):
            raise ValueError("Gaze must have one X/Y pixel pair per input world frame")

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
        gaze_offset = 0
        
        # Iterate over the paired metadata and frame-buffer files.
        for file_num, (metadata_path, data_path) in file_iterator:
            frame_count = len(np.load(data_path, mmap_mode="r"))
            chunk_gaze = (gaze if gaze is None or gaze.shape == (2,)
                          else gaze[gaze_offset:gaze_offset + frame_count])
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

            # Run the transformation pipeline
            data_and_processing_metadata: dict[str, object] = world_util.world_transformation_pipeline(
                frame_buffer, agc_settings, chunk_gaze, age, n_workers=n_workers, demosaic_pool=demosaic_pool,
            )
            transformed_data: np.ndarray = data_and_processing_metadata["data"]

            # Save calibration/settings from the first processed chunk only.
            # Exclude AGC here because its per-frame values are appended separately.
            if processing_metadata is None:
                processing_metadata = {
                    key: value for key, value in data_and_processing_metadata["metadata"].items()
                    if key != "agc_settings"
                }
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

def _mean_camera_rgb_by_ms_window(
    timestamps: np.ndarray,
    world_files: list[tuple[str, str]],
    world_output_path: str,
    required: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average camera RGB radiance within minispectrometer midpoint windows.

    Args:
        timestamps: Strictly increasing finite MS timestamps in seconds, shape (N,).
        world_files: Ordered raw world (metadata_path, data_path) pairs. Only the
            metadata paths are read here, to obtain frame timestamps in nanoseconds.
        world_output_path: World output directory containing transformed.hdf5,
            or legacy world_chunk<index>.mat files when the HDF5 file is absent.
        required: Boolean array of shape (N,); only True windows are accumulated.

    Returns:
        Tuple of mean RGB radiance (N, 3), frame counts (N,), and window bounds
        in seconds (N, 2). Midpoints belong to the earlier window (left-open,
        right-closed). Outer bounds are infinite. Empty or unrequested windows
        have NaN means and zero counts. No temporal offset is applied.

    Raises:
        ValueError: Timestamps are invalid, saved timestamps do not match raw
            metadata, an HDF5 file is incomplete, or required RGB values are invalid.
        FileNotFoundError: A required legacy world chunk is missing.

    Notes:
        Supports /data/frames and earlier /data datasets in native HDF5 files.
        Legacy MATLAB files can supply precomputed mean_camera_rgb or RGB frames.
        Native HDF5 /data/average_frames is preferred to avoid reading images.
        If absent, frames are read individually so the full video need not fit in memory."""
    if not np.all(np.isfinite(timestamps)) or np.any(np.diff(timestamps) <= 0):
        raise ValueError("Minispect timestamps must be finite and strictly increasing.")
    if not len(timestamps):
        return np.empty((0, 3)), np.empty(0, dtype=np.int64), np.empty((0, 2))

    # Adjacent MS timestamps define midpoint boundaries for the frame bins.
    boundaries: np.ndarray = timestamps[:-1] + np.diff(timestamps) / 2
    edges: np.ndarray = np.r_[-np.inf, boundaries, np.inf]
    sums: np.ndarray = np.zeros((len(timestamps), 3))
    counts: np.ndarray = np.zeros(len(timestamps), dtype=np.int64)

    for index, (metadata_path, _) in enumerate(world_files):
        world_t = np.load(metadata_path)[:, 0]
        if not np.all(np.isfinite(world_t)) or np.any(np.diff(world_t) <= 0):
            raise ValueError("World timestamps must be finite and strictly increasing.")
        bins = np.searchsorted(boundaries, world_t / 1e9, side="left")
        wanted = required[bins]
        if not np.any(wanted):
            continue
        hdf5_path: str = os.path.join(world_output_path, "transformed.hdf5")
        if os.path.isfile(hdf5_path):
            with h5py.File(hdf5_path, "r") as saved:
                if not saved.attrs.get("completed", False):
                    raise ValueError(f"Incomplete world output: {hdf5_path}")
                saved_t: np.ndarray = saved["agcSettings/timestamps"][:]
                # Locate raw timestamps so selected chunk ranges remain supported.
                positions: np.ndarray = np.searchsorted(saved_t, world_t)
                if np.any(positions >= len(saved_t)) or not np.array_equal(saved_t[positions], world_t):
                    raise ValueError(f"World output timestamps do not match raw metadata: {hdf5_path}")
                rgb = np.zeros((len(world_t), 3))
                # Prefer the current path, retaining support for earlier cached files.
                averages_path: str = "data/average_frames" if "data/average_frames" in saved else "averages/mean_camera_rgb"
                if averages_path in saved:
                    # Read only the small precomputed RGB rows needed by MS windows.
                    averages: h5py.Dataset = saved[averages_path]
                    if averages.shape != (len(saved_t), 3):
                        raise ValueError(f"Camera averages/timestamps shape mismatch: {hdf5_path}")
                    rgb[wanted] = averages[positions[wanted]]
                else:
                    # Older files remain usable, but need spatial means recomputed.
                    frames: h5py.Dataset = saved["data/frames"] if isinstance(saved["data"], h5py.Group) else saved["data"]
                    for frame in np.flatnonzero(wanted):
                        rgb[frame] = frames[positions[frame]].mean(axis=(0, 1))
            path = hdf5_path
        else:
            path = os.path.join(world_output_path, f"world_chunk{index}.mat")
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"Joint spectra require the overlapping world chunk: {path}"
                )

            with h5py.File(path, "r") as saved:
                saved_t = np.asarray(saved["metadata/timestamps"]).ravel()
                if not np.array_equal(world_t, saved_t):
                    raise ValueError(
                        f"World output timestamps do not match raw metadata: {path}"
                    )
                if "mean_camera_rgb" in saved["metadata"]:
                    rgb = np.asarray(saved["metadata/mean_camera_rgb"]).T
                else:
                    data = saved["data"]
                    if (
                        data.ndim != 4
                        or data.shape[0] != 3
                        or data.shape[-1] != len(world_t)
                    ):
                        raise ValueError(f"Reprocess legacy non-RGB world output: {path}")
                    rgb = np.zeros((len(world_t), 3))
                    for frame in np.flatnonzero(wanted):
                        rgb[frame] = np.asarray(data[..., frame]).mean(axis=(1, 2))
        if rgb.shape != (len(world_t), 3) or not np.all(np.isfinite(rgb[wanted])):
            raise ValueError(f"Invalid camera RGB means in {path}")

        # Accumulate across chunk boundaries without retaining camera images.
        np.add.at(sums, bins[wanted], rgb[wanted])
        np.add.at(counts, bins[wanted], 1)

    # Empty windows remain NaN so they cannot be mistaken for dark scenes.
    means: np.ndarray = np.full_like(sums, np.nan)
    np.divide(sums, counts[:, None], out=means, where=counts[:, None] > 0)
    return means, counts, np.column_stack((edges[:-1], edges[1:]))


def _process_raw_ms_helper(
    ms_files: list[tuple[str, str]],
    world_files: list[tuple[str, str]],
    output_path: str,
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
        ms_files: Ordered (timestamp .npy path, reading-buffer .npy path) pairs.
            Each timestamp array flattens to (N,) in seconds on the world
            recording clock. Each reading buffer has shape (N, channels), with
            the first seven columns containing AS7341 F1-F7 measurements.
        output_path: MS destination directory, typically GKA/<attempt>/M.
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
        world_files: Raw world metadata/data pairs required for camera alignment.
            Requires the sibling W/transformed.hdf5 with /data/average_frames.

    Returns:
        None. Writes transformed.hdf5 with these /data fields (NumPy shapes):
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
        Every MS reading must have camera coverage and finite averaged RGB values.
        Empty buffers are skipped; empty averaging windows raise AssertionError.

    Raises:
        AssertionError: Any MS reading has no world frames or nonfinite averaged RGB.
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
        Per-reading arrays use extendible chunks, gzip level 4, and shuffle.
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
                fitted_radiance, sampling, fitted_value, fitted_errors = ms_util.estimate_joint_radiance_spectrum(
                    reading_buffer[has_world_frames], rgb[has_world_frames],
                    base_alpha=base_alpha, reference_intensity=reference_intensity, beta=beta,
                )
                radiance[has_world_frames] = fitted_radiance
                fit_value[has_world_frames] = fitted_value
                fit_errors[has_world_frames] = fitted_errors

            # 5. Keep fixed run settings separate from arrays that grow with each chunk.
            # These settings describe the actual fit above; calibration coefficients
            # come directly from ms_util, rather than being copied by hand here.
            processing_settings: dict[str, object] = {
                "world_frame_index_bounds_convention": "zero-based world row indices, start inclusive, stop exclusive",
                "estimator": "estimateJointRadianceSpectrum",
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
    path_to_raw: str,
    output_path: str,
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
        path_to_raw: Existing nonempty directory of raw sensor chunk files.
        output_path: Parent output directory. World data go to W/transformed.hdf5;
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
        None. Writes native HDF5 world and MS output files.
        Pupil files may be discovered but are not processed by this function.

    Raises:
        AssertionError: The input directory is missing/empty or world metadata
            has an unsupported column count.
        ValueError: World output is incomplete, selected world chunks are empty,
            or data validation/alignment fails in either processing helper.
        OSError: Input or output files cannot be accessed.

    Notes:
        Selected world chunks must cover the requested MS windows. MS uses saved
        per-frame RGB averages when available. Script callers using worker
        processes need an if __name__ == "__main__" guard."""
    # Ensure the raw recording is an existing, non-empty directory.
    assert (
        os.path.isdir(path_to_raw) and len(os.listdir(path_to_raw)) > 0
    ), f"{path_to_raw} must be an existing, non-empty directory"

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
        files_by_sensor["W"], os.path.join(world_output_path, "transformed.hdf5"),
        age=age, gaze_angles_xy=gaze_angles_xy, chunk_range=chunk_ranges["W"],
        overwrite_existing=overwrite_existing, verbose=verbose,
        n_workers=n_workers, demosaic_n_workers=demosaic_n_workers,
    )

    # This still runs when world processing skips an already completed file.
    _process_raw_ms_helper(
        ms_files=files_by_sensor["M"], world_files=files_by_sensor["W"],
        output_path=os.path.join(output_path, "M"),
        chunk_range=chunk_ranges["M"], overwrite_existing=overwrite_existing,
        verbose=verbose,
    )


def _plot_raw_world_helper(
    path_to_raw: str,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Plot raw world channel means and settings from a chunk directory.

    Args:
        path_to_raw: Directory containing world data and metadata chunks.
        axes: Two axes for channel means and settings, or None for a local figure.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and two primary axes; time is in seconds with frame numbers above.
    """
    if axes is None:
        figure, axes = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
    else:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if axes.size != 2:
            raise ValueError("Raw world plotting requires two axes.")
        figure = axes[0].figure

    # Stored rows stay aligned with the raw frames, including acquisition gaps.
    metadata = world_util.world_metadata_from_chunks(
        path_to_raw, fill_missing_frames=False, verbose=verbose,
    )
    world_util.plot_world_camera_channels(path_to_raw, ax=axes[0], verbose=verbose)
    _, (gain_axis, exposure_axis) = world_util.plot_world_camera_settings(
        metadata, ax=axes[1], time_origin=0,
    )
    gain_axis.set_ylim(0, 11)
    exposure_axis.set_ylim(0, 8333)
    gain_axis.set_xlabel("time [s]")
    lines = gain_axis.get_lines() + exposure_axis.get_lines()
    gain_axis.legend(lines, [line.get_label() for line in lines], loc="upper right")
    return figure, axes


def _plot_raw_ms_helper(
    path_to_raw: str,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Collect raw MS chunks and plot AS, TS, LS, and temperature readings.

    Args:
        path_to_raw: Directory containing paired MS metadata and data chunks.
        axes: A (2, 2) axes array, or None to create a local figure.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and (2, 2) axes array populated by ``ms_util.plot_chunk_readings``.

    Raises:
        ValueError: Chunks are empty, data and timestamps disagree, or axes are invalid.
    """
    timestamp_buffer: list[np.ndarray] = []
    reading_buffer: list[np.ndarray] = []

    # Use the MS timestamps directly; they are already in seconds.
    for metadata_path, data_path in tqdm(group_sensors_files(path_to_raw)["M"], desc="Loading raw MS chunks", disable=not verbose):
        timestamps = np.load(metadata_path).reshape(-1)
        readings = np.load(data_path)
        if readings.ndim != 2 or len(readings) != len(timestamps):
            raise ValueError(f"MS readings and timestamps disagree: {data_path}")
        if len(readings):
            timestamp_buffer.append(timestamps)
            reading_buffer.append(readings)

    if not reading_buffer:
        raise ValueError(f"No nonempty MS chunks found: {path_to_raw}")
    timestamps = np.concatenate(timestamp_buffer)
    readings = np.concatenate(reading_buffer)

    if axes is None:
        figure, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    else:
        axes = np.asarray(axes, dtype=object)
        if axes.shape != (2, 2):
            raise ValueError("Raw MS plotting requires a (2, 2) axes array.")
        figure = axes.flat[0].figure

    ms_util.plot_chunk_readings(timestamps, readings, axes=axes, verbose=verbose)
    return figure, axes


def plot_raw_recording(path_to_raw: str, verbose: bool = False) -> tuple[Figure, np.ndarray]:
    """Plot raw world and MS chunks together in a two-row, three-column figure.

    Args:
        path_to_raw: Directory directly containing world and MS .npy chunks.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and (2, 3) primary axes. World means/settings occupy the left
        column; four MS chip panels occupy the remaining columns. Times use
        each sensor's own time vector in seconds. Call ``plt.show()`` to display.
    """
    path_to_raw = os.path.abspath(os.path.expanduser(path_to_raw))
    figure, axes = plt.subplots(2, 3, figsize=(20, 10), constrained_layout=True)

    # Each helper draws on its assigned panels instead of creating a figure.
    try:
        _plot_raw_world_helper(path_to_raw, axes[:, 0], verbose=verbose)
        _plot_raw_ms_helper(path_to_raw, axes[:, 1:], verbose=verbose)
    except Exception:
        plt.close(figure)
        raise

    figure.suptitle(f"Raw recording: {os.path.basename(path_to_raw)}")
    return figure, axes


def _plot_processed_world_helper(
    world_directory: str,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Plot RGB means, opponent means, and AGC settings for processed world data.

    Args:
        world_directory: W directory containing ``transformed.hdf5``.
        axes: Three Matplotlib axes in plot order, or None to create a
            local three-row figure. The AGC panel adds a right exposure axis.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and the three primary axes. The figure is not shown or saved.
        Times use the world recording's own timestamps converted to seconds.
        RGB and opponent values are calibrated radiance and use automatic limits.

    Raises:
        ValueError: Axes count or frame/metadata shapes are incompatible.
        OSError: The HDF5 file cannot be opened.
        KeyError: A required dataset is missing.
    """
    # Read scalar metadata and RGB means; reduce opponent images in small batches.
    with h5py.File(os.path.join(world_directory, "transformed.hdf5"), "r") as saved:
        time: np.ndarray = saved["agcSettings/timestamps"][:] * WORLD_NANOSECONDS_TO_SECONDS
        rgb: np.ndarray = saved["data/average_frames"][:]
        frames: h5py.Dataset = saved["data/frames"]
        settings: np.ndarray = np.column_stack([
            saved[f"agcSettings/{name}"][:] for name in ("Again", "Dgain", "exposure")
        ])

        # Each saved row must describe the same world frame.
        if (
            time.ndim != 1
            or rgb.shape != (len(time), 3)
            or frames.ndim != 4
            or frames.shape[0] != len(time)
            or frames.shape[-1] != 3
            or settings.shape != (len(time), 3)
        ):
            raise ValueError("World frames, RGB means, settings, and timestamps must align.")

        # Count the batches first, including a final partially filled chunk.
        chunk_size: int = 1000
        num_chunks: int = (len(time) + chunk_size - 1) // chunk_size
        opponent: np.ndarray = np.full((len(time), 3), np.nan)

        for chunk_index in tqdm(
            range(num_chunks), desc="Loading world chunks", disable=not verbose,
        ):
            start: int = chunk_index * chunk_size
            stop: int = min(start + chunk_size, len(time))
            buffer: np.ndarray = frames[start:stop]

            # Reduce each frame separately so the inner bar reports frame progress.
            # Frames without finite coverage retain their initialized NaN values.
            for frame_index in tqdm(
                range(len(buffer)), desc="Averaging opponent frames",
                disable=not verbose, leave=False,
            ):
                frame: np.ndarray = buffer[frame_index]
                finite: np.ndarray = np.isfinite(frame)
                counts: np.ndarray = finite.sum(axis=(0, 1))
                totals: np.ndarray = np.where(finite, frame, 0).sum(axis=(0, 1))
                np.divide(
                    totals, counts, out=opponent[start + frame_index], where=counts > 0,
                )

    # Create a local figure only when the caller has not supplied axes.
    if axes is None:
        figure, axes = plt.subplots(
            3, 1, figsize=(12, 10), sharex=True, constrained_layout=True
        )
    else:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if len(axes) != 3:
            raise ValueError("World plotting requires three axes.")
        figure = axes[0].figure

    # These are radiance values, not pixel counts, so retain automatic y limits.
    for channel, (label, color) in enumerate(tqdm(
        zip(("R", "G", "B"), ("red", "green", "blue"))
    , desc="Plotting RGB channels", disable=not verbose)):
        axes[0].plot(time, rgb[:, channel], color=color, label=label)
    axes[0].set_title("World: average RGB radiance")
    axes[0].set_ylabel("Average Channel Value")

    # Plot the spatial mean of each opponent channel on its own trace.
    for channel, (label, color) in enumerate(tqdm(
        zip(("LM", "L-M", "S"), ("black", "darkorange", "purple"))
    , desc="Plotting opponent channels", disable=not verbose)):
        axes[1].plot(time, opponent[:, channel], color=color, label=label)
    axes[1].set_title("World: average opponent radiance")
    axes[1].set_ylabel("Average Channel Value")

    # Add stored-frame coordinates above both radiance panels.
    for axis in tqdm(axes[:2], desc="Formatting world panels", disable=not verbose):
        axis.set_xlabel("time [s]")
        axis.grid(True, alpha=0.25)
        axis.legend(loc="upper right")
        frame_timestamps = np.asarray(time, dtype=float)
        if (frame_timestamps.ndim != 1 or frame_timestamps.size == 0
                or not np.all(np.isfinite(frame_timestamps)) or np.any(np.diff(frame_timestamps) <= 0)):
            raise ValueError("Frame frame_timestamps must be finite, nonempty, and increasing.")

        # Interpolate frame positions so dropped frames do not distort the time scale.
        frame_indices: np.ndarray = np.arange(len(frame_timestamps), dtype=float)
        if len(frame_timestamps) == 1:
            frame_timestamps = np.array([frame_timestamps[0], frame_timestamps[0] + 1])
            frame_indices = np.array([0., 1.])

        def interpolate(values: np.ndarray, source: np.ndarray, target: np.ndarray) -> np.ndarray:
            """Interpolate positions, extending the end slopes into plot margins."""
            result = np.interp(values, source, target)
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

        frame_axis = axis.secondary_xaxis(
            "top", functions=(
                lambda times: interpolate(times, frame_timestamps, frame_indices),
                lambda frames: interpolate(frames, frame_indices, frame_timestamps),
            ),
        )
        frame_axis.set_xlabel("Frame number")
        frame_axis.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))

    # Use the same metadata reader and settings plot as raw recordings.
    metadata = world_util.world_metadata_from_chunks(
        os.path.join(world_directory, "transformed.hdf5"), fill_missing_frames=False, verbose=verbose,
    )
    _, (gain_axis, exposure_axis) = world_util.plot_world_camera_settings(
        metadata, ax=axes[2], time_origin=0, verbose=verbose,
    )
    gain_axis.set_ylim(0, 11)
    exposure_axis.set_ylim(0, 8333)
    gain_axis.set_xlabel("time [s]")
    lines = gain_axis.get_lines() + exposure_axis.get_lines()
    gain_axis.legend(lines, [line.get_label() for line in lines], loc="upper right")

    return figure, axes


def _plot_processed_ms_helper(
    ms_directory: str,
    axes: np.ndarray | None = None,
    verbose: bool = False,
) -> tuple[Figure, np.ndarray]:
    """Plot spectral radiance, fit value, and signed fit errors over MS time.

    Args:
        ms_directory: M directory containing ``transformed.hdf5``.
        axes: Three axes in plot order, or None to create a local three-row
            figure. The radiance panel adds a colorbar.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and three primary axes, without showing or saving the figure.
        Radiance is a time-by-wavelength heatmap; fit errors are individual
        F1-F7 and camera RGB residual traces. Timestamps retain their saved
        seconds without adjusting them to the world time vector.

    Raises:
        ValueError: Axes count, timestamps, or dataset shapes are incompatible.
        OSError: The HDF5 file cannot be opened.
        KeyError: A required dataset is missing.
    """
    # Load the spectral estimates and their per-reading fit diagnostics.
    with h5py.File(os.path.join(ms_directory, "transformed.hdf5"), "r") as saved:
        time: np.ndarray = saved["timestamps"][:]
        radiance: np.ndarray = saved["data/radiance"][:]
        sampling: np.ndarray = saved["data/S"][:].reshape(-1)
        fit_value: np.ndarray = saved["data/fVal"][:]
        fit_errors: np.ndarray = saved["data/fitErrors"][:]

    # The heatmap requires ordered times and one spectrum per timestamp.
    if (
        time.ndim != 1
        or len(time) == 0
        or not np.all(np.isfinite(time))
        or np.any(np.diff(time) <= 0)
    ):
        raise ValueError("MS timestamps must be nonempty, finite, and increasing.")

    if (
        sampling.size != 3
        or sampling[1] <= 0
        or radiance.shape != (len(time), int(sampling[2]))
        or fit_value.shape != (len(time),)
        or fit_errors.shape != (len(time), ms_util.MS_JOINT_CHANNEL_COUNT + 3)
    ):
        raise ValueError("MS radiance, wavelength sampling, fit values, and errors must align.")

    # S stores the starting wavelength, step size, and sample count.
    wavelengths: np.ndarray = sampling[0] + sampling[1] * np.arange(int(sampling[2]))

    # Create a local figure only when the caller has not supplied axes.
    if axes is None:
        figure, axes = plt.subplots(
            3, 1, figsize=(12, 10), sharex=True, constrained_layout=True
        )
    else:
        axes = np.asarray(axes, dtype=object).reshape(-1)
        if len(axes) != 3:
            raise ValueError("MS plotting requires three axes.")
        figure = axes[0].figure

    # Keep every wavelength visible instead of collapsing the spectrum to a mean.
    mesh = axes[0].pcolormesh(
        time, wavelengths, np.ma.masked_invalid(radiance.T),
        shading="auto", cmap="viridis",
    )
    figure.colorbar(mesh, ax=axes[0], label="Spectral radiance")
    axes[0].set_title("MS: spectral radiance over time")
    axes[0].set_ylabel("Wavelength [nm]")

    # Show the overall fit value separately from the signed channel errors.
    axes[1].plot(time, fit_value, color="tab:blue")
    axes[1].set_title("MS: fit value")
    axes[1].set_ylabel("Fit Value")

    # The saved residual columns are F1-F7 followed by camera R, G, and B.
    labels: list[str] = [
        f"F{index + 1}" for index in range(ms_util.MS_JOINT_CHANNEL_COUNT)
    ] + ["R", "G", "B"]
    for channel, label in enumerate(tqdm(labels, desc="Plotting fit errors", disable=not verbose)):
        axes[2].plot(
            time, fit_errors[:, channel],
            color=plt.get_cmap("tab10")(channel), label=label,
        )
    axes[2].set_title("MS: signed fit errors")
    axes[2].set_ylabel("Fit Error")
    axes[2].legend(loc="upper right", ncol=2)

    # Apply common time labels and panel formatting.
    for axis in tqdm(axes, desc="Labeling MS panels", disable=not verbose):
        axis.set_xlabel("time [s]")
    for axis in tqdm(axes[1:], desc="Formatting MS panels", disable=not verbose):
        axis.grid(True, alpha=0.25)

    return figure, axes


def plot_processed_recording(recording_directory: str, verbose: bool = False) -> tuple[Figure, np.ndarray]:
    """Plot processed world and MS data together in a 3-by-2 figure.

    Args:
        recording_directory: Parent directory containing W/transformed.hdf5
            and M/transformed.hdf5.

        verbose: Show tqdm progress bars while loading and plotting data.

    Returns:
        Figure and a (3, 2) array of primary axes. The left column contains
        world RGB means, opponent means, and AGC settings. The right contains
        MS radiance, fit values, and fit errors. Added exposure/colorbar axes
        are available through the figure. Call ``plt.show()`` to display or
        ``figure.savefig(...)`` to save the result.

    Raises:
        OSError: Either processed sensor file cannot be opened.
        KeyError: A required dataset is missing.
        ValueError: Saved dataset shapes or timestamps are incompatible.
    """
    # Place world and MS panels side by side with independent time axes.
    recording_directory = os.path.abspath(os.path.expanduser(recording_directory))
    figure, axes = plt.subplots(
        3, 2, figsize=(18, 12), constrained_layout=True
    )

    # Each helper plots its sensor's own time vector on the supplied axes.
    try:
        _plot_processed_world_helper(os.path.join(recording_directory, "W"), axes[:, 0], verbose=verbose)
        _plot_processed_ms_helper(os.path.join(recording_directory, "M"), axes[:, 1], verbose=verbose)
    except Exception:
        # Avoid leaving a partially populated figure open if loading fails.
        plt.close(figure)
        raise

    figure.suptitle(f"Processed recording: {os.path.basename(recording_directory)}")

    return figure, axes


def main() -> None:
    """Placeholder script entry point; performs no processing and returns None."""
    pass


if __name__ == "__main__":
    main()
