"""Chunk-based I/O for light logger sensor recordings.

Provides parsers for world camera, pupil camera, and minispect/sunglasses
sensor data stored as chunked numpy files, along with utilities for grouping
chunk files by sensor and locating chunks by timestamp.
"""

import gc
import os
from natsort import natsorted
import numpy as np
import pathlib
import sys
from tqdm.auto import tqdm

# Import the sensor utility libraries 
light_logger_analysis_path: str = str(pathlib.Path(__file__).parents[3])
sensor_utility_path: str = os.path.join(light_logger_analysis_path, "library", "sensor_utility")
for path in (sensor_utility_path,):
    sys.path.append(path)

import world_util 
import ms_util 
import pupil_util

def minispect_and_sunglasses_chunk_parser(chunk_paths: tuple[str],
                                          apply_digital_gain: bool=False, 
                                          convert_time_units: bool=False,
                                          use_mean_frame: bool=False, 
                                          convert_to_float: bool=False, 
                                          apply_phase_offset: bool=False, 
                                          apply_RGB_correction: bool=False, 
                                          apply_fielding_function: bool=False,
                                          debayer_images: bool=False,
                                          mean_axes: tuple[int] = (1, 2),
                                          contains_AGC_metadata: bool=False, 
                                          password="1234"
                                         ) -> dict[str, dict]:
    """Parse one mixed minispect/sunglasses chunk into project-standard data.

    A single MS chunk stores multiple logical sensors together. This parser
    loads the shared timestamp vector and packed value matrix, derives the
    accelerometer timebase used by the minispect utilities, splits the
    trailing sunglasses samples away from the multisensor payload, and then
    delegates the minispect decoding to ``ms_util.parse_readings``.

    Args:
        chunk_paths: Tuple of (metadata_path, values_path) for the chunk,
            or an empty tuple if no readings exist.
        apply_digital_gain: Unused placeholder so this parser shares a
            common call signature with the camera parsers.
        convert_time_units: Unused because these timestamps are already in
            seconds.
        use_mean_frame: Unused because these values are not image frames.
        convert_to_float: Whether to convert values to float64.
        apply_phase_offset: Whether to apply the MS time offset correction.
        apply_RGB_correction: Unused for this non-image sensor family.
        apply_fielding_function: Unused for this non-image sensor family.
        debayer_images: Unused for this non-image sensor family.
        mean_axes: Unused placeholder for signature compatibility.
        contains_AGC_metadata: Unused because AGC metadata is not expected
            in MS chunks.
        password: Reserved for future encrypted-chunk support; currently
            unused.

    Returns:
        Dictionary with 'M' (minispect) and 'S' (sunglasses) keys, each
        containing 't' (timestamps) and 'v' (values) sub-dictionaries.
    """

    # Initialize return variable
    chunk_dict: dict[str, dict] = {'M': {}, 'S': {}}

    # If there were no readings for this sensor 
    if(len(chunk_paths) == 0):
        # Initialize the M field as empty 
        chunk_dict['M']['t'] = {sensor_name: np.array([], dtype=np.float64)
                                for sensor_name in ms_util.MS_SENSOR_NAMES
                               } 

        chunk_dict['M']['v'] = {sensor_name: np.empty((0, num_channels), dtype=np.float64) 
                                for sensor_name, num_channels in ms_util.MS_NAMES_AND_CHANNELS.items() 
                                }
        
        # Initialize the sunglasses field as empty
        chunk_dict['S']['t'] = np.array([], dtype=np.float64)
        chunk_dict['S']['v'] = np.array([], dtype=np.float64)
    
        return chunk_dict
    
    # Otherwise, we must process the chunks, so first retrieve the desired chunk 
    (t_path, v_path) = chunk_paths
    
    # Load in the metadata and the values vector for this chunk
    t: np.ndarray = np.load(t_path) 
    v: np.ndarray = np.load(v_path) 

    # Next, generate the accelerometer specific t vector from this t vector 
    accelerometer_t: np.ndarray = ms_util.generate_accelerometer_t(t)

    # We will now split between the minispect and the sunglasses data. 
    # The sunglasses data is a 12 bit (unpacked to 16 bit) value stored as 
    # the last 2 bytes of every MS reading 
    m_v: np.ndarray = np.ascontiguousarray(v[:, :v.shape[1]-ms_util.SUNGLASSES_UNCOMPRESSED_FRAME_SHAPE[0]])
    s_v: np.ndarray = np.ascontiguousarray(v[:, -int(ms_util.SUNGLASSES_UNCOMPRESSED_FRAME_SHAPE[0]):]).flatten() #.astype(np.uint16).flatten()

    # Now, send the MS buffers to the minispect util library to parse them. 
    # They are dataframes converted to numpy arrays
    AS_df, TS_df, LS_df, LS_temp_df = ms_util.parse_readings(m_v)

    # Apply the desired transformations 
    if(convert_time_units is True):
        # MS/Sunglasses timestamps are already in seconds, no need to convert 
        pass  
    
    # Apply the time offset correction if desired 
    if(apply_phase_offset is True):
        t += ms_util.MS_TIME_OFFSET

    # No use for converting mean frame as these data are not in the shape of frames 
    # but individual sensor readings 

    # First, let's save the sunglasses information since this is relatively esay 
    chunk_dict['S']['t'] = t 
    chunk_dict['S']['v'] = s_v

    # Next, we will save 
    # Save the finalized values while converting to float if desired
    chunk_dict['M']['t'] = {sensor_name: t if sensor_name != "LS" else accelerometer_t
                            for sensor_name in ms_util.MS_SENSOR_NAMES
                           } 

    chunk_dict['M']['v'] = {sensor_name: df.to_numpy() if convert_to_float is False else df.to_numpy().astype(np.float64)
                            for sensor_name, df in zip(ms_util.MS_SENSOR_NAMES, (AS_df, TS_df, LS_df, LS_temp_df))
                            }

    return chunk_dict

def pupil_chunk_parser(chunk_paths: tuple[str],
                       apply_digital_gain: bool=False, 
                       convert_time_units: bool=False,
                       use_mean_frame: bool=False, 
                       convert_to_float: bool=False,
                       apply_phase_offset: bool=False, 
                       apply_RGB_correction: bool=False, 
                       apply_fielding_function: bool=False,
                       debayer_images: bool=False,
                       mean_axes: tuple[int] = (1, 2),
                       contains_AGC_metadata: bool=False, 
                       password="1234"
                      ) -> dict[str, dict]:
    """Parse one pupil-camera chunk into timestamps, frames, and settings.

    The parser loads the metadata matrix and frame stack, separates the
    timestamp column from any AGC setting columns, and then applies the
    subset of transformations that are meaningful for pupil data: optional
    phase-offset correction, optional digital-gain scaling, optional frame
    averaging, and optional conversion to ``float64``.

    Args:
        chunk_paths: Tuple of (metadata_path, values_path) for the chunk,
            or an empty tuple if no readings exist.
        apply_digital_gain: Whether to apply digital gain scaling per frame.
        convert_time_units: Unused because pupil timestamps are already in
            seconds.
        use_mean_frame: Whether to reduce frames by averaging over axes.
        convert_to_float: Whether to convert values to float64.
        apply_phase_offset: Whether to apply the pupil time offset.
        apply_RGB_correction: Unused for the pupil parser.
        apply_fielding_function: Unused for the pupil parser.
        debayer_images: Unused because pupil frames are not handled as raw
            Bayer images here.
        mean_axes: Axes over which to compute the mean when reducing.
        contains_AGC_metadata: Whether AGC metadata columns are present.
        password: Reserved for future encrypted-chunk support; currently
            unused.

    Returns:
        Dictionary with 'P' key containing 't' (timestamps), 'v' (values),
        and 'settings' (AGC metadata) sub-dictionaries.
    """
    # Initialize return variable
    chunk_dict: dict[str, dict] = {'P': {}}

    # If there were no readings for this sensor 
    if(len(chunk_paths) == 0):
        # initialize t and v as empty numpy arrays.
        chunk_dict['P']['t'] = np.array([], dtype=np.float64)
        chunk_dict['P']['v'] = np.array([], dtype=np.float64)

        # Initialize empty settings dict 
        # Save this as a settings dict 
        chunk_dict['P']['settings'] = { setting: np.array([], dtype=np.float64)
                                        for setting in pupil_util.PUPIL_AGC_METADATA_COLS
                                        }


        return chunk_dict
    
    # Otherwise, we must process the chunks, so first retrieve the desired chunk 
    (t_path, v_path) = chunk_paths

    # Determine if we must decompress + decrypt these files

    # Read the numpy arrays in and extract the t from the metadata 
    metadata: np.ndarray = np.load(t_path) 
    t: np.ndarray = metadata.flatten() if len(metadata.shape) == 1 or metadata.shape[1] == 1 else np.ascontiguousarray(metadata[:, 0])
    v: np.ndarray = np.load(v_path) 
    assert(len(t) == len(v))

    # Extract the AGC metadata (if it exists). 
    chunk_dict['P']['settings'] = { setting: np.array([], dtype=np.float64) if len(metadata.shape) == 1 or metadata.shape[1] == 1 else np.ascontiguousarray(metadata[:, 1+setting_col_num])
                                    for setting_col_num, setting in enumerate(pupil_util.PUPIL_AGC_METADATA_COLS) 
                                    }

    # Apply the desired transformations 
    if(convert_time_units is True):
        # Pupil timestamps are already in seconds, no need to convert 
        pass  

    # Apply the time offset if desired (calculated in seconds)
    if(apply_phase_offset is True):
        t += pupil_util.PUPIL_TIME_OFFSET 

    if(apply_digital_gain is True):
        # Retrieve the digital gain scalars per frame 
        digital_gain_scalars: np.ndarray = chunk_dict['P']['settings']['Dgain'] if len(chunk_dict['P']['settings']['Dgain']) > 0 else np.full((len(t),), 1) 

        # Apply these scalars elementwise to the value vector (frame buffer)
        v = v.astype(np.float64, copy=False) * digital_gain_scalars[:, np.newaxis, np.newaxis]
    
    if(use_mean_frame is True):
        v_mean = np.mean(v, axis=mean_axes)
        del v
        v = v_mean
        del v_mean

    if(v.dtype != np.uint8):
        v = np.clip(v, 0, 255) if use_mean_frame is True else np.round(np.clip(v, 0, 255)).astype(np.uint8)

    if(convert_to_float is True):
        t = t.astype(np.float64)
        v = v.astype(np.float64)

    # Save the finalized values
    chunk_dict['P']['t'] = t
    chunk_dict['P']['v'] = v

    return chunk_dict

def world_chunk_parser(chunk_paths: tuple[str],
                       apply_digital_gain: bool=False,
                       convert_time_units: bool=False,
                       use_mean_frame: bool=False,
                       convert_to_float: bool=False,
                       apply_phase_offset: bool=False,
                       apply_RGB_correction: bool=False,
                       apply_fielding_function: bool=False,
                       debayer_images: bool=False,
                       mean_axes: tuple[int] = (1, 2),
                       contains_AGC_metadata: bool=False,
                       password="1234",
                       differentiate_color: bool=False, 
                       clip_data: bool=False 
                       ) -> dict[str, dict]:
    """Parse one world-camera chunk with optional calibration transforms.

    This parser is the main Python entry point for raw world-camera chunks.
    It loads the timestamp/AGC metadata matrix and Bayer frame stack, then
    conditionally applies timestamp conversion, optional phase correction,
    optional digital gain, optional fielding and Bayer RGB scaling, optional
    debayering, and optional mean reductions. When
    ``differentiate_color=True`` together with ``use_mean_frame=True``, the
    spatial reduction is performed separately for the Bayer ``R``, ``G``,
    and ``B`` pixel classes rather than on the frame as a whole.

    Args:
        chunk_paths: Tuple of (metadata_path, values_path) for the chunk,
            or an empty tuple if no readings exist.
        apply_digital_gain: Whether to apply per-frame digital gain scaling.
        convert_time_units: Whether to convert nanosecond timestamps to
            seconds.
        use_mean_frame: Whether to reduce frames by averaging over axes.
        convert_to_float: Whether to convert values to float64.
        apply_phase_offset: Whether to apply the world camera time offset.
        apply_RGB_correction: Whether to apply per-channel RGB scalar
            corrections on the Bayer pattern.
        apply_fielding_function: Whether to apply the spatial fielding
            function correction.
        debayer_images: Whether to debayer the raw Bayer frames.
        mean_axes: Axes over which to compute the mean when reducing.
        contains_AGC_metadata: Whether AGC metadata columns are present.
        password: Reserved for future encrypted-chunk support; currently
            unused.
        differentiate_color: When True and use_mean_frame is True, computes
            separate per-color-channel means from the Bayer pattern.

    Returns:
        Dictionary with a ``'W'`` entry containing timestamps at
        ``['W']['t']``, processed frame data or reduced statistics at
        ``['W']['v']``, and per-frame AGC arrays at ``['W']['settings']``.
    """

    if(differentiate_color is True):
        assert use_mean_frame is True, f"To differentiate color, the mean frame parameter must be true"

    # Initialize a return value
    chunk_dict: dict = {'W': {}} 
    
    # If there was no chunk paths passed in, 
    # (e.g.) If there were no readings for this sensor 
    # Initialize values as empty. 
    if(len(chunk_paths) == 0):
        # initialize t and v as empty numpy arrays.
        chunk_dict['W']['t'] = np.array([], dtype=np.float64)
        chunk_dict['W']['v'] = np.array([], dtype=np.float64)

        # Initialize empty settings dict 
        # Save this as a settings dict 
        chunk_dict['W']['settings'] = { setting: np.array([], dtype=np.float64)
                                        for setting in world_util.WORLD_AGC_METADATA_COLS
                                        }

        return chunk_dict
    
    # Otherwise, we must process the chunks, so first retrieve the desired chunk 
    (t_path, v_path) = chunk_paths
    
    # Retrieve the metadata, the t specifically and the frame buffer (v)
    metadata: np.ndarray = np.load(t_path)
    t: np.ndarray = metadata.flatten() if len(metadata.shape) == 1 or metadata.shape[1] == 1 else np.ascontiguousarray(metadata[:, 0])
    v: np.ndarray = np.load(v_path)

    # Assert we have metadata for every frame   
    assert len(t) == len(v), "Length of metadata vector must match length of value vector"

    # Extract the AGC metadata (if it exists). 
    chunk_dict['W']['settings'] = { setting: np.array([], dtype=np.float64) if len(metadata.shape) == 1 or metadata.shape[1] == 1 else np.ascontiguousarray(metadata[:, 1+setting_col_num])
                                    for setting_col_num, setting in enumerate(world_util.WORLD_AGC_METADATA_COLS) 
                                    }

    # Apply the desired transformations 
    # Convert time units to seconds 
    if(convert_time_units is True): 
        # World Timestamps are expressed in nanoseconds since system boot. Convert to seconds 
        t = np.ascontiguousarray(t, dtype=np.float64) / (10**9) # 1 sec = 10^9 nanoseconds

        # Apply the time offset if so desired
        # Calculated in seconds, so only can be applied 
        # when convert time units is True
        if(apply_phase_offset is True):
            t += world_util.WORLD_TIME_OFFSET


    # Let's convert to float at the start before we do all of the transformations 
    # if it is not float already 
    v = v.astype(np.float64, copy=False)


    # Apply digital gain as calculated from the AGC
    if(apply_digital_gain is True):
        digital_gain_scalars: np.ndarray = chunk_dict['W']['settings']['AGCDgain'] if len(chunk_dict['W']['settings']['AGCDgain']) > 0 else np.full((len(t),), 1)
        v *= digital_gain_scalars[:, np.newaxis, np.newaxis]

    # Apply the fielding function to each color plane
    if(apply_fielding_function is True):
        fielding_function: np.ndarray = world_util.WORLD_FIELDING_FUNCTIONS[v.shape[1:3]]
        v *= fielding_function

    # Apply correction scalars for each of the RGB channels
    if(apply_RGB_correction is True):
        bayer_RGB_mask: np.ndarray = world_util.generate_RGB_mask(v[0])
        bayer_RGB_pixel_locations: list[np.ndarray] = [ np.argwhere(bayer_RGB_mask == color) for color in "RGB" ]
        assert bayer_RGB_mask.shape == v.shape[1:], f"Bayer RGB mask shape: {bayer_RGB_mask.shape} not equal to frame shape: {v.shape[1:]}"

        assert len(bayer_RGB_pixel_locations) == len(world_util.WORLD_RGB_SCALARS)
        for (pixels, weight) in zip(bayer_RGB_pixel_locations, world_util.WORLD_RGB_SCALARS):
            rows: np.ndarray = pixels[:, 0]
            cols: np.ndarray = pixels[:, 1]
            v[:, rows, cols] *= weight

    # Next, we will debayer the images if desired. Naturally, 
    # this involves clipping and rounding, so could cause (e.g.)
    # the mean to deviate from what you expect 
    if(debayer_images is True):
        v = v if v.dtype == np.uint8 else np.round(np.clip(v, 0, 255)).astype(np.uint8)

        debayered_v: np.ndarray = np.empty(list(v.shape) + [3], dtype=np.uint8)
        for frame_num, frame in enumerate(v):
            debayered_v[frame_num] = world_util.debayer_image(frame)
        v = debayered_v


    # Next, if we want the mean frame, let's do those operations 
    if(use_mean_frame is True):
        # If we we want to differentiate color in the mean, let's do so now 
        if(differentiate_color is True):
            assert tuple(sorted(mean_axes)) == (1, 2), f"To differentiate color, it has to be on a per frame basis"

            # First, let's get the bayer mask for this frame shape
            bayer_mask: np.ndarray = world_util.generate_RGB_mask(np.zeros(v.shape[1:3], dtype=np.uint8), marker="num")

            # Find where each of the individual pixels are
            r_mask: np.ndarray = bayer_mask == 0
            g_mask: np.ndarray = bayer_mask == 1
            b_mask: np.ndarray = bayer_mask == 2

            # Compute per-color and frame means one at a time to avoid
            # materializing large intermediate slices simultaneously
            r_means: np.ndarray = v[:, r_mask].mean(axis=1, dtype=np.float64)
            g_means: np.ndarray = v[:, g_mask].mean(axis=1, dtype=np.float64)
            b_means: np.ndarray = v[:, b_mask].mean(axis=1, dtype=np.float64)
            frame_means: np.ndarray = v.mean(axis=(1, 2), dtype=np.float64)

            v = np.stack([r_means, g_means, b_means, frame_means], axis=1)

        # Otherwise, simply mean the desired axes
        else:
            v = np.mean(v, axis=mean_axes)

    # As we are returning camera data, ensure 
    # it is in the 0-255 range
    if(clip_data is True and v.dtype != np.uint8):
        v = np.round(np.clip(v, 0, 255)).astype(np.uint8)

    # Convert result to float or not
    if(convert_to_float is True):
        t = t.astype(np.float64, copy=False)
        v = v.astype(np.float64, copy=False)

    # Save the finalized values
    chunk_dict['W']['t'] = t
    chunk_dict['W']['v'] = v

    return chunk_dict


def group_sensors_files(recording_path: str) -> dict[str, list[tuple]]:
    """Pair sensor metadata chunks with their matching value files.

    The light-logger recording format stores each chunk as a metadata array
    plus a value file whose name differs only by the presence of the
    ``"_metadata"`` marker. This helper groups those pairs by sensor and
    preserves natural-sort ordering so chunk indices line up with the
    recording timeline.

    Args:
        recording_path: Recording directory to inspect.

    Returns:
        Dictionary keyed by ``"W"``, ``"P"``, and ``"M"`` whose values are
        ordered ``(metadata_path, value_path)`` tuples.
    """
    def group_sensor_files(sensor_name: str) -> list:
        """Collect ordered chunk pairs for one sensor namespace.

        Args:
            sensor_name: Filename fragment identifying the sensor stream.

        Returns:
            List of ``(metadata_path, value_path)`` tuples for that sensor.
        """
        return [ ( os.path.join(recording_path, file), os.path.join(recording_path, file.replace("_metadata", "") ) ) 
                   for file in natsorted(os.listdir(recording_path)) 
                   if sensor_name in file and "metadata" in file
               ]

    return {sensor[0].upper(): group_sensor_files(sensor) # n chunks = [ (metadata_path, frame_buffer_path), ...  ]
            for sensor in ("world", "pupil", "ms")
           }


def find_chunks_by_timestamp(recording_path: str,
                             sensor_name: str, 
                             timestamp_range: tuple[float | None] = (0, None)
                            ) -> tuple[int | None]:
    """Translate a relative time window into the corresponding chunk range.

    The function first establishes the recording's relative time axis from
    the first and last chunk timestamps, then scans the chunk boundaries
    until it finds the chunk containing the requested start time and the
    chunk containing the requested end time. The returned end index follows
    normal Python slice semantics and is therefore exclusive.

    Args:
        recording_path: Path to the recording directory.
        sensor_name: Single-character sensor key ('W', 'P', or 'M').
        timestamp_range: Tuple of (start_seconds, end_seconds) relative
            to recording start. None for end means the recording end.

    Returns:
        ``(start_chunk_index, end_chunk_index)`` with an exclusive end
        bound. Either entry can remain ``None`` if that edge is never found.
    """

    # Retrieve the files for this sensor
    chunks_paths: dict[str, list[tuple]] = group_sensors_files(recording_path)[sensor_name]
    
    # Now, we have to remember that all timestamps in the video are calculated off of 
    # the initial timestamp of the sensor. Therefore, let's read in the first chunk 
    # for this, and the last chunk so we know if the timestamp is in range 
    first_chunk_metadata_matrix_path, first_chunk_v_path = chunks_paths[0]
    first_chunk_metadata_matrix: np.memmap = np.load(first_chunk_metadata_matrix_path, mmap_mode='r')
    last_chunk_metadata_matrix_path, last_chunk_v_path = chunks_paths[-1]
    last_chunk_metadata_matrix: np.memmap = np.load(last_chunk_metadata_matrix_path, mmap_mode='r')

    # Define a variable for time unit conversion 
    conversion_factor: float = (10 ** 9 if sensor_name == "W" else 1)

    # Retrieve the first/last timestamp. (Sensors may optionally have other metadata in the row, so handle accordingly)
    sensor_start_time: float = first_chunk_metadata_matrix[0] if len(first_chunk_metadata_matrix.shape) == 1 else first_chunk_metadata_matrix[0, 0]
    sensor_end_time: float = last_chunk_metadata_matrix[-1] if len(last_chunk_metadata_matrix.shape) == 1 else last_chunk_metadata_matrix[-1, 0]

    # Construct the timestamps of the produced video from these chunks 
    relative_end_time: float = (sensor_end_time - sensor_start_time) / conversion_factor  # World Timestamps are saved in terms of nanoseconds. Convert to seconds 

    # Extract the desired start, end time 
    start, end = timestamp_range
    end = end or relative_end_time

    assert(start >= 0 and start <= relative_end_time)
    assert(end >= 0 and end <= relative_end_time)

    # Now, let's initialize a return variable 
    # and find the chunks 
    chunk_range: tuple[int | None] = [None, None]

    # Now, we will iterate through all the chunks and find a start and end 
    for chunk_num, (chunk_metadata_path, chunk_frame_buffer_path) in enumerate(chunks_paths):
        # Peer into the chunk 
        chunk_metadata_matrix: np.memmap = np.load(chunk_metadata_path, mmap_mode='r')

        # Find the start and end of this chunk 
        chunk_start: float = chunk_metadata_matrix[0] if len(chunk_metadata_matrix.shape) == 1 else chunk_metadata_matrix[0][0]
        chunk_end: float = chunk_metadata_matrix[-1] if len(chunk_metadata_matrix.shape) == 1 else chunk_metadata_matrix[-1][0]

        # Convert the time to relative timings 
        relative_chunk_start: float = (chunk_start - sensor_start_time) / conversion_factor
        relative_chunk_end: float = (chunk_end - sensor_start_time) / conversion_factor

        # Determine if it is our start 
        if(start >= relative_chunk_start and start <= relative_chunk_end):
            chunk_range[0] = chunk_num
        
        # Determine if it is our end 
        if(end >= relative_chunk_start and end <= relative_chunk_end):
            chunk_range[-1] = chunk_num + 1

    return tuple(chunk_range) 

def parse_chunks(recording_path: str,
                 apply_digital_gain: bool=False,
                 use_mean_frame: bool=False,
                 convert_time_units: bool=False,
                 convert_to_float: bool=False,
                 apply_phase_correction: bool=False,
                 apply_RGB_correction: bool=False,
                 apply_fielding_function: bool=False,
                 debayer_images: bool=False,
                 time_ranges: tuple[float | None] | None={sensor_name: (None, None)
                                                          for sensor_name in "WPM"
                                                         },
                 chunk_ranges: tuple[int] | None={sensor_name: (0, None)
                                                  for sensor_name in "WPM"
                                                 },
                 mean_axes: dict[str, int] = {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                 contains_agc_metadata_dict: dict[str, bool] = {'W': True, 'P': False, 'M': False},
                 password: str="1234",
                 differentiate_color: bool=False,
                 verbose: bool=False
                ) -> list[dict[str, dict]]:
    """Parse an entire recording across all supported sensor streams.

    This is the orchestration layer above the individual sensor parsers. It
    groups files by sensor, resolves either explicit chunk ranges or
    timestamp-derived chunk ranges, iterates across the union of chunk
    positions present in the recording, and asks the appropriate parser to
    decode each sensor at each chunk index. When a sensor has no data for a
    selected chunk position, that parser's empty-output convention is used
    so the output structure remains regular.

    Args:
        recording_path: Path to the recording directory.
        apply_digital_gain: Whether to apply digital gain correction.
        use_mean_frame: Whether to reduce frames by averaging.
        convert_time_units: Whether to convert timestamps to seconds.
        convert_to_float: Whether to convert values to float64.
        apply_phase_correction: Whether to apply sensor phase offsets.
        apply_RGB_correction: Whether to apply RGB scalar corrections.
        apply_fielding_function: Whether to apply the fielding function.
        debayer_images: Whether to debayer Bayer-pattern frames.
        time_ranges: Optional per-sensor relative time windows used to
            select chunks.
        chunk_ranges: Optional per-sensor explicit chunk-index ranges.
        mean_axes: Per-sensor axes passed through when frame means are
            requested.
        contains_agc_metadata_dict: Per-sensor flags describing whether AGC
            columns exist in the metadata arrays.
        password: Reserved for future encrypted-chunk support.
        differentiate_color: Whether to compute per-color-channel means
            for the world camera.
        verbose: Whether to display a progress bar.

    Returns:
        List of per-chunk dictionaries containing the parsed outputs from
        the sensor-specific parsers.
    """

    # Ensure the time_ranges and chunk ranges are tuples of ints 
    for parameter_dict in (time_ranges, chunk_ranges):
        for sensor, range_tuple in parameter_dict.items():
            parameter_dict[sensor] = tuple([ None if val is None else int(val) for val in range_tuple])

    # Ensure the mean axes are tuples of ints if not single numbers 
    for sensor, axes in mean_axes.items():
        if(isinstance(axes, float) or isinstance(axes, int)): continue 
        mean_axes[sensor] = tuple( [ int(axis) for axis in axes ] )

    # Retrieve each of the sensor's files grouped together 
    sensors_chunks_paths: dict[str, list[tuple]] = group_sensors_files(recording_path)

    # Let's now calculate the max number of chunks 
    max_num_chunks: int = max( [ len(chunk_paths) 
                                 for (sensor_name, chunk_paths) 
                                 in sensors_chunks_paths.items() ] 
                                 + [ 0 ]  
                             )
    
    # Assert we have at one of the sensors has readings
    assert max_num_chunks > 0, "No chunks found for any sensor"

    ## Load in the config file for the recording
    #config_data: dict | None = None
    #with open(os.path.join(recording_path, "config.pkl"), 'rb') as f:
    #    config_data = dill.load(f)
 
    # Save a dictionary of the finalized splices of chunks 
    selected_sensor_ranges: dict[str, tuple[int]] = chunk_ranges.copy() 

    # Now, let's iterate over the sensors 
    for sensor_name in "WPM":
        # Let's retrieve the chunks for this sensor 
        sensor_chunks: list[tuple] = sensors_chunks_paths[sensor_name]

        # Let's also retrieve the desired time range or chunks to use 
        time_range: tuple[float | None] = tuple(time_ranges[sensor_name]) if sensor_name in time_ranges else (None, None)
        chunk_range: tuple[int | None] = tuple(chunk_ranges[sensor_name]) if sensor_name in chunk_ranges else (None, None)
        
        # If the user demanded a range that is out of bounds 
        # throw a helpful error 
        if(chunk_range[-1] is not None and chunk_range[-1] > len(sensor_chunks)): 
            raise Exception(f"ERROR: Chunk range: {chunk_range} is out of bounds for sensor: {sensor_name} with num chunks: {len(sensor_chunks)}")

        # We will only allow selection using a single method 
        # Values compared here are the defaults for each selection type
        if(time_range != (None, None) and chunk_range not in ((None, None), (0, None))):
            raise Exception(f"ERROR: Sensor: {sensor_name} Can only specify by either time or chunk number, not both")

        # If we entered a time range and want this sensor 
        # let's gather the chunk_range through the timestamps 
        if(time_range != (None, None) and chunk_range != (None, None)):
            chunk_range = find_chunks_by_timestamp(recording_path, sensor_name, time_range)

        # Save this finalized chunk range 
        selected_sensor_ranges[sensor_name] = chunk_range

    # Initialize the return list of parsed chunk dicts. 
    parsed_chunks: list[dict[str, dict]] = []

    # Define a dictionary for the sensors and their associated parsers 
    parsers: dict = {sensor_name: parser 
                     for sensor_name, parser 
                     in zip("WPM", (world_chunk_parser,
                                           pupil_chunk_parser,
                                           minispect_and_sunglasses_chunk_parser
                                          )
                           )
                    }

    # Populate the chunks
    chunk_iterator = range(max_num_chunks) if verbose is False else tqdm(range(max_num_chunks), desc="Parsing chunks")
    for chunk_num in chunk_iterator:
        # Initialize the dict for this chunk 
        chunk_dict: dict = {}

        # Populate the chunk dict. Note: Use sensors_chunks_paths dict
        # here instead of chunk dict because the sunglasses data is stored in the MS data, 
        # so I just want to iterate over W, P, M
        for sensor_name in sensors_chunks_paths.keys():    
            # Determine if this sensor should have a value 
            # for this chunk
            start, end = selected_sensor_ranges[sensor_name]
            end = end or len(sensors_chunks_paths[sensor_name])

            # Gather the metadata and frame buffer to send to 
            # the parser.
            chunk: tuple[str] = () if start is None or chunk_num not in range(start, end) else (sensors_chunks_paths[sensor_name][chunk_num])

            # Parse the chunk and save it in the dictionary
            sensor_kwargs: dict = {}
            if sensor_name == "W":
                sensor_kwargs["differentiate_color"] = differentiate_color

            parsed_sensor_dict: dict = parsers[sensor_name](chunk,
                                                            apply_digital_gain,
                                                            convert_time_units,
                                                            use_mean_frame,
                                                            convert_to_float,
                                                            apply_phase_correction,
                                                            apply_RGB_correction,
                                                            apply_fielding_function,
                                                            debayer_images,
                                                            mean_axes[sensor_name],
                                                            contains_agc_metadata_dict[sensor_name],
                                                            password,
                                                            **sensor_kwargs
                                                           )
            chunk_dict |= parsed_sensor_dict
            gc.collect()

        # Append the parsed chunk dict
        parsed_chunks.append(chunk_dict)

    return parsed_chunks
