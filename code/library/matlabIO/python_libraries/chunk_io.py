import os 
from natsort import natsorted
import numpy as np
import pathlib
import sys

# Import the sensor utility libraries 
light_logger_analysis_path: str = str(pathlib.Path(__file__).parents[3])
sensor_utility_path: str = os.path.join(light_logger_analysis_path, "library", "sensor_utility")
for path in (sensor_utility_path,):
    sys.path.append(path)

import world_util 
import ms_util 
import pupil_util

"""Define the parser for the MS and sunglasses data"""
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

"""Define the parser for the pupil camera data"""
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
    v: np.ndarray = np.load(v_path) if os.path.splitext(v_path)[1]
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
        v = np.clip(np.round(v.astype(np.float64) * digital_gain_scalars[:, np.newaxis, np.newaxis]), 0, 255).astype(np.uint8)
    
    if(use_mean_frame is True):
        v = np.mean(v, axis=mean_axes)

    if(convert_to_float is True):
        t = t.astype(np.float64)
        v = v.astype(np.float64) 

    # Save the finalized values 
    chunk_dict['P']['t'] = t 
    chunk_dict['P']['v'] = v

    return chunk_dict

"""Define the individual per chunk parser for the world camera data"""
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
                       password="1234"
                       ) -> dict[str, dict]:
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

    # First, convert to float for float multiplications
    # if not in float already
    v = v.astype(np.float64, copy=False)

    # Apply digital gain as calculated from the AGC 
    if(apply_digital_gain is True):
        assert v.dtype == np.float64, f"Frame buffer dtype must be float64 to apply dgain. Current dtype is {v.dtype}"

        # Retrieve the digital gain scalars per frame 
        digital_gain_scalars: np.ndarray = chunk_dict['W']['settings']['Dgain'] if len(chunk_dict['W']['settings']['Dgain']) > 0 else np.full((len(t),), 1) 

        # Multiply the frames by the dgains
        v *= digital_gain_scalars[:, np.newaxis, np.newaxis]

    # Apply the fielding function to each color plane 
    if(apply_fielding_function is True):
        assert v.dtype == np.float64 and v.ndim == 3, f"Frame buffer dtype must be float64 and in bayer format to apply fielding function. Current dtype is {v.dtype}, current ndim: {v.ndim}"
        
        # Get the fielding function for this frame size
        fielding_function: np.ndarray = world_util.WORLD_FIELDING_FUNCTIONS[v.shape[1:3]]

        # Apply fielding function
        v *= fielding_function

    # Apply correction scalars for each of the RGB channels
    if(apply_RGB_correction is True):
        assert v.ndim == 3 and v.dtype == np.float64, f"To apply color weights, buffer must be in bayer form np.float64."

        # Apply the radiometric correction
        bayer_RGB_mask: np.ndarray = world_util.generate_RGB_mask(v[0])
        bayer_RGB_pixel_locations: list[np.ndarray] = [ np.argwhere(bayer_RGB_mask == color) for color in "RGB" ]
        assert bayer_RGB_mask.shape == v.shape[1:], f"Bayer RGB mask shape: {bayer_RGB_mask.shape} not equal to frame shape: {v.shape[1:]}" 

        # Apply the color corrections IN-PLACE and vectorized
        assert len(bayer_RGB_pixel_locations) == len(world_util.WORLD_RGB_SCALARS)
        for (pixels, weight) in zip(bayer_RGB_pixel_locations, world_util.WORLD_RGB_SCALARS):
            rows: np.ndarray = pixels[:, 0]
            cols: np.ndarray = pixels[:, 1]

            # Apply the weight to the specified pixels 
            v[:, rows, cols] *= weight

    # Convert back to np.uint8 to debayer images
    if(v.dtype != np.uint8):
        v = np.round(np.clip(v, 0, 255)).astype(np.uint8)

    #  Debayer images if desired 
    if(debayer_images is True):
        debayered_v: np.ndarray = np.empty(list(v.shape) + [3], dtype=np.uint8)
        for frame_num, frame in v:
            debayered_v[frame_num] = world_util.debayer_image(frame)

        # Replace v with the debayered version
        v = debayered_v

    # Just return the mean frame 
    if(use_mean_frame is True):
        v = np.mean(v, axis=mean_axes)

    # Convert result to float or not 
    if(convert_to_float is True):
        t = t.astype(np.float64)
        v = v.astype(np.float64) 

    # Save the finalized values 
    chunk_dict['W']['t'] = t 
    chunk_dict['W']['v'] = v

    return chunk_dict


"""Inspect a recording directory and return a dictionary of sensor names
   with all of the paths to the chunks for that sensor in order
"""
def group_sensors_files(recording_path: str) -> dict[str, list[tuple]]:    
    """Find all of the chunk filepaths of the sensors and group them into tuples of (t, frames)
    Do this by iterating over the files over the recording in numerically sorted order (to order the chunks properly),
    then find the files that denote they are the metadata of the sensors. Then, simply remove the metadata part of 
    the filename and you arrive at the value file"""
    def group_sensor_files(sensor_name: str) -> list:
        return [ ( os.path.join(recording_path, file), os.path.join(recording_path, file.replace("_metadata", "") ) ) 
                   for file in natsorted(os.listdir(recording_path)) 
                   if sensor_name in file and "metadata" in file
               ]

    return {sensor[0].upper(): group_sensor_files(sensor) # n chunks = [ (metadata_path, frame_buffer_path), ...  ]
            for sensor in ("world", "pupil", "ms")
           }


"""Given a timestamp range from a sensor video, find the range of  
   chunks that compose this range
"""
def find_chunks_by_timestamp(recording_path: str, 
                             sensor_name: str, 
                             timestamp_range: tuple[float | None] = (0, None)
                            ) -> tuple[int | None]:
    
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
                                                          for sensor_name in SENSOR_NAMES
                                                         }, 
                 chunk_ranges: tuple[int] | None={sensor_name: (0, None)
                                                  for sensor_name in SENSOR_NAMES
                                                 }, 
                 mean_axes: dict[str, int] = {'W': (1, 2), 'P': (1, 2), 'M': (0,)},
                 contains_agc_metadata_dict: dict[str, bool] = {'W': True, 'P': False, 'M': False},
                 password: str="1234"
                ) -> list[dict[str, dict]]:

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
    for sensor_name in SENSOR_NAMES:
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
                     in zip(SENSOR_NAMES, (world_chunk_parser,
                                           pupil_chunk_parser,
                                           minispect_and_sunglasses_chunk_parser
                                          )
                           )
                    }

    # Populate the chunks
    for chunk_num in range(max_num_chunks):
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
                                                            password
                                                           )
            chunk_dict |= parsed_sensor_dict

        # Append the parsed chunk dict 
        parsed_chunks.append(chunk_dict)


    return parsed_chunks
