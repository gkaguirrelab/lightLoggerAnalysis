import os
import numpy as np 
from natsort import natsorted 
from tqdm.auto import tqdm
from typing import Iterable
import sys
import matplotlib.pyplot as plt
import pandas as pd

"""Define Characteristics of the minispectrometer device"""
MS_FPS: int = 1 # Capture at APPROXIMATELY 1 FPS. The MS is more variable than others and does not maintain as tight of an FPS
MS_COMPRESSED_FRAME_SHAPE: np.ndarray = np.array([1552], dtype=np.uint16) # MS data is read in packets of 148 bytes 
MS_COM_PORT: str = '/dev/ttyACM0' if sys.platform.startswith('linux') else '/dev/tty.usbmodem142201' # Port nomenclature changes between operating systems 
MS_BAUDRATE: int = 115200 # Baudrate for communication with the on-board Arduino 
MS_AS_CHANNELS: int = 10 # Number of channels in the MS AS sensor readings 
MS_AS_DTYPE: object = np.uint16 # Define the type of the AS sensor readings
MS_TS_CHANNELS: int = 2 # Number of channels in the MS TS sensor readings
MS_TS_DTYPE: object = np.uint16 # Define the type of the TS sensor readings 
MS_LS_CHANNELS: int = 6 # Number of channels from the LS chip (X, Y, R) for
                         # acceleration and (ΩP, ΩR, ΩY) for rotation 
MS_LS_DTYPE: object = np.int16 # Define the data type of the MS accelerometer readings
MS_LS_BUFFER_SIZE: int = 127  # Number of readings from all LS channels in a single MS packet/reading
MS_TEMP_CHANNELS: int = 1 # Number of channels in the MS temp sensor readings
MS_TEMP_DTYPE: object = np.float32 # Define the data type of the MS temperature readings
MS_UNCOMPRESSED_FRAME_SHAPE: np.ndarray = np.array([ MS_AS_CHANNELS + 
                                                     MS_TS_CHANNELS + 
                                                     (MS_LS_CHANNELS * MS_LS_BUFFER_SIZE) +
                                                     MS_TEMP_CHANNELS
                                                   ], 
                                                    dtype=np.uint16
                                                   ) # MS data is then uncompressed and parsed into numpy as 73 numbers
MS_MODE_DICT: dict[str, str] = {"wait": "W", # Define a mapping between semantic names for modes and their arduino names
                                "calibration": "C", 
                                "science": "S"
                                } 

# Define the list of the MS sensor names
MS_SENSOR_NAMES: list[str] = ["AS", "TS", "LS", "TEMP"]

# Define a mapping between names of sensor and the number of channels they have 
MS_NAMES_AND_CHANNELS: dict[str, int] = {sensor_name: num_channels
                                         for sensor_name, num_channels 
                                         in zip(MS_SENSOR_NAMES, (MS_AS_CHANNELS, MS_TS_CHANNELS, MS_LS_CHANNELS, MS_TEMP_CHANNELS))
                                        }


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""Define Characteristics of the sunglasses detector device"""
SUNGLASSES_SENSOR_ADDR: bytes = 0x6b # The address where the hall sunglasses sensor exists. This is what we will read from to get readings 
SUNGLASSES_UNCOMPRESSED_FRAME_SHAPE = np.array([1], dtype=np.uint16) # The sunglasses shape will be 2 bytes (as it is a 12 bit number). 
                                                                    # We will define nothing more in regards to this, as we will simply attach its readings to the MS 
                                                                    # readings. We will parse it first into its 16 bit form, so it takes up one slot
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""Define Characteristics of the merged readings, that is, the datatype and shape"""
MS_SUNGLASSES_FRAME_DTYPE: object = np.float32  # Define the datatype of the saved readings from the MS and the sunglasses merged together 
MS_SUNGLASSES_FRAME_SHAPE: np.ndarray = MS_UNCOMPRESSED_FRAME_SHAPE + SUNGLASSES_UNCOMPRESSED_FRAME_SHAPE # Define the shape of the MS readings and the sunglasses readings merged 
MS_SUNGLASSES_METADATA_DTYPE: object = np.float64 # Define the datatype of the merged sensors' metadata (timestamps)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""
    Generate a t vector for the accelerometer readings. 
    We cannot use the default t vector associated with the MS 
    buffer reading metadata because accelerometer readings 
    are sent as buffers of n readings. 
"""
def generate_accelerometer_t(t: np.ndarray | list | tuple) -> np.ndarray:
    # Convert to numpy array if passed in a list, for instance, 
    # like from MATLAB 
    if(isinstance(t, np.ndarray) is False):
        t = np.array(t)
    
    # If t is empty, just return t 
    if(len(t) == 0):
        return t

    # Ensure that the input vector is flattened 
    t = t.flatten()

    # Define an empty array that is MS_LS_BUFFER_SIZE times the original 
    LS_t: np.ndarray = np.empty(MS_LS_BUFFER_SIZE * len(t), dtype=MS_SUNGLASSES_METADATA_DTYPE)

    # Define a variable to keep track of where we should be inserting values 
    insertion_pointer: int = 0 

    # By the time the first timestamp arrives, there are already MS_LS_BUFFER_SIZE 
    # readings that all called that point in time. This is also true for every other timestamp
    # as well. That first start time is inherently impossible to know, so we can approximate it 
    # by taking the median of all the other ones 
    median_gap_length: float = np.median(np.diff(t))

    # We will fill in the first gap using the median gap length 
    gap_t: np.ndarray = np.linspace(start=t[0] - median_gap_length, stop=t[0],
                                    num=MS_LS_BUFFER_SIZE).flatten()
    
    # Save this filled gap into the newly built t vector 
    LS_t[insertion_pointer: insertion_pointer + MS_LS_BUFFER_SIZE ] = gap_t
    insertion_pointer += MS_LS_BUFFER_SIZE

    # Next, we will iterate over the timestamps we do have 
    for timestamp_num in range(1, len(t)):
        # We will retrieve the current and last reading timestamps. 
        # This will give us the gap that we will evenly fill in 
        current_timestamp: np.float64 = t[timestamp_num]
        previous_timestamp: np.float64 = t[timestamp_num - 1]

        # We will now generate evenly spaced timestamps to fill in this gap 
        gap_t: np.ndarray = np.linspace(start=previous_timestamp, stop=current_timestamp,
                                        num=MS_LS_BUFFER_SIZE).flatten()

        # Save this filled gap into the newly built t vector 
        LS_t[insertion_pointer: insertion_pointer + MS_LS_BUFFER_SIZE] = gap_t
        insertion_pointer += MS_LS_BUFFER_SIZE

    # Then return the correct timestamps 
    return LS_t


"""Plot a single recording buffer of the MS (nSeconds, MS_SUNGLASSES_FRAME_SHAPE)
   onto a MATPLOT lib axis and return said axis. 
"""
def plot_chunk_readings(t_vector: np.ndarray, readings_buffer: np.ndarray, 
                        axes: plt.Axes=None
                       ) -> None | tuple:
    # First, we will generate a figure + associated axis to plot 
    # if this has not been passed in 
    fig: None | plt.figure = None
    if(axes is None):
        fig, axes = plt.subplots(2, 2, figsize=(10,10)) # Size is 2,2 because we have 4 sensors 

    # Next, we will parse the value (v) array into the appropriate chips 
    # and channels 
    chips_and_channels: dict[str, np.ndarray] = {chip_name: channels_df
                                                 for chip_name, channels_df 
                                                 in zip(("AS", "TS", "LS", "TEMP"),
                                                        parse_readings(readings_buffer)
                                                       ) 
                                                }

    # Next, we will plot all of the different chips onto the axes 
    for chip_num, (chip_name, chip_channels_df) in enumerate(chips_and_channels.items()):
        # Convert the channels to a numpy array for simplicity, but also retrieve the column names 
        chip_channels: np.ndarray = chip_channels_df.to_numpy()
        channel_labels: list = chip_channels_df.columns

        # Retrieve the axes to plot on 
        ax = axes[chip_num % 2, chip_num // 2]

        # Label the graph 
        ax.set_title(chip_name, fontweight="bold")
        ax.set_ylabel("Count")
        ax.set_xlabel("Time [s]")

        # Unpack the t vector if needed 
        t: np.ndarray = t_vector if chip_name != "LS" else generate_accelerometer_t(t_vector)

        # Next, plot the channels of the chip onto the graph 
        for channel_num in range(chip_channels.shape[1]):
            # Plot the channel with its associated label
            ax.plot(t, chip_channels[:, channel_num], label=channel_labels[channel_num])


        # Display a legend for the plot 
        ax.legend()

        # Display the plot 
        plt.show() 

    # If we passed in axes, no need to return. 
    if(not fig): return None 

    # Otherwise, return the axes and figure 
    return fig, axes

"""After receiving a packet of data from the MS over serial, 
   parse it into raw buffers for each of the sensors 
"""
def parse_packet(packet: bytes) -> tuple:
    AS_channels: np.ndarray = np.frombuffer(packet[0:20], dtype=np.uint16)
    TS_channels: np.ndarray = np.frombuffer(packet[20:24], dtype=np.uint16)
    LS_channels: np.ndarray = np.frombuffer(packet[24 : len(packet) - 4 ], dtype=np.int16)
    TEMP_channels: np.ndarray = np.frombuffer(packet[-4:], dtype=np.float32)

    return AS_channels, TS_channels, LS_channels, TEMP_channels

"""Parse a single recording buffer of MS readings (nSeconds, MS_SUNGLASSES_FRAME_SHAPE) and return them as a tuple of
   either numpy arrays or pd.DataFrames. Returning as numpy can sometimes help crashing in MATLAB
   as there can be a version mismatch of MATLAB pandas and Python libraries."""
def parse_readings(readings_buffer: np.ndarray) -> tuple:
    # First, we must parse the np.array of bytes into the respective sensors' bytes
    AS_channels, TS_channels, LS_channels, LS_temp = parse_SERIAL(readings_buffer)

    # Parse the readings into dataframes and return them
    AS_df: pd.DataFrame = reading_to_df(AS_channels, MS_AS_DTYPE, sensor_name='A')
    TS_df: pd.DataFrame = reading_to_df(TS_channels, MS_TS_DTYPE, sensor_name='T')
    LS_df: pd.DataFrame = reading_to_df(LS_channels, MS_LS_DTYPE, sensor_name='L')
    LS_temp_df: pd.DataFrame = reading_to_df(LS_temp, MS_TEMP_DTYPE, sensor_name='c')

    return AS_df, TS_df, LS_df, LS_temp_df

"""Extract the different sensor's information out of the readings buffer"""
def parse_SERIAL(readings_buffer: np.ndarray) -> tuple:
    # Define a variable to keep track of where we are in the buffer 
    current_buffer_pos: int = 0 
    
    # Splice out the values and return them to their intended types 
    AS_channels: np.ndarray = np.ascontiguousarray(readings_buffer[:, current_buffer_pos:current_buffer_pos+MS_AS_CHANNELS])
    AS_channels = AS_channels.astype(MS_AS_DTYPE)
    current_buffer_pos += MS_AS_CHANNELS

    TS_channels: np.ndarray = np.ascontiguousarray(readings_buffer[:, current_buffer_pos:current_buffer_pos+MS_TS_CHANNELS])
    TS_channels = TS_channels.astype(MS_TS_DTYPE)
    current_buffer_pos += MS_TS_CHANNELS

    LS_channels: np.ndarray = np.ascontiguousarray(readings_buffer[:, current_buffer_pos: current_buffer_pos + (MS_LS_CHANNELS * MS_LS_BUFFER_SIZE)])
    LS_channels = LS_channels.astype(MS_LS_DTYPE)
    current_buffer_pos += MS_LS_CHANNELS

    LS_temp: np.ndarray = np.ascontiguousarray(readings_buffer[:, current_buffer_pos:current_buffer_pos+MS_TEMP_CHANNELS])
    LS_temp = LS_temp.astype(MS_TEMP_DTYPE)

    return AS_channels, TS_channels, LS_channels, LS_temp 


"""Parse a reading buffer and return the resulting dataframe with labeled cols"""
def reading_to_df(readings: np.ndarray, 
                  channel_type : type, 
                  sensor_name: str
                 ) -> pd.DataFrame:
    
    """Unpack a buffer of accelerometer buffers. We will have (nRows, 60 cols)
       where these 60 cols are buffers of 10 accelerometer readings from the 6 channels
    """
    def unpack_accel_readings(buffer_of_buffers: np.ndarray) -> np.ndarray:
        # Allocate a new array that is the correct shape 
        unpacked_buffer: np.ndarray = np.full((buffer_of_buffers.shape[0] * MS_LS_BUFFER_SIZE, MS_LS_CHANNELS), -1, dtype=np.int16)
        assert(buffer_of_buffers.shape[1] == MS_LS_CHANNELS * MS_LS_BUFFER_SIZE)

        # Iterate over the original buffer
        current_reading_num: int = 0 
        for packet_num in range(buffer_of_buffers.shape[0]):    
            # For a given packet, let's extract the 10 readings
            # from all the channels. They are stored in the format 
            # [r1a, r1b...r2a, r2b...]
            for channel_num in range(MS_LS_CHANNELS):
                # Assign 10 readings to the given channel by starting at the current channel num then counting by the number of channels in a reading until the end 
                unpacked_buffer[current_reading_num : current_reading_num + MS_LS_BUFFER_SIZE, channel_num] = buffer_of_buffers[packet_num, channel_num::MS_LS_CHANNELS]

            # Increment the current reading number by the number of readings we just filled in 
            current_reading_num += MS_LS_BUFFER_SIZE 

        return unpacked_buffer

    # Determine if this is the accelerometer readings or not 
    is_accelerometer: bool = sensor_name[0] == 'L'

    # Each packet from the MS included a buffer of accelerometer values, so if we are parsing 
    # the accelerometer, we need to unpack it first 
    readings = readings if is_accelerometer is False else unpack_accel_readings(readings)

    # Create a CSV from the numpy array of readings
    df: pd.DataFrame = pd.DataFrame(readings)
    
    # Form column names and their associated types for the readings   
    columns : list =  [ str(i) for i in range(df.shape[1]) ] if is_accelerometer is False else ['X', 'Y', 'Z', 'ΩP', 'ΩR', 'ΩY']
    types :list = [ channel_type for i in range(df.shape[1]) ]
    
    # Reformat the DataFrame with col names and types
    df.columns = columns 
    df = df.astype({col:type_ for col, type_ in zip(columns, types)})

    # Return the DF
    return df 

def ms_data_from_chunks(recording_path: str,
                        verbose: bool=False
                       ) -> np.ndarray:
    # First, we need to find all of the MS chunks 
    ms_data_chunks: list[str] = natsorted([os.path.join(recording_path, filename)
                                           for filename in os.listdir(recording_path)
                                           if filename.startswith("ms")
                                           and "metadata" not in filename
                                          ]
                                        )
    ms_metadata_chunks: list[str] = natsorted([os.path.join(recording_path, filename)
                                           for filename in os.listdir(recording_path)
                                           if filename.startswith("ms")
                                           and "metadata" in filename
                                          ]
                                        )
    assert len(ms_data_chunks) > 0, f"0 MS Chunks found @ {recording_path}"
    assert len(ms_data_chunks) == len(ms_metadata_chunks), f"Recording @ {recording_path} has different number of chunks/metadata"

    # Initialize the output variables 
    # that will have all the chunk data joined 
    data: list[np.ndarray] | np.ndarray = []
    metadata: list[np.ndarray] | np.ndarray = []

    # Once we have the paths to the MS chunks, we will simply read them in 
    path_iterator: Iterable = range(len(ms_data_chunks)) if verbose is False else tqdm(range(len(ms_data_chunks)), desc="Loading MS chunks")
    for chunk_num in path_iterator:
        # Retrieve the path to this chunk 
        data_chunk_path: str = ms_data_chunks[chunk_num]
        metadata_chunk_path: str = ms_metadata_chunks[chunk_num]

        # Parse the buffer, selecting only the AS chip
        data_chunk: np.ndarray = parse_readings(data_chunk_path)[0].to_numpy() 
        metadata_chunk: np.ndarray = np.load(metadata_chunk_path)

        # Save this data in the running ararys 
        data.append(data_chunk)
        metadata.append(metadata_chunk)

    # Convert to np.ndarray 
    data = np.array(data)
    metadata = np.array(metadata)

    return data, metadata

def main():
    pass 

if(__name__ == "__main__"):
    main()