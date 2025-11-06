import numpy as np 
import natsort
import os 

"""Find the start times of each of the sensors in the recording"""
def find_sensor_start_end_times(path_to_recording: str) -> tuple:
    # Group metadata files by sensor 
    sensor_files: dict[str, tuple[str]] = {sensor: [os.path.join(path_to_recording, file) 
                                                    for file in natsort.natsorted(os.listdir(path_to_recording))
                                                    if file.startswith(sensor) and "metadata" in file
                                                   ] 
                                           for sensor in ("ms", "pupil", "world")
                                          }
    

    # Initialize the output dict 
    sensor_start_ends: dict[str, tuple] = {sensor: (-1, -1)
                                           for sensor in sensor_files
                                          }

    # Iterate over the sensors and fill in the start end time
    # if they existed
    for sensor, files in sensor_files.items():
        # Skip ununused snesors 
        if(len(files) == 0): 
            continue 

        # Load in the first and last files 
        start_chunk_path, end_chunk_path = files[0], files[-1]
        start_chunk: np.ndarray = np.load(start_chunk_path) 
        end_chunk: np.ndarray = np.load(end_chunk_path)

        # Find the start and end time 
        start: float = start_chunk[0] if start_chunk.ndim == 1 else start_chunk[0, 0]
        end: float = end_chunk[-1] if end_chunk.ndim == 1 else end_chunk[-1, 0]
        
        # Define the start by 10**9 to covnvert to seconds on the wolrd camera 
        if(sensor == "world"):
            start = start / (10 ** 9)
            end = end / (10 ** 9)

        # Save the appropriate start/end
        sensor_start_ends[sensor] = (start, end)

    return sensor_start_ends


def main():
    pass 

if(__name__ == "__main__"):
    pass 