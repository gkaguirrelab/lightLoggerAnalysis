import numpy as np 
import os 
import sys
from typing import Iterable
import natsort


# Import relevant custom libraries with helper functions and constants 
light_logger_analysis_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLoggerAnalysis")
video_io_util_path: str = os.path.join(light_logger_analysis_dir_path, "code", "library", "matlabIO", "python_libraries")
for path in (light_logger_analysis_dir_path, video_io_util_path):
    assert os.path.exists(path), f"Expected path: {path} does not exist"
    sys.path.append(path)

import video_io 

"""Find the start times of each of the sensors in the recording"""
def find_sensor_start_end_times(path_to_recording: str) -> dict[str, tuple]:
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

"""Convert world camera frame number to pupil camera frame number"""
def world_to_pupil(world_frame_numbers: Iterable, 
                   path_to_recording_folder: str, 
                   path_to_recording_chunks: str,
                   pupil_phase_offset: float=0.05
                  ) -> np.ndarray:

    # Read in the first and last timestamp of the pupil and world cameras 
    sensor_start_ends: dict[str, tuple] = find_sensor_start_end_times(path_to_recording_chunks)
    
    # Generate timestamps for the two videos
    path_to_world_video: str = os.path.join(path_to_recording_folder, "W.avi")
    path_to_pupil_video: str = os.path.join(path_to_recording_folder, "P.avi")

    world_t: np.ndarray = np.linspace(*sensor_start_ends["world"], video_io.inspect_video_frame_count(path_to_world_video), endpoint=True)
    pupil_t: np.ndarray = np.linspace(*sensor_start_ends["pupil"], video_io.inspect_video_frame_count(path_to_pupil_video), endpoint=True) + pupil_phase_offset

    pupil_frame_indices: list = []
    for world_frame_number in world_frame_numbers:
        # Find the nearest neighbor in pupil camera 
        time_differences_from_target: np.ndarray = np.abs(world_t[int(world_frame_number)] - pupil_t)
        pupil_index: int = int(np.argmin(time_differences_from_target))
        
        # Append it to the list 
        pupil_frame_indices.append(pupil_index)

    return np.array(pupil_frame_indices)


"""Convert pupil camera frame number to world camera frame number"""
def pupil_to_world() -> np.ndarray:
    pass 


def main():
    pass 

if(__name__ == "__main__"):
    main() 