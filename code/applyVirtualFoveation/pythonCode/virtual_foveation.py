import numpy as np 
import natsort
import os 
import sys
import matplotlib.pyplot as plt
import subprocess
from typing import Literal

egocentric_video_mapper_path: str = os.path.join(os.path.dirname(__file__), "egocentric_video_mapper")
assert os.path.exists(egocentric_video_mapper_path), "Could not find egocentric video mapper. Ensure it is cloned"
sys.path.append(egocentric_video_mapper_path)
from pupil_labs.egocentric_video_mapper.utils import show_videos_preview

import pupil_labs.egocentric_video_mapper.__main__ as egocentric_main


class Args:
    def __init__(self, 
                 neon_timeseries_dir: str, 
                 alternative_vid_path: str, 
                 output_dir: str, 
                 mapping_choice: Literal["Fixations", "Gaze", "Both"], 
                 optic_flow_algorithm: Literal["Lucas-Kanade", "Gunnar Farneback"], 
                 image_matcher: Literal["Efficient_LOFTR", "LOFTR_indoor"],
                 refresh_time_threshold_sec: float, 
                 render_video_comparison: bool, 
                 render_video: bool
                ):
        
        
        self.neon_timeseries_dir = neon_timeseries_dir
        self.alternative_vid_path = alternative_vid_path
        self.output_dir = output_dir
        self.mapping_choice = mapping_choice
        self.optic_flow_choice = optic_flow_algorithm
        self.matcher = image_matcher
        try:
            self.refresh_time_thrshld = refresh_time_threshold_sec
        except NameError:
            self.refresh_time_thrshld = None
        try:
            self.render_comparison_video = render_video_comparison
        except NameError:
            self.render_comparison_video = False
        try:
            self.render_video = render_video
        except NameError:
            self.render_video = False


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

        # Let's load in the first chunk. This is simple
        start_chunk_path = files[0]
        start_chunk: np.ndarray = np.load(start_chunk_path) 

        # If there is a single file for the recording, this is simple.
        # The last chunk cannot be empty 
        end_chunk: np.ndarray | None = None
        if(len(files) == 1):
            end_chunk_path: str = files[-1]
            end_chunk = np.load(end_chunk_path)

        # Otherwise, let's check the last two chunks 
        else:
            for chunk_num in range(-2, 0):
                end_chunk_path: str = files[chunk_num]
                temp_chunk = np.load(end_chunk_path)
                if(len(temp_chunk) != 0):
                    end_chunk = temp_chunk

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

def convert_avi_to_mp4(input_path: str, output_path: str) -> None:
    cmd = [
        "ffmpeg",
        "-y",                     # force overwrite
        "-i", input_path,
        "-c:v", "libx264",
        "-crf", "18",
        "-preset", "slow",
        "-c:a", "aac",
        "-b:a", "192k",
        output_path
    ]

    subprocess.run(cmd, check=True)


"""
Wrapper function for the Pupil Labs Egocentric Video Mapper 
so that one does not need to fill out and run the jupyter notebook 
every time
"""
def run_egocentric_video_mapper(neon_timeseries_dir: str, 
                                alternative_vid_path: str,  
                                output_dir: str, 
                                mapping_choice: Literal["Fixations", "Gaze", "Both"]='Both',
                                refresh_time_threshold_sec: float=0.5 ,
                                render_video: bool= False,
                                render_video_comparison: bool= False,
                                optic_flow_algorithm: Literal["Lucas-Kanade", "Gunnar Farneback"] = "Lucas-Kanade", 
                                image_matcher: Literal["Efficient_LOFTR", "LOFTR_indoor"] = "Efficient_LOFTR",
                                show_video_preview: bool=False
                               ) -> None:

    # Show a video preview of the two videos if desired 
    if(show_video_preview is True):
        show_videos_preview(neon_timeseries_dir, alternative_vid_path)

    # First, we need to convert our video to .mp4 if it is not already in 
    # that format due to the pupil labs code requiring a .avi video 
    alternate_vid_dir: str = os.path.dirname(alternative_vid_path)
    temp_output_path: str = os.path.join(alternate_vid_dir, "temp.mp4")
    convert_avi_to_mp4(alternative_vid_path, temp_output_path)

    # Package the args for the egocentric video mapper 
    args: Args = Args(neon_timeseries_dir, 
                      temp_output_path, 
                      output_dir, 
                      mapping_choice, 
                      optic_flow_algorithm, 
                      image_matcher,
                      refresh_time_threshold_sec, 
                      render_video_comparison, 
                      render_video
                    )

    # Run the egocentric video mapper
    egocentric_main.main(args)

    # Remove the temporary video 
    os.remove(temp_output_path)

    return 


def main():
    pass 

if(__name__ == "__main__"):
    pass 