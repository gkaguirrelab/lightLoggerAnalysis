import numpy as np
import multiprocessing as mp
import queue
import os
import sys

# Import relevant custom libraries with helper functions and constants 
light_logger_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLogger")
Pi_util_dir_path: str = os.path.join(light_logger_dir_path, "raspberry_pi_firmware", "utility")
pupil_util_path: str = os.path.join(light_logger_dir_path, "pupil")
for path in (Pi_util_dir_path, pupil_util_path):
    sys.path.append(path)

# Import Pi_util for helper functions 
import Pi_util 

# Import pupil_util for constants 
import pupil_util

"""Find the background image of a video.
   Calculated as the average from of all the videos
"""
def find_background_image(video: str | np.ndarray,
                          start_frame: int, end_frame: int | None=None,
                          is_grayscale: bool=False
                         ) -> np.ndarray:
    # If given a np.ndarray, simply take the average frame 
    if(isinstance(video, np.ndarray)):
        return np.average(video[start_frame:end_frame], axis=(0,))
    
    # Otherwise, if given a video, need to stream frames 
    # and take the average later 

    # Allocate the sum frame 
    frame_size: tuple[int] = Pi_util.inspect_video_framesize(video)
    frame_sum: np.ndarray = np.zeros(frame_size, dtype=np.float64)

    # Stream frames in and add them to the per pixel sum 
    frame_queue: mp.Queue = mp.Queue(maxsize=5)
    stop_event: object = mp.Event()
    frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                              args=(video, 0, float("inf") if end_frame is None else end_frame, 
                                                    is_grayscale, frame_queue, stop_event)
                                             )
    # Start the subprocess 
    frame_stream_process.start()

    # Parse frames from the video, gathering their features
    frame_num: int = 1
    while(True):
        # Attempt to retrieve a frame from the frame queue 
        try:
            frame: np.ndarray | None = frame_queue.get(timeout=5)
        except queue.Empty:
            # Stop the subprocess
            stop_event.set()
            frame_stream_process.join()   
            raise Exception("Did not recieve frame in time")

        # If no frame arrived, then we are done 
        if(frame is None):
            break
        
        # Otherwise, we have a frame, 
        # so we let's add it to the running average
        frame_sum += frame

        # Increment the frame number 
        frame_num += 1

    # Join the subprocess 
    frame_stream_process.join()   

    return ( np.clip(frame_sum / frame_num, 0, 255) ).astype(np.uint8) 


"""Extract the position of gaze calibration targets 
   from a video per frame 
"""
def extract_gaze_stimulus(video: str | np.ndarray, 
                          start_frame: int=0, end_frame: int | None=None
                         ) -> np.ndarray:
    
    background: np.ndarray = find_background_image(video, 
                                                   start_frame, 
                                                   end_frame
                                                  )


    return background

    return 

def main() -> None:
    pass 

if(__name__ == "__main__"):
    main() 