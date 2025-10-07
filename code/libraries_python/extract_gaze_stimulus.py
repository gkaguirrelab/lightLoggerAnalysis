import numpy as np
import multiprocessing as mp
import queue
from skimage import color, filters, feature, transform, util
import matplotlib.pyplot as plt
import os
import sys
from skimage.feature import canny
from skimage.transform import hough_circle, hough_circle_peaks
from typing import Iterable
import math
import heapq
import matplotlib.patches as patches
from scipy.signal import find_peaks

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

"""Given a dict of detected gaze targets and an image to display 
   on, visualize the results 
"""
def visualize_targets(background: np.ndarray, targets: dict[int, np.ndarray]) -> object:
    # Initialize a figure 
    fig, ax = plt.subplots() 

    # Imshow the background 
    ax.imshow(background, cmap="gray")

    # Draw the targets onto the background 
    for target_num, (cx, cy, r) in targets.items():
        circle = plt.Circle((cx, cy), r, fill=False, color='red', linewidth=2, label=f"Target: {target_num}")
        ax.add_patch(circle)

    # Ensure the figure is pretty 
    ax.set_title("Gaze Calibration Targets")
    ax.legend()

    # Display the figure 
    plt.show()

    return fig

"""Calculate the euclidean distance between circles in pixel space"""
def calculate_euclidean_distance(a: tuple[int], b: tuple[int]):
    return math.sqrt( (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)

"""Find the start and end of the stimulus period
   Denoted by the screen flashing a single frame of RED 
"""
def find_stimulus_period(video: str | np.ndarray,
                         visualize_results: bool=False, 
                         peak_height: float=100,
                         peak_prominence: int=5,
                         peak_min_distance: int=240
                        ) -> tuple[int]:
    # If np.ndarray, not yet implemeneted 
    if(isinstance(video, np.ndarray)):
        raise NotImplementedError
    
    # Otherwise, we will go through the frames 
    # of the video and find the peaks in average red channel 
    # intensity of each frame
     # Stream frames in and add them to the per pixel sum 
    frame_queue: mp.Queue = mp.Queue(maxsize=5)
    stop_event: object = mp.Event()
    frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                              args=(video, 0, float("inf"), 
                                                    False, frame_queue, stop_event)
                                             )
    # Start the subprocess 
    frame_stream_process.start()

    # Parse frames from the video, finding their average red 
    # channel and saving their frame num + average red channel
    red_intensities: list[int] = []
    frame_num: int = 0 
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
        # Let's find its red intensity and 
        # save it to the list # Note (cv2 reads in as BGR)
        red_intensities.append( [frame_num, np.mean(frame[:, :, 2])] )

        # Increment the frame number 
        frame_num += 1

    # Join the subprocess 
    frame_stream_process.join()     

    # Convert red intensities to np.ndarray 
    red_intensities = np.array(red_intensities)

    # Find the positive peaks in the signal corresponding 
    # to the start/end
    peaks, props = find_peaks(red_intensities[:, 1],
                              height=peak_height,                
                              prominence=peak_prominence,              
                              distance=peak_min_distance             
                             )
    
    # Find the approximate start/end (add a small delta)
    # to get out of the peak positions
    delta: int = 200 
    peaks += delta
    start, end = red_intensities[peaks, 0][:2] 

    # Visualize results if desired 
    if(visualize_results is True):
        fig, ax = plt.subplots() 
        ax.plot(red_intensities[:, 0], red_intensities[:, 1], '-x', label="Red Intensity by Frame")
        ax.plot(peaks, red_intensities[peaks][:, 1], 'x', label="Peaks")
        ax.plot(peaks[:2], red_intensities[peaks[:2], 1], 'x', label="Start/End")
        ax.set_title("Stimulus Start/End")
        ax.set_xlabel("Frame #")
        ax.set_ylabel("Intensity")
        ax.legend() 
        
        plt.show()

        return start, end, fig

    return start, end


"""Find the background image of a video.
   Calculated as the average from of all the videos
"""
def find_background_image(video: str | np.ndarray,
                          start_frame: int, end_frame: int | None=None,
                          is_grayscale: bool=False,
                          visualize_results: bool=False, 
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
                                              args=(video, start_frame, float("inf") if end_frame is None else end_frame, 
                                                    is_grayscale, frame_queue, stop_event)
                                             )
    # Start the subprocess 
    frame_stream_process.start()

    # Parse frames from the video, gathering their features
    frame_num: int = 0
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

        # If frame is an empty frame (such as when chunks)
        # are being written, discard this 
        if( (frame == 0).all()):
            continue
        
        # Otherwise, we have a frame, 
        # so we let's add it to the running average
        frame_sum += frame

        # Increment the frame number 
        frame_num += 1

    # Join the subprocess 
    frame_stream_process.join()   

    # Calculate the backgroudn image 
    background_img: np.ndarray = ( np.clip(frame_sum / frame_num, 0, 255) ).astype(np.uint8) 

    # Visualize the results if desired 
    if(visualize_results is True):
        fig, ax = plt.subplots() 
        ax.imshow(background_img, cmap='gray')
        ax.set_title("Background")
        plt.show()

        return background_img, fig
    
    return background_img 

"""From a given frame, threshold the frame and return 
   the circles from the frame
"""
def extract_target_circles(video: str, 
                           background_img: np.ndarray, 
                           start_frame: int=0, end_frame: int | None=None,
                           is_grayscale: bool=False,
                           threshold_value: int=127,
                           radius_range: Iterable=range(4, 7, 1),
                           min_intercircle_distance: float=8,
                           visualize_results: bool=False 
                          ) -> np.ndarray:
    # If given a np.ndarray, simply take the average frame 
    if(isinstance(video, np.ndarray)):
        raise NotImplementedError
    
    # Otherwise, if given a video, need to stream frames 
    # and take the average later 

    # Stream frames in and add them to the per pixel sum 
    frame_queue: mp.Queue = mp.Queue(maxsize=5)
    stop_event: object = mp.Event()
    frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                              args=(video, start_frame, float("inf") if end_frame is None else end_frame, 
                                                    is_grayscale, frame_queue, stop_event)
                                             )
    # Start the subprocess 
    frame_stream_process.start()

    # Define the circles array 
    circles: dict[tuple | None] = {}

    # Parse frames from the video, gathering their features
    num_detected_circles: int = 0
    frame_num: int = 0
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

        # Subtract the background from the frame
        background_subtracted: np.ndarray = np.clip(frame.astype(np.float64) - background_img.astype(np.float64), 0, 255).astype(np.uint8)

        # Then, threshold to just leave the circle remaining 
        thresholded: np.ndarray = background_subtracted > threshold_value

        if(frame_num % 120 == 0 and frame_num >= 120):
            plt.imshow(thresholded)
            plt.title(f"Frame Num: {frame_num}")
            plt.show()


        # Increment the frame num 
        frame_num += 1

        # Skip NULL frames 
        if(thresholded == 0).all():
            continue

        edges: np.ndarray = canny(thresholded, sigma=2.0)
        hough_res: object = hough_circle(edges, range(2, 5, 1))
        accums, cxs, cys, radii_found = hough_circle_peaks(hough_res, radius_range, total_num_peaks=1)

        # Compact the circles together for easy plotting in the future 
        frame_circles: list[tuple] =  [ tuple([int(cx), int(cy), int(radius)])
                                        for cx, cy, radius in zip(cxs, cys, radii_found)
                                      ] 
        
        # Skip frames that have no circles 
        if(len(frame_circles) == 0):
            continue        
        most_prominent_circle: tuple = frame_circles[0]  

        # If we have detected no circles before  
        if(len(circles) == 0):
            circles[num_detected_circles] = np.array(most_prominent_circle) 
            num_detected_circles += 1 
            continue
        
        plt.imshow(thresholded)
        plt.show()


        # If this circle has basically been detected before, 
        # then we skip 
        previously_detected_circles: list = list(circles.items())

        # First, check if this circle has not been seen before 
        if(all( calculate_euclidean_distance(previously_detected_circle, most_prominent_circle) > min_intercircle_distance for circle_num, previously_detected_circle in previously_detected_circles)):
            circles[num_detected_circles] = np.array(most_prominent_circle) 
            num_detected_circles += 1 
            continue

        # Otherwise, average the ones its close to 
        for circle_num, previously_detected_circle in previously_detected_circles:
            # Determine if the previously detected circle is approximately 
            # equal to the current circle 

            # If they are different circles, simply just add to the dict 
            if(calculate_euclidean_distance(previously_detected_circle, most_prominent_circle) < min_intercircle_distance):
                # Otherwise, average them 
                circles[circle_num] = np.array([ int(v) for v in (previously_detected_circle + np.array(most_prominent_circle)) / 2 ])

    # Join the subprocess 
    frame_stream_process.join()   

    # Visualize the results if desired
    fig: object | None = None
    if(visualize_results is True):
        fig = visualize_targets(background_img, circles)
        return circles, fig

    return circles 

"""Extract the position of gaze calibration targets 
   from a video per frame 
"""
def extract_gaze_stimulus(video: str | np.ndarray, 
                          start_frame: int=0, end_frame: int | None=None,
                          is_grayscale: bool=False,
                          threshold_value: int=127,
                          radius_range: Iterable=range(4, 7, 1),
                          min_intercircle_distance: float=8,
                          visualize_results: bool=False 
                         ) -> dict:
    
    # First find the background image. We will use this to background subtract from 
    # the entire video 
    background = find_background_image(video, 
                                                        start_frame, 
                                                        end_frame
                                                       )
    
    # Next extract the unique circles from the video
    # (location averaged for nearby circles)
    circles: dict = extract_target_circles(video,
                                           background,
                                           start_frame,
                                           end_frame,
                                           is_grayscale,
                                           threshold_value,
                                           radius_range,
                                           min_intercircle_distance,
                                           visualize_results)


    return circles 

def main() -> None:
    pass 

if(__name__ == "__main__"):
    main() 