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
import matplotlib.patches as patches
from scipy.signal import find_peaks
from typing import Literal
import cv2
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import shutil

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
                        ) -> tuple[int | object]:
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
                         ) -> np.ndarray | tuple[np.ndarray, object]:

    # Allocate the sum frame 
    frame_size: tuple[int] = Pi_util.inspect_video_framesize(video)
    frame_sum: np.ndarray = np.zeros(frame_size if is_grayscale is True else tuple(list(frame_size)+[3]), dtype=np.float64)

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
        ax.imshow(background_img, cmap='gray' if is_grayscale is True else None)
        ax.set_title("Background")
        plt.show()

        return background_img, fig
    
    return background_img 

"""Given a matplotlib figure, render the figure as an image
   figure and return it as a numpy array
"""
def rasterize_figure(matplotlib_figure: object) -> np.ndarray:
    # Rasterize the image
    canvas: object = FigureCanvas(matplotlib_figure)
    canvas.draw()
    
    # Retrieve the rasterized image
    width, height = matplotlib_figure.canvas.get_width_height()
    visualized_rgb: np.ndarray = np.frombuffer(matplotlib_figure.canvas.tostring_rgb(), dtype=np.uint8)
    visualized_rgb: np.ndarray = visualized_rgb.reshape((height, width, 3))

    # If you’re writing video with OpenCV, convert to BGR
    visualized_bgr = cv2.cvtColor(visualized_rgb, cv2.COLOR_RGB2BGR).astype(np.uint8)

    return visualized_bgr

"""From a given frame, threshold the frame and return 
   the circles from the frame
"""
def extract_target_circles(video: str, 
                           background_img: np.ndarray, 
                           start_frame: int=0, end_frame: int | None=None,
                           is_grayscale: bool=False,
                           threshold_value: int=127,
                           radius_range: Iterable=range(3, 7, 1),
                           min_intercircle_distance: float=8,
                           visualize_results: Literal["None", "Circles", "Video"]="None", 
                           visualization_output_path: str | None=None
                          ) -> list[tuple] | tuple[list, object]:
  

    # Subfunction to init the output stream 
    def _init_output_stream(output_path: str, output_fps: float, output_frame_size_w_h: tuple[int]) -> cv2.VideoWriter:
        # Initialize the video writer
        output_stream: cv2.VideoWriter = cv2.VideoWriter(output_path,
                                                         0,
                                                         output_fps,
                                                         output_frame_size_w_h,
                                                         isColor=True
                                                        )
        
        return output_stream

    
    # Stream frames in 
    read_queue: mp.Queue = mp.Queue(maxsize=5)
    read_stop: object = mp.Event()
    read_process: object = mp.Process(target=Pi_util.destruct_video, 
                                              args=(video, start_frame, float("inf") if end_frame is None else end_frame, 
                                                    is_grayscale, read_queue, read_stop)
                                             )

    # Start the subprocesses 
    read_process.start()

    # Define the circles array 
    circles: dict[tuple | None] = {}

    # Initialize the video visualization writer variables. This will be 
    # assigned a value later if we want to actually use video visualization
    video_writer: cv2.VideoWriter | None = None 
    visualization_frame_size: tuple[int] | None = None 
    visualization_fps: int = Pi_util.inspect_video_FPS(video)
    visualized_bgr: np.ndarray | None = None

    # Parse frames from the video, gathering their features
    num_detected_circles: int = 0
    frame_num: int = 0
    while(True):
        # Attempt to retrieve a frame from the frame queue 
        try:
            frame: np.ndarray | None = read_queue.get(timeout=5)
        except queue.Empty:
            # Stop the subprocess
            read_stop.set()
            read_process.join()   
            raise Exception("Did not recieve frame in time")

        # If no frame arrived, then we are done 
        if(frame is None):
            break

        # Subtract the background from the frame
        background_subtracted: np.ndarray = np.clip(frame.astype(np.float64) - background_img.astype(np.float64), 0, 255).astype(np.uint8)

        # Then, threshold to just leave the circle remaining 
        thresholded: np.ndarray = background_subtracted > threshold_value

        # Initialize the plots for visualization 
        fig: object | None = None 
        axes: np.ndarray | None = None
        if(visualize_results == "Video"):
            fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=300)
            axes = axes.flatten() 
            
            # Generate supra title for the figure
            fig.suptitle(f"Frame: {frame_num}", fontsize=12, y=0.8)

            # Display the original frame 
            axes[0].set_title(f"Orignal Frame")
            axes[0].imshow(frame, cmap='gray')

            # Display the background subtracted image and thresholded 
            axes[1].set_title(f"Background Subtracted + Thresholded")
            axes[1].imshow(background_subtracted, cmap='gray')
            axes[1].imshow(thresholded, cmap='Reds', alpha=0.25)


        # Increment the frame num 
        frame_num += 1

        # Skip NULL frames 
        if(thresholded == 0).all():
            # write to the video now if no circle detected
            if(visualize_results == "Video"):
                # Rasterize the figure image 
                visualized_bgr: np.ndarray = rasterize_figure(fig)
                height, width = visualized_bgr.shape[:2]

                # Initialize the video writer if not already initialized 
                if(video_writer is None):
                    video_writer = _init_output_stream(visualization_output_path if visualization_output_path is not None else "./visualized_targets.avi",
                                                       visualization_fps,
                                                       (width, height)
                                                      )
                    visualization_frame_size = visualized_bgr.shape 

                # Assert the frame size is equal to what we expect 
                if(visualization_frame_size != visualized_bgr.shape):
                    print("BIG PROBLEM!!!!", flush=True)
                    raise RuntimeError(f"Change in shape of visualized output. Expected {visualization_frame_size} got {visualized_bgr.shape}")

                # Write the frame to the video
                video_writer.write(visualized_bgr)


                # Close the figure now that we are done with it
                plt.close(fig)

            continue



        edges: np.ndarray = canny(thresholded, sigma=1.5)
        hough_res: object = hough_circle(edges, radius_range)
        accums, cxs, cys, radii_found = hough_circle_peaks(hough_res, radius_range, total_num_peaks=1)

        # Compact the circles together for easy plotting in the future 
        frame_circles: list[tuple] =  [ tuple([int(cx), int(cy), int(radius)])
                                        for cx, cy, radius in zip(cxs, cys, radii_found)
                                      ] 
        
        # Skip frames that have no circles 
        if(len(frame_circles) == 0):
            # write to the video now if no circle detected
            if(visualize_results == "Video"):
                # Rasterize the figure image 
                visualized_bgr: np.ndarray = rasterize_figure(fig)
                height, width = visualized_bgr.shape[:2]

                # Initialize the video writer if not already initialized 
                if(video_writer is None):
                    video_writer = _init_output_stream(visualization_output_path if visualization_output_path is not None else "./visualized_targets.avi",
                                                       visualization_fps,
                                                       (width, height)
                                                      )
                    visualization_frame_size = visualized_bgr.shape 

                # Assert the frame size is equal to what we expect 
                if(visualization_frame_size != visualized_bgr.shape):
                    print("BIG PROBLEM!!!!", flush=True)
                    raise RuntimeError(f"Change in shape of visualized output. Expected {visualization_frame_size} got {visualized_bgr.shape}")


                # write the frame to the video 
                video_writer.write(visualized_bgr)

                # Close the figure now that we are done with it
                plt.close(fig)

            continue        
        
        most_prominent_circle: tuple = frame_circles[0] 

        # Apply most prominent circle to the visualization 
        # video if desired 
        if(visualize_results == "Video"):
            circle_patch = patches.Circle(most_prominent_circle[:2], most_prominent_circle[2], 
                                          edgecolor='green', 
                                          facecolor='none', 
                                          linewidth=1,
                                          fill=False,
                                          label="Detected Circle"
                                        )
            axes[1].add_patch(circle_patch)
            axes[1].legend() 

            # Rasterize the figure image 
            visualized_bgr: np.ndarray = rasterize_figure(fig)
            height, width = visualized_bgr.shape[:2]

            # Initialize the video writer if not already initialized 
            if(video_writer is None):
                video_writer = _init_output_stream(visualization_output_path if visualization_output_path is not None else "./visualized_targets.avi",
                                                    visualization_fps,
                                                    (width, height)
                                                    )
                visualization_frame_size = visualized_bgr.shape 

            # Assert the frame size is equal to what we expect 
            if(visualization_frame_size != visualized_bgr.shape):
                print("BIG PROBLEM!!!!", flush=True)
                raise RuntimeError(f"Change in shape of visualized output. Expected {visualization_frame_size} got {visualized_bgr.shape}")

            # Write the frame to the video
            video_writer.write(visualized_bgr)

            # Close the figure now that we are done with it
            plt.close(fig)

        # If we have detected no circles before  
        if(len(circles) == 0):
            circles[num_detected_circles] = np.array(most_prominent_circle) 
            num_detected_circles += 1 
            continue

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
    read_process.join()   

    # Close video visualizaiton if desired
    # and remove the tempdir 
    if(visualize_results == "Video"):
        video_writer.release()  

    # Visualized only the detected circles if desired
    if(visualize_results == "Circles"):
        fig: object | None = None
        fig = visualize_targets(background_img, circles)

        # If output path given, save the figure 
        if(visualization_output_path is not None):
            plt.savefig(visualization_output_path)

        return circles.items(), fig
    
    return circles.items()  

def main() -> None:
    pass 

if(__name__ == "__main__"):
    main() 