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
import argparse
import scipy.io
import tqdm 
import h5py


# Import relevant custom libraries with helper functions and constants 
light_logger_analysis_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLoggerAnalysis")
video_io_util_path: str = os.path.join(light_logger_analysis_dir_path, "code", "library", "matlabIO", "python_libraries")
world_util: str = os.path.join(light_logger_analysis_dir_path, "code", "library", "sensor_utility")
for path in (light_logger_analysis_dir_path, video_io_util_path, world_util):
    assert os.path.exists(path), f"Expected path: {path} does not exist"
    sys.path.append(path)

import video_io 
import world_util

def parse_args() -> str:
    # Initialize argument parser 
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Extract gaze targets from world gaze calibration video")

    # Add the path to video argument 
    parser.add_argument("path_to_video", type=str, help="Path to world camera video from gaze calibration")
    parser.add_argument("path_to_output_file", type=str, help="Path to numpy output file")
    parser.add_argument("--intended_gaze_targets_path", "--ig", default=None, required=False, type=str, help="GAZE CAL ONLY | Path to the gaze targets' positions in degrees when they were displayed on screen")
    parser.add_argument("--frame_indices", "--fi", nargs="*", default=[], type=int, help="APRIL TAG ONLY | List of frames to use for selection")

    # Parse the args
    args: object = parser.parse_args()

    return args.path_to_video,  args.path_to_output_file, args.intended_gaze_targets_path, args.frame_indices

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
                         peak_min_distance: int=240,
                         verbose: bool=False
                        ) -> tuple[int | object]:
    # If np.ndarray, not yet implemeneted 
    if(isinstance(video, np.ndarray)):
        raise NotImplementedError
    
    # Define the RGB mask variable. We will use this to calculate 
    # just the red pixels. It will be calculated based on the size 
    # of the images in the video 
    RGB_mask: np.ndarray = world_util.generate_RGB_mask(video_io.inspect_video_framesize(video))
    red_pixel_locations: tuple[np.ndarray] = np.where(RGB_mask == "R")
    red_pixel_rows: np.ndarray = red_pixel_locations[0]
    red_pixel_cols: np.ndarray = red_pixel_locations[1] 

    # Find the number of frames in the video 
    num_video_frames: int = video_io.inspect_video_frame_count(video)

    # Open a reader to the video 
    frame_reader: cv2.VideoCapture = cv2.VideoCapture(video)

    # Parse frames from the video, finding their average red 
    # channel and saving their frame num + average red channel
    red_intensities: np.ndarray = np.zeros((num_video_frames, 2), dtype=np.float64)
    frame_iterator: Iterable = range(num_video_frames) if verbose is False else tqdm.tqdm(range(num_video_frames), desc="Finding stimulus start/end")
    for frame_num in frame_iterator:
        # Read a frame 
        ret, frame = frame_reader.read()

        # If no frame arrived, then we are done 
        if(not ret):
            break

        # Convert the frame to grayscale 
        frame: np.ndarray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        # Otherwise, we have a frame, 
        # Let's find its red intensity and 
        red_intensities[frame_num, :] = [frame_num, np.mean(frame[red_pixel_rows, red_pixel_cols])]

    # Cleanly close the frame reader
    frame_reader.release()

    # Find the peaks
    peaks, props = find_peaks(
        red_intensities[:, 1],
        height=peak_height,
        prominence=peak_prominence,
        distance=peak_min_distance
    )

    # Get frame numbers of the peak indices
    peak_frames = red_intensities[peaks, 0]

    # Add delta to frame numbers (NOT to indices)
    delta = 200
    start, end = [ int(item) for item in (peak_frames + delta)[:2]]

    # Visualization
    if visualize_results is True:
        fig, ax = plt.subplots()

        # Main intensity curve
        ax.plot(
            red_intensities[:, 0],
            red_intensities[:, 1],
            '-x',
            label="Red Intensity by Frame"
        )

        # Mark detected peaks
        ax.plot(
            red_intensities[peaks, 0],
            red_intensities[peaks, 1],
            'x',
            label="Peaks"
        )

        # Mark start/end (shifted frames)
        ax.plot(
            peak_frames[:2] + delta,
            red_intensities[peaks[:2], 1],
            'x',
            label="Start/End"
        )

        ax.set_title("Stimulus Start/End")
        ax.set_xlabel("Frame #")
        ax.set_ylabel("Intensity")
        ax.legend()

        plt.show()

    return start, end, fig


"""Find the background image of a video.
   Calculated as the average from of all the videos
"""
def find_background_image(video: str | np.ndarray,
                          start_frame: int, end_frame: int | None=None,
                          is_grayscale: bool=False,
                          visualize_results: bool=False,
                          verbose: bool=False
                         ) -> np.ndarray | tuple[np.ndarray, object]:

    if(isinstance(video, np.ndarray)):
        NotImplementedError("Array support not yet added")

    # Allocate the sum frame 
    frame_size: tuple[int] = video_io.inspect_video_framesize(video)
    frame_sum: np.ndarray = np.zeros(frame_size if is_grayscale is True else tuple(list(frame_size)+[3]), dtype=np.float64)

    # Open a frame reader to the video at the start frame 
    frame_reader: cv2.VideoCapture = cv2.VideoCapture(video)
    frame_reader.set(cv2.CAP_PROP_POS_FRAMES, start_frame)

    # Parse frames from the video, gathering their features
    frame_iterator: Iterable = range(start_frame, end_frame) if verbose is False else tqdm.tqdm(range(start_frame, end_frame), desc="Finding background image")
    for frame_num in frame_iterator:
        # Read a frame from the video 
        ret, frame = frame_reader.read()

        # Break if no frame returned 
        if(not ret):
            break 

        # If frame is an empty frame (such as when chunks)
        # are being written, discard this 
        if( (frame == 0).all()):
            continue
        
        # If grayscale, covert frame to gray 
        if(is_grayscale is True):
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            
        # Otherwise, we have a frame, 
        # so we let's add it to the running average
        frame_sum += frame

        # Increment the frame number 
        frame_num += 1

    # Cleanly close the video reader 
    frame_reader.release() 

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
def auto_extract_target_circles(video: str, 
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
    read_process: object = mp.Process(target=video_io.destruct_video, 
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
    visualization_fps: int = video_io.inspect_video_FPS(video)
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

"""Manually extract"""
def manual_extract_target_circles(background_img: np.ndarray, path_to_intended_targets: str | None=None) -> np.ndarray:    
    
    if(path_to_intended_targets is not None):
        # First, let's load in the intended order of the points 
        intended_positions: np.ndarray | None = None
        with h5py.File(path_to_intended_targets, "r") as f:
            intended_positions = np.transpose(np.array(f["taskData"]["gaze_target_positions_deg"][()]))

        # FIGURE 1 — Intended order
        fig1, ax1 = plt.subplots()
        ax1.set_title("Gaze targets were presented in the following order")
        for i, (x, y) in enumerate(intended_positions, start=1):
            ax1.plot(x, y, 'o')
            ax1.text(x, y, str(i), ha='left', va='top')

        plt.show(block=False)  

  
    fig2, ax2 = plt.subplots()
    ax2.imshow(background_img, cmap=None if background_img.ndim == 3 else "gray")
    ax2.set_title("Pick gaze target points in order [Enter to quit]")

    # Click until Enter is pressed
    points: list[tuple] = fig2.ginput(0, timeout=0)
    plt.close(fig2)

    # FIGURE 3 — Verification
    fig3, ax3 = plt.subplots()
    ax3.imshow(background_img)
    ax3.set_title("Verify order gaze target points were clicked")

    for i, (x, y) in enumerate(points, start=1):
        ax3.plot(x, y, 'o')
        ax3.text(x + 5, y - 5, str(i), ha='left', va='top')

    plt.show()

    return np.array(points)


def extract_gazecal_stimulus(video_path: str, intended_gaze_targets_path: str) -> np.ndarray:
    print("---Finding stimulus start/end---") 
    start, end, fig = find_stimulus_period(video_path, visualize_results=True, verbose=True)

    # Calculate the background image 
    print(f"---Finding background image of frames: [{start}, {end})---")
    background_gray: np.ndarray = find_background_image(video_path, 
                                                        is_grayscale=True,
                                                        start_frame=start,
                                                        end_frame=end,
                                                        visualize_results=False,
                                                        verbose=True
                                                       )
    
    # Manually select gaze target poinst 
    print("---Selecting targets---")
    gaze_target_points: np.ndarray = manual_extract_target_circles(background_gray, intended_gaze_targets_path)

    return gaze_target_points

def extract_april_tag_stimulus(source_video: str,
                               target_frames_idx: Iterable,
                              ) -> np.ndarray:

    # First, extract the target frames from the video 
    target_frames: np.ndarray = video_io.extract_frames_from_video(source_video, target_frames_idx, is_grayscale=True)

    # Construct the avg frame to display the targets in a ghosting format 
    background_img: np.ndarray = np.mean(target_frames, axis=(0,))

    # Retrieve the target points in screen space by clicking the image
    target_points: np.ndarray = manual_extract_target_circles(background_img)

    return target_points

def main() -> None:
    # Parse arguments
    print("---Parsing Args---")
    video_path, path_to_output, intended_gaze_targets_path, frame_indices = parse_args()
    assert os.path.exists(video_path), f"Video path: {video_path} does not exist"

    # Determine which type of stimulus we are selecting 
    gaze_target_points: np.ndarray | None = None
    if(intended_gaze_targets_path is None):
        assert len(frame_indices) > 0, "No frame indices passed"
        gaze_target_points = extract_april_tag_stimulus(video_path, frame_indices)
    else:
        assert intended_gaze_targets_path is not None, "No intended gaze targets path passed"
        assert os.path.exists(intended_gaze_targets_path), f"Gaze targets path: {intended_gaze_targets_path} does not exist"
        gaze_target_points = extract_gazecal_stimulus(video_path, intended_gaze_targets_path)
    
    print("---Saving selections---")
    assert gaze_target_points is not None, "Gaze target points is None"
    if(path_to_output.strip().endswith(".mat")):
        scipy.io.savemat(path_to_output, {"gaze_targets": gaze_target_points})
    else:
        np.save(path_to_output, gaze_target_points)


if(__name__ == "__main__"):
    main() 