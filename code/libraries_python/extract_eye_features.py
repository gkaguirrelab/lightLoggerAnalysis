import sys 
import os
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pupil_detectors import Detector2D
from pye3d.detector_3d import CameraModel, Detector3D, DetectorMode
from typing import Literal
import queue
import multiprocessing as mp
from tqdm import tqdm

# Import relevant custom libraries with helper functions and constants 
light_logger_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLogger")
Pi_util_dir_path: str = os.path.join(light_logger_dir_path, "raspberry_pi_firmware", "utility")
pupil_util_path: str = os.path.join(light_logger_dir_path, "pupil")
for path in (Pi_util_dir_path, pupil_util_path):
    sys.path.append(path)

sys.path.append(os.path.join(os.path.dirname(__file__), "pylids"))
import pylids

# Import Pi_util for helper functions 
import Pi_util 

# Import pupil_util for constants 
import pupil_util

"""Helper function to analyze a video with pupil labs for pupil features
   given either a path to the video or a series of frames
"""
def pupil_labs_analyze_video(video: str | np.ndarray, 
                             is_grayscale: bool=False
                             ) -> list[dict]:
    # Create 2D detector from pupil labs 
    detector_2d: object = Detector2D()
    
    # Create pye3D detector from pupil labs
    camera: object = CameraModel(focal_length=pupil_util.PUPIL_FOCAL_LENGTH, resolution=pupil_util.PUPIL_FRAME_SHAPE.tolist())
    detector_3d: object = Detector3D(camera=camera, long_term_mode=DetectorMode.blocking)

    # Initialize return value 
    eye_features_by_frame: list[dict] = []
    
    # Iterate over the frames of the video 
    # when passed in as a numpy array 
    if(isinstance(video, np.ndarray)):
        # If the video is not grayscale, convert to grayscale 
        if(len(video.shape) > 3):
            video = np.array([cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
                              for frame in video
                             ], 
                             dtype=np.uint8
                            )

        for (frame_num, frame) in tqdm(enumerate(video), desc="Visualizing"):
            # Run a 2D detector on the video frame 
            result_2d = detector_2d.detect(frame)
            result_2d["timestamp"] = frame_num / pupil_util.PUPIL_CAM_FPS

            # Pass 2D detection result to 3D detector
            result_3d = detector_3d.update_and_detect(result_2d, frame)

            # Append the eye features for this frame to the growing list 
            eye_features_by_frame.append(result_3d)        

        return eye_features_by_frame
    
    # Otherwise, we must iterate over the frames of the video 
    # only having one in memory at a time, so that we do not blow up our memory 
    frame_queue: mp.Queue = mp.Queue(maxsize=5)

    # Initialize a second process to stream frames from the video 
    stop_event: object = mp.Event()
    frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                              args=(video, 0, float("inf"), is_grayscale, frame_queue, stop_event)
                                             )
    # Start the subprocess 
    frame_stream_process.start()

    # Inspect the video frame count to know how long to iterate for
    num_frames: int = Pi_util.inspect_video_frame_count(video)

    # Parse frames from the video, gathering their features
    for frame_num in tqdm(range(num_frames), desc="Analyzing"):
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
        
        # Otherwise, analyze the frame

        # Run a 2D detector on the video frame 
        result_2d = detector_2d.detect(frame)
        result_2d["timestamp"] = frame_num / pupil_util.PUPIL_CAM_FPS

        # Pass 2D detection result to 3D detector
        result_3d = detector_3d.update_and_detect(result_2d, frame)

        # Append the eye features for this frame to the growing list 
        eye_features_by_frame.append(result_3d) 

    # Join the subprocess 
    frame_stream_process.join()   

    return eye_features_by_frame

"""Helper function to analyze a video with pylids for pupil or eyelids
   given a path to the video
"""
def pylids_analyze_video(video: str, 
                         target: Literal["pupil", "eyelid"]
                        ) -> list[dict]:
    # Ensure the video argument is a string
    assert isinstance(video, str), "Using pylids method, video argument must be string path to video"
    
    # Perform pylids feature extraction
    pylids_out: dict[str, list] = pylids.analyze_video(eye_vid=video, 
                                                       estimate_eyelids=True if target=="eyelid" else False,  
                                                       model_name=f'pytorch-{target}-v1', 
                                                       save_vid=False
                                                      )

    # Next, convert format from dict[str, list[dict]] to a list of dictionaries
    pylid_features: list[dict] = []

    # First, let's find out how many frames we had 
    num_frames: int = len(pylids_out['dlc_confidence'])
    assert all(len(value) == num_frames for key, value in pylids_out.items()), "Pylids output shape mismatch. Unable to determine num frames"

    # Next, let's iterate over these frames and recompose into new dictionary 
    for frame_num in range(num_frames):
        # Initialize a frame dictionary for all the keys and values for this frame  
        frame_pylid_features: dict[str, object] = {key: value_arr[frame_num]
                                                   for key, value_arr in pylids_out.items()
                                                  } 
        
        # Append per frame eyelid features to the list 
        pylid_features.append(frame_pylid_features)

    return pylid_features

"""Plot detected pupil on a given frame"""
def visualize_pupil(frame: np.ndarray, 
                    frame_pupil_features: dict[str, object], 
                    ax: plt.Axes=None,
                    method: Literal["pupil-labs", "pylids"]="pupil-labs"
                   ) -> None | tuple:
    
    # Assert the frame is 2D grayscale
    assert(len(frame.shape) == 2)

    # Generate a figure if not given an axis 
    fig: object = None
    if(ax is None):
        fig, ax = plt.subplots()

    # Show the frame on the axis
    ax.imshow(frame, cmap='gray')

    # Plot the pupil on the frame
    if(method == "pylids"):
        # Gather ellipse parameters to plot the pupil
        center_x, center_y = frame_pupil_features['ellipse']['center']
        width, height = frame_pupil_features['ellipse']['axes']
        theta: float = frame_pupil_features['ellipse']['angle']

        # Generate pupil ellipse
        ellipse_patch: object = patches.Ellipse((center_x, center_y), width, height, fill =False, angle=theta, edgecolor='red')
        
        # Add ellipse to the image
        ax.add_patch(ellipse_patch)
    else:
        # Retrieve the ellipse parameters from the frame feature dict
        ellipse: dict = frame_pupil_features["ellipse"]

        # Retrieve the ellipse parameters from the ellipse dictionary 
        cx, cy = tuple(int(v) for v in ellipse["center"])
        major, minor = ellipse['axes']
        angle: float = ellipse['angle']

        # Generate the ellipse patch 
        ellipse_patch: object = patches.Ellipse((cx, cy),             # (x_center, y_center)
                                                width=major,          # full length of major axis
                                                height=minor,         # full length of minor axis
                                                angle=angle,          # degrees, counterclockwise
                                                edgecolor="red",      # outline color
                                                facecolor="none",     # no fill
                                                linewidth=2
                                               )
        # Add the patch to the image 
        ax.add_patch(ellipse_patch)

    # Pretty the plot 
    ax.set_title("Pupil")

    # Show the plot if we generated the ax here 
    if(fig is not None):
        plt.show()

    # Return the figure if generated from here
    return None if fig is None else (fig, ax)

"""Given a path to a video or frames of a video, extract
   pupil features of the eye per frame
"""
def extract_pupil_features(video: str | np.ndarray,
                           is_grayscale: bool=False, 
                           visualize_results: bool=False, 
                           method: Literal["pupil-labs", "pylids"]="pupil-labs"
                          ) -> list[dict]:

    # Initialize eye features variable 
    pupil_features: list[dict]

    # If desired method is pylids, then its quite easy, 
    # just call the simple analysis function 
    if(method == "pylids"):
        pupil_features = pylids_analyze_video(video, "pupil")
    # Otherwise, we want to use the pupil labs method,
    else:
        pupil_features = pupil_labs_analyze_video(video, is_grayscale)

    # Check to ensure that the video is well formed
    assert Pi_util.inspect_video_frame_count(video) == len(pupil_features), "Video framecount not equal to analyzed frames by pylids. Some frames may be corrupted"

    # If we want to visualize the results 
    if(visualize_results is True):
        # If the video is an np.ndarray, simply 
        # iterate over the frames and visualize
        if(isinstance(video, np.ndarray)):
            # If the video is not grayscale, convert to grayscale 
            if(len(video.shape) > 3):
                video = np.array([cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
                                for frame in video
                                ], 
                                dtype=np.uint8
                                )        
        
            # Iterate over the frames of the video and display the results 
            for (frame_num, (frame, frame_pupil_features)) in tqdm(enumerate( zip(video, pupil_features) ), desc="Visualizing"):
                # Generate a figure to display the results on 
                fig, axes = plt.subplots()
                fig.suptitle(f"Frame: {frame_num} | Eye Features")

                # Plot the pupil features 
                visualize_pupil(frame, 
                                frame_pupil_features, 
                                method=method,
                                ax=axes
                              )
                
                # Show the figure 
                plt.show()

                # Close the figure generated 
                plt.close(fig)

            return pupil_features     
    
        # Otherwise, we must iterate over the frames of the video 
        # only having one in memory at a time, so that we do not blow up our memory 
        frame_queue: mp.Queue = mp.Queue(maxsize=5)

        # Initialize a second process to stream frames from the video 
        stop_event: object = mp.Event()
        frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                                  args=(video, 0, float("inf"), is_grayscale, 
                                                        frame_queue, stop_event
                                                       )
                                                 )
        # Start the subprocess 
        frame_stream_process.start()

        # Parse frames from the video, gathering their features
        num_frames: int = Pi_util.inspect_video_frame_count(video)
        for frame_num in tqdm(range(num_frames), desc="Visualizing"):
            # Attempt to retrieve a frame from the frame queue
            try: 
                frame: np.ndarray | None = frame_queue.get(timeout=5)
            except queue.Empty:
                # Close the subprocess on error 
                stop_event.set() 
                frame_stream_process.join() 
                raise Exception("Did not receive frame in time")

            # If no frame arrived, then we are done 
            if(frame is None):
                break
            
            # Otherwise, visualize the frame and pupil features 
            # Generate a figure to display the results on 
            fig, axes = plt.subplots()
            fig.suptitle(f"Frame: {frame_num} | Eye Features")

            # Plot the pupil features 
            visualize_pupil(frame, 
                            pupil_features[frame_num], 
                            method=method,
                            ax=axes
                        )
            
            # Show the desired figure 
            plt.show()

            # Close the figure generated 
            plt.close(fig)

        # Join the subprocess 
        frame_stream_process.join()   

    return pupil_features

"""Plot detected eyelids on a given frame"""
def visualize_eyelids(frame: np.ndarray, 
                      frame_eyelid_features: dict[str, object], 
                      ax: plt.Axes=None
                     ) -> None:
    
    # Assert the frame is 2D grayscale
    assert(len(frame.shape) == 2)

    # Generate a figure if not given an axis 
    fig: object = None
    if(ax is None):
        fig, ax = plt.subplots()

    # Show the frame on the axis
    ax.imshow(frame, cmap='gray')

    # Plot the eyelids on the frame
    ax.plot(frame_eyelid_features['eyelid_x'], frame_eyelid_features['eyelid_lo_y'], 'r')
    ax.plot(frame_eyelid_features['eyelid_x'], frame_eyelid_features['eyelid_up_y'], 'r')

    # Pretty the plot 
    ax.set_title("Eyelids")

    # Show the plot if we generated the ax here 
    if(fig is not None):
        plt.show()

    # Return the figure if generated from here
    return None if fig is None else (fig, ax)

"""Given a path to a video or frames of a video, extract
   eyelid  features of the eye per frame
"""
def extract_eyelid_features(video: str | np.ndarray,
                            is_grayscale: bool=False,
                            visualize_results: bool=False, 
                           )-> list[dict]:
    # Extract eyelid features with pylids
    eyelid_features: dict[str, dict] = pylids_analyze_video(video, "eyelid")

    # Check to ensure that the video is well formed
    assert Pi_util.inspect_video_frame_count(video) == len(eyelid_features), "Video framecount not equal to analyzed frames by pylids. Some frames may be corrupted"

    # Visualize results if desired 
    if(visualize_results is True):
        # If passed in a numpy array of frames, 
        # simply iterate over them and visualize 
        if(isinstance(video, np.ndarray)):
            # If the video is not grayscale, convert to grayscale 
            if(len(video.shape) > 3):
                video = np.array([cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
                                for frame in video
                                ], 
                                dtype=np.uint8
                                )        
                

            # Iterate over the frames of the video and display the results 
            for (frame_num, (frame, frame_eyelid_features)) in tqdm(enumerate( zip(video, eyelid_features) ), desc="Visualizing"):
                # Generate a figure to display the results on 
                fig, axes = plt.subplots()
                fig.suptitle(f"Frame: {frame_num} | Eye Features")

                # Plot the pupil features 
                visualize_pupil(frame, 
                                frame_eyelid_features, 
                                ax=axes
                               )
                
                # Close the figure generated 
                plt.close(fig)

            return eyelid_features
        
        # Otherwise, we must iterate over the frames of the video 
        # only having one in memory at a time, so that we do not blow up our memory 
        frame_queue: mp.Queue = mp.Queue(maxsize=5)

        # Initialize a second process to stream frames from the video 
        stop_event: object = mp.Event()
        frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                                  args=(video, 0, float("inf"), is_grayscale, 
                                                        frame_queue, stop_event
                                                       )
                                                 )
        # Start the subprocess 
        frame_stream_process.start()

        # Parse frames from the video, gathering their features
        num_frames: int = Pi_util.inspect_video_frame_count(video)
        for frame_num in tqdm(range(num_frames), desc="Visualizing"):
            # Attempt to retrieve a frame from the frame queue 
            try:
                frame: np.ndarray | None = frame_queue.get(timeout=5)
            except queue.Empty:
                # Close the subprocess on error 
                stop_event.set()
                frame_stream_process.join()
                raise Exception("Did not receive frame in time")

            # If no frame arrived, then we are done 
            if(frame is None):
                break
            
            # Otherwise, visualize the frame and pupil features 
            # Generate a figure to display the results on 
            fig, axes = plt.subplots()
            fig.suptitle(f"Frame: {frame_num} | Eye Features")

            # Plot the pupil features 
            visualize_eyelids(frame, 
                              eyelid_features[frame_num], 
                              ax=axes
                             )
            
            # Show the desired figure 
            plt.show()

            # Close the figure generated 
            plt.close(fig)

        # Join the subprocess 
        frame_stream_process.join()   

    return eyelid_features

"""Given a path to a pupil video or frames of a video, 
   extract features of the eye 
   per frame from the video using either pylids or 
   pupil labs for pupil detection, and pylids for 
   eyelid detection. 
"""
def extract_eye_features(video: str | np.ndarray, 
                         is_grayscale: bool=True,
                         visualize_results: bool=False,
                         pupil_feature_method: Literal["pupil-labs", "pylids"]="pupil-labs"
                        ) -> list[dict]:
    
    # First, extract the pupil features of the video 
    pupil_features: list[dict] = extract_pupil_features(video, 
                                                        is_grayscale, # Do not visualize single features if we want all features  
                                                        not visualize_results if visualize_results is True else visualize_results, 
                                                        method=pupil_feature_method
                                                       )
    
    # Then, extract the eyelid features 
    eyelid_features: list[dict] = extract_eyelid_features(video, 
                                                          is_grayscale,# Do not visualize single features if we want all features  
                                                          not visualize_results if visualize_results is True else visualize_results,
                                                         )
    
    # Assert we have the same number of features for both eyelid and pupil features 
    assert len(pupil_features) == len(eyelid_features), "Frame number difference between pupil features and eyelid features"

    # Combine pupil and eyelid feature dictionaries 
    eye_features: dict[str, list[dict]] = [ {'pupil': frame_pupil_features, 'eyelids': frame_eyelid_features} 
                                            for frame_pupil_features, frame_eyelid_features
                                            in zip(pupil_features, eyelid_features)
                                          ]
    
    # Ensure the analysis was properly done (e.g. there were no silently corrupted frames not explictly caught with error by Pylids)
    assert Pi_util.inspect_video_frame_count(video) == len(eye_features), "Video framecount not equal to analyzed frames by pylids. Some frames may be corrupted"


    # Visualize features if desired 
    if(visualize_results is True):

        # If passed in a numpy array of frames, iterate 
        # over them as normal 
        if(isinstance(video, np.ndarray)):
            # If the video is not grayscale, convert to grayscale 
            if(len(video.shape) > 3):
                video = np.array([cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
                                for frame in video
                                ], 
                                dtype=np.uint8
                                )        
                

            # Iterate over the frames of the video and display the results 
            for (frame_num, (frame, frame_eye_features)) in tqdm(enumerate( zip(video, eye_features) ), desc="Visualizing"):
                # Generate a figure to display the results on 
                fig, axes = plt.subplots(1, 2)
                fig.suptitle(f"Frame: {frame_num} | Eye Features")

                # Plot the pupil features 
                visualize_pupil(frame, 
                                frame_eye_features["pupil"], 
                                method=pupil_feature_method,
                                ax=axes[0]
                            )

                # Plot the eyelid features 
                visualize_eyelids(frame, 
                                frame_eye_features["eyelids"],
                                ax=axes[1]
                                )

            # Show the plot 
            plt.show()

            # Close the figure and free the memory 
            plt.close(fig)

            return eye_features

        # Otherwise, we must iterate over the frames of the video 
        # only having one in memory at a time, so that we do not blow up our memory 
        frame_queue: object = mp.Queue(maxsize=5)

        # Initialize a second process to stream frames from the video 
        stop_event: object = mp.Event()
        frame_stream_process: object = mp.Process(target=Pi_util.destruct_video, 
                                                  args=(video, 0, float("inf"), is_grayscale, 
                                                        frame_queue, stop_event
                                                       )
                                                 )
        # Start the subprocess 
        frame_stream_process.start()

        # Parse frames from the video, gathering their features
        num_frames: int = Pi_util.inspect_video_frame_count(video)
        for frame_num in tqdm(range(num_frames), desc="Visualizing"):
            # Attempt to retrieve a frame from the frame queue 
            try:
                frame: np.ndarray | None = frame_queue.get(timeout=5)
            except queue.Empty:
                # Close the subprocess on error 
                stop_event.set()
                frame_stream_process.join()
                raise Exception("Did not receive frame in time")

            # If no frame arrived, then we are done 
            if(frame is None):
                break
            
            # Otherwise, visualize the frame and pupil features 
            # Generate a figure to display the results on 
            fig, axes = plt.subplots(1, 2)
            fig.suptitle(f"Frame: {frame_num} | Eye Features")

            # Plot the pupil features 
            visualize_pupil(frame, 
                            eye_features[frame_num]["pupil"], 
                            method=pupil_feature_method,
                            ax=axes[0]
                        )

            # Plot the eyelid features 
            visualize_eyelids(frame, 
                              eye_features[frame_num]["eyelids"],
                              ax=axes[1]
                             )
            
            # Make sure no empty space 
            plt.tight_layout()

            # Show the desired figure 
            plt.show()

            # Close the figure generated 
            plt.close(fig)

        # Join the subprocess 
        frame_stream_process.join()   

    return eye_features

def main():
    pass 

if(__name__ == "__main__"):
    main()