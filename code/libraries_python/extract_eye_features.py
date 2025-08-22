import sys 
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pupil_detectors import Detector2D
from pye3d.detector_3d import CameraModel, Detector3D, DetectorMode


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

"""Given a path to a pupil video, extract features of the eye 
   per frame from the video using pupil-labs library
"""
def extract_eye_features(video: str | np.ndarray, 
                         is_grayscale: bool=True,
                         visualize_results: bool=False,
                        ) -> list[dict]:
    # First, read in the video and extract it as frames (if not already provided as frames)
    if(isinstance(video, str)):
        video: np.ndarray = Pi_util.destruct_video(video, is_grayscale=is_grayscale)

    # If video is not np.array, make it one 
    if(not isinstance(video, np.ndarray)):
        video: np.ndarray = np.array(video, dtype=np.uint8)
        
    # If the video is not grayscale, convert to grayscale 
    if(len(video.shape) > 3):
        pass 

    # Create 2D detector from pupil labs 
    detector_2d: object = Detector2D()
    
    # Create pye3D detector from pupil labs
    camera: object = CameraModel(focal_length=pupil_util.PUPIL_FOCAL_LENGTH, resolution=pupil_util.PUPIL_FRAME_SHAPE.tolist())
    detector_3d: object = Detector3D(camera=camera, long_term_mode=DetectorMode.blocking)

    # Initialize return value 
    eye_features_by_frame: list[dict] = []
    
    # Iterate over the frames of the video 
    for frame_idx, frame in enumerate(video):
        print(f"Processing frame: {frame_idx+1}/{len(video)}")

        # Run a 2D detector on the video frame 
        result_2d = detector_2d.detect(frame)
        result_2d["timestamp"] = frame_idx / pupil_util.PUPIL_CAM_FPS

        # Pass 2D detection result to 3D detector
        result_3d = detector_3d.update_and_detect(result_2d, frame)
        ellipse_3d = result_3d["ellipse"]

        # Generate a figure to plot the result 
        if(visualize_results):
            fig, ax = plt.subplots()

            # Plot the resulting image onto the axis 
            plt.imshow(frame, cmap='gray')

            # Create an ellipse patch
            ellipse_patch = patches.Ellipse(    
                xy=ellipse_3d["center"],  # (x, y)
                width=ellipse_3d["axes"][0],
                height=ellipse_3d["axes"][1],
                angle=ellipse_3d["angle"],
                edgecolor="lime",   # green outline
                facecolor="none",   # transparent fill
                linewidth=2
            )

            # Add ellipse to the axes
            ax.add_patch(ellipse_patch)

            # Show the result
            plt.show()

            # Clear the old figure 
            plt.close(fig)

        # Append the eye features for this frame to the growing list 
        eye_features_by_frame.append(result_3d)        

    return eye_features_by_frame


def main():
    pass 

if(__name__ == "__main__"):
    main()