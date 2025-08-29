import sys 
import os
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pupil_detectors import Detector2D
from pye3d.detector_3d import CameraModel, Detector3D, DetectorMode
from typing import Literal


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
                             is_grayscale: bool=False, 
                             roi: None | np.ndarray=None
                             ) -> list[dict]:
    # First, read in the video and extract it as frames (if not already provided as frames)
    if(isinstance(video, str)):
        video: np.ndarray = Pi_util.destruct_video(video, is_grayscale=is_grayscale)

    # If video is not np.array, make it one 
    if(not isinstance(video, np.ndarray)):
        video: np.ndarray = np.array(video, dtype=np.uint8)

    # If the video is not grayscale, convert to grayscale 
    if(len(video.shape) > 3):
        video = np.array([cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) 
                          for frame in video
                         ], 
                         dtype=np.uint8
                        )
        
    # Next, splice out the ROI of the video if desired 
    if(roi is not None):
        # Format of the defined RR should be [ (BOTTOM LEFT y/x), TOP RIGHT (y/x)]
        bottom_left, top_right = roi
        
        # Splice out the desired content from all frames 
        v = v[:, bottom_left[0]:top_right[0], bottom_left[1]:top_right[1]]

    # Create 2D detector from pupil labs 
    detector_2d: object = Detector2D()
    
    # Create pye3D detector from pupil labs
    camera: object = CameraModel(focal_length=pupil_util.PUPIL_FOCAL_LENGTH, resolution=pupil_util.PUPIL_FRAME_SHAPE.tolist())
    detector_3d: object = Detector3D(camera=camera, long_term_mode=DetectorMode.blocking)

    # Initialize return value 
    eye_features_by_frame: list[dict] = []
    
    # Iterate over the frames of the video 
    for frame_idx, frame in enumerate(video):
        # Run a 2D detector on the video frame 
        result_2d = detector_2d.detect(frame)
        result_2d["timestamp"] = frame_idx / pupil_util.PUPIL_CAM_FPS

        # Pass 2D detection result to 3D detector
        result_3d = detector_3d.update_and_detect(result_2d, frame)

        # Append the eye features for this frame to the growing list 
        eye_features_by_frame.append(result_3d)        

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
                                                       estimate_eyelids=False, 
                                                       model_name=f'pytorch-{target}-v1', 
                                                       save_vid=False
                                                      )
    

    # Next, convert format from dict[str, list[dict]] to a list of dictionaries

    # First, let's find out how many frames we had 
    num_frames: int = len(pylids_out['dlc_confidence'])
    assert all(len(value) == num_frames for key, value in pylids_out.items()), "Pylids output shape mismatch. Unable to determine num frames"

    return


"""Given a path to a video or frames of a video, extract
   pupil features of the eye per frame
"""
def extract_pupil_features(video: str | np.ndarray,
                           is_grayscale: bool=False, 
                           visualize_results: bool=False, 
                           method: Literal["pupil-labs", "pylids"]="pupil-labs", 
                           roi: None | np.ndarray = None
                          ) -> list[dict]:

    # Initialize eye features variable 
    pupil_features: list[dict]

    # If desired method is pylids, then its quite easy, 
    # just call the simple analysis function 
    if(method == "pylids"):
        pupil_features = pylids_analyze_video(video, "pupil")
    # Otherwise, we want to use the pupil labs method,
    else:
        pupil_features = pupil_labs_analyze_video(video, is_grayscale, roi)

    # If we want to visualize the results 
    if(visualize_results is True):
        raise NotImplementedError("TODO: Implement visualization")
        pass 

    return pupil_features


"""Given a path to a video or frames of a video, extract
   eyelid  features of the eye per frame
"""
def extract_eyelid_features(video: str | np.ndarray,
                            visualize_results: bool=False, 
                           )-> list[dict]:
    # Extract eyelid features with pylids
    eyelid_features: dict[str, dict] = pylids_analyze_video(video, "eyelid")


    # Next, convert format from dict[str, list[dict]] to a list of dictionaries

    # Visualize results if desired 
    if(visualize_results is True):
        raise NotImplementedError("TODO: Implement roi splicing")

    return [ ]

"""Given a path to a pupil video or frames of a video, 
   extract features of the eye 
   per frame from the video using either pylids or 
   pupil labs for pupil detection, and pylids for 
   eyelid detection. 
"""
def extract_eye_features(video: str | np.ndarray, 
                         is_grayscale: bool=True,
                         visualize_results: bool=False,
                         pupil_feature_method: Literal["pupil-labs", "pylids"]="pupil-labs", 
                         pupil_roi_box: None | np.ndarray=None
                        ) -> list[dict]:
    
    # First, extract the pupil features of the video 
    print("EXTRACTING PUPIL FEATURES")
    pupil_features: list[dict] = extract_pupil_features(video, 
                                                        is_grayscale, # Do not visualize single features if we want all features  
                                                        not visualize_results if visualize_results is True else visualize_results, 
                                                        method=pupil_feature_method, 
                                                        roi=pupil_roi_box
                                                       )
    
    # Then, extract the eyelid features 
    print("EXTRACTING EYELID FEATURES")
    eyelid_features: list[dict] = extract_eyelid_features(video, # Do not visualize single features if we want all features  
                                                          not visualize_results if visualize_results is True else visualize_results,
                                                         )
    
    # Combine pupil and eyelid feature dictionaries 


    # Visualize features if desired 
    if(visualize_results is True):
        raise NotImplementedError("TODO: Implement roi splicing")


def main():
    pass 

if(__name__ == "__main__"):
    main()