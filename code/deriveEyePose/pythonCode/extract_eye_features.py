"""Eye feature extraction from pupil camera video.

Provides functions for detecting and visualizing pupil and eyelid features
from eye-tracking video using Pupil Labs and pylids backends.
"""

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
import warnings

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

def generate_playable_video(video: str | np.ndarray, eye_features: list[dict],
                            output_path: str="visualized_eyefeatures.avi",
                            is_grayscale: bool=False,
                            safe_execution: bool=True,
                            keypoint_threshold: float=0.65
                           ) -> None:
    """Visualize eye features and write them into a playable video.

    Streams frames from the input video, overlays pupil and eyelid
    visualizations, and writes the annotated frames to an output video file.

    Args:
        video: Path to the video file or numpy array of frames.
        eye_features: Per-frame eye feature dictionaries containing 'pupil'
            and/or 'eyelids' keys.
        output_path: Filepath for the output visualization video.
        is_grayscale: Whether the input video is grayscale.
        safe_execution: When True, asserts that the video frame count matches
            the number of analyzed frames.
        keypoint_threshold: Minimum confidence for rendering keypoints.
    """

    # Capture the actual framecount of the video. This will be used for safeguarding 
    # Silent errors when running with safe execution 
    if(isinstance(video, str)): assert os.path.exists(video), f"Invalid path: {video}"
    actual_video_framecount: int = Pi_util.inspect_video_frame_count(video) if isinstance(video, str) else len(video)

    # Guard that we have the matching number of frames in the video and those analyzed
    if(safe_execution is True):
        assert actual_video_framecount == len(eye_features), f"Video framecount ({actual_video_framecount}) not equal to analyzed frames ({len(eye_features)}). Some frames may be corrupted"
    else:
        warnings.warn("Safe execution disabled. Skipping video frame/eye feature length check. Visualization may become out of sync")


    # Otherwise, if given a path to a video, we will stream frames, 
    # apply the visualization to the frame, and then write it to a video
  
    # Otherwise, we must iterate over the frames of the video 
    # only having one in memory at a time, so that we do not blow up our memory 
    read_queue: mp.Queue = mp.Queue(maxsize=5)
    write_queue: mp.Queue = mp.Queue(maxsize=5)
    read_stop_event: object = mp.Event()
    write_stop_event: object = mp.Event()

    # Initialize a second process to stream frames from the video 
    read_process: object = mp.Process(target=Pi_util.destruct_video, 
                                      args=(video, 0, float("inf"), is_grayscale, read_queue, read_stop_event)
                                     )
    
    # Initialize a process to stream frames to be written
    video_fps: float = Pi_util.inspect_video_FPS(video)
    write_process: object = mp.Process(target=Pi_util.frames_to_video, 
                                       args=(write_queue, output_path, video_fps, write_stop_event, 60)
                                      )
    
    # Assemble processes/stop events into lists
    processes: list[object] = [read_process, write_process]
    stop_events: list[object] = [read_stop_event, write_stop_event]

    # Start the subprocesses
    for process in processes:
        process.start() 

    # Parse frames from the video, gathering their features
    frame_num: int = 0 
    while(True):
        # Attempt to retrieve a frame from the frame queue
        try: 
            frame: np.ndarray | None = read_queue.get(timeout=10)
        except queue.Empty:
            # Close the subprocesses on error 
            for process, stop_event in zip(processes, stop_events):
                stop_event.set() 
                process.join() 

            raise Exception("Did not receive frame in time")

        # If no frame arrived, then we are done 
        if(frame is None):
            write_queue.put(None)
            break    

        # Otherwise, draw the visualizations on the image
        assert frame.dtype == np.uint8, f"Frame read with wrong dtype. Is {frame.dtype}, should be np.uint8"

        # Retrieve the eye features for this frame 
        try:
            frame_eye_features: dict = eye_features[frame_num]

        # If safe mode was disabled and there are more frames 
        # than those analyzed, simply quit now
        except IndexError as e:
            write_queue.put(None)
            
            # Close the subprocesses on error 
            for process, stop_event in zip(processes[::-1], stop_events[::-1]):
                stop_event.set() 

            break 

        # Apply the ransformations onto the frame 
        visualized_frame = frame.copy() 
        if('pupil' in frame_eye_features):
            # Visualize the pupil image on this frame 
            visualized_frame = visualize_pupil(visualized_frame, 
                                               frame_eye_features['pupil'], 
                                               keypoint_threshold=keypoint_threshold
                                              )

        if('eyelids' in frame_eye_features):
            visualized_frame = visualize_eyelids(visualized_frame, 
                                                 frame_eye_features["eyelids"]
                                                )

        # Send the frame to the writer to be written 
        assert visualized_frame.dtype == np.uint8, f"Visualized frame is wrong dtype. Is {visualized_frame.dtype}, should be np.uint8"
        write_queue.put(visualized_frame)
     
        # Increment the frame num 
        frame_num += 1 

    # Close the subprocesses 
    for process, stop_event in zip(processes, stop_events):
        if(stop_event.set() is False):
            stop_event.set() 

        process.join() 

    # Ensure the queues are empty for clean exit
    for q in (read_queue, write_queue):
        Pi_util.clear_mp_queue(q)

    return 

def pupil_labs_analyze_video(video: str | np.ndarray,
                             is_grayscale: bool=False
                             ) -> list[dict]:
    """Analyze a video with Pupil Labs detectors for pupil features.

    Runs 2D and 3D pupil detection on each frame using the Pupil Labs
    pipeline. Supports both in-memory numpy arrays and file paths with
    streaming via multiprocessing.

    Args:
        video: Path to the video file or numpy array of frames.
        is_grayscale: Whether the input video is already grayscale.

    Returns:
        List of per-frame 3D pupil detection result dictionaries.
    """
        
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

        for frame_num, frame in enumerate(video):
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
        
        # Otherwise, analyze the frame

        # Run a 2D detector on the video frame 
        result_2d = detector_2d.detect(frame)
        result_2d["timestamp"] = frame_num / pupil_util.PUPIL_CAM_FPS

        # Pass 2D detection result to 3D detector
        result_3d = detector_3d.update_and_detect(result_2d, frame)

        # Append the eye features for this frame to the growing list 
        eye_features_by_frame.append(result_3d) 

        # Increment the frame number 
        frame_num += 1

    # Join the subprocess 
    frame_stream_process.join()   

    return eye_features_by_frame

def pylids_analyze_video(video: str,
                         target: Literal["pupil", "eyelid"]
                        ) -> list[dict]:
    """Analyze a video with pylids for pupil or eyelid features.

    Args:
        video: Path to the video file. Must be a string path.
        target: Which feature to extract, either "pupil" or "eyelid".

    Returns:
        List of per-frame feature dictionaries with keypoint coordinates
        and confidence values.
    """
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

def visualize_pupil(frame: np.ndarray,
                    frame_pupil_features: dict[str, object],
                    plot_output: bool=False,
                    keypoint_threshold: float = 0.65
                   ) -> None | tuple:
    """Plot detected pupil ellipse and keypoints on a given frame.

    Draws the fitted pupil ellipse and, if available, DLC keypoints with
    confidence labels onto the frame image.

    Args:
        frame: Single video frame as a numpy array.
        frame_pupil_features: Dictionary of pupil features for this frame,
            including 'ellipse' parameters and optionally DLC keypoints.
        plot_output: If True, also generates a matplotlib figure.
        keypoint_threshold: Minimum confidence to render a keypoint in green
            rather than red.

    Returns:
        The annotated frame as a uint8 numpy array if plot_output is False,
        otherwise a tuple of (annotated_frame, fig, ax).
    """
    
    # Convert to color if not colored already 
    frame_colored: np.ndarray = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR) if frame.ndim == 2 else frame.copy()

    # Retrieve the information about the ellipse to plot 
    cx = cy = None 
    major = minor = None 
    angle: float = None

    # Gather ellipse parameters to plot the pupil
    cx, cy = tuple([ int(v) for v in frame_pupil_features['ellipse']['center']])
    major, minor = [ int(v) for v in frame_pupil_features['ellipse']['axes'] ]
    angle: float = frame_pupil_features['ellipse']['angle']

    
    # Only print valid ellipses
    if(all(axis > 0 and axis < (2 * max(frame_colored.shape[:2])) for axis in (major, minor) ) ):
        cv2.ellipse(frame_colored,
                    center=(int(round(cx)), int(round(cy))),
                    axes=(int(round(major/2)), int(round(minor/2))),  # semi-axes
                    angle=angle,                              
                    startAngle=0, endAngle=360,
                    color=(0, 0, 255), # BGR format
                    thickness=2,
                    lineType=cv2.LINE_AA
                    )
        
        # If there are key points (like if you used pupil labs), then let's also visualize them 
        if("dlc_kpts_x" in frame_pupil_features):
            # Visualize and label the keypoints
            for point_num, (x, y) in enumerate(zip(frame_pupil_features["dlc_kpts_x"], frame_pupil_features["dlc_kpts_y"])):
                confidence: int = int(frame_pupil_features["dlc_confidence"][point_num] * 100)

                # We will draw low confidence points with different colors
                color: tuple[int] = (0, 255, 0) if confidence >= keypoint_threshold else (0, 0, 255)
                cv2.circle(frame_colored, center=(int(x), int(y)), radius=3, color=color, thickness=-1, lineType=cv2.LINE_AA)
                
            # Draw their confidence measurements over the circles
            for point_num, (x, y) in enumerate(zip(frame_pupil_features["dlc_kpts_x"], frame_pupil_features["dlc_kpts_y"])):
                confidence: int = int(frame_pupil_features["dlc_confidence"][point_num] * 100)
                cv2.putText(frame_colored, str(confidence), (int(x + 6), int(y - 6)), cv2.FONT_HERSHEY_SIMPLEX, 0.25, (255, 255, 255), 1, cv2.LINE_AA)

    # Ensure we are back in uint8 space
    frame_colored = frame_colored.astype(np.uint8)

    if(plot_output is True):
        # Generate a figure ot show the plto
        fig, ax = plt.subplots()

        # Show the frame on the axis
        ax.imshow(ax.imshow(cv2.cvtColor(frame_colored, cv2.COLOR_BGR2RGB)))

        # Pretty the plot 
        ax.set_title("Pupil")

    # Return the figure if generated from here
    return frame_colored if plot_output is False else (frame_colored, fig, ax)

def extract_pupil_features(video: str | np.ndarray,
                           is_grayscale: bool=False,
                           visualize_results: bool=False,
                           method: Literal["pupil-labs", "pylids"]="pupil-labs",
                           visualization_output_filepath: str="visualized_pupilfeatures.avi",
			               safe_execution: bool=True,
                           keypoint_threshold: float=0.65
                          ) -> list[dict]:
    """Extract pupil features from a video on a per-frame basis.

    Supports two detection backends: Pupil Labs (2D/3D detector) and pylids
    (DLC-based keypoint detector). Optionally generates a visualization video
    and produces perimeter info formatted for MATLAB consumption.

    Args:
        video: Path to the video file or numpy array of frames.
        is_grayscale: Whether the input video is grayscale.
        visualize_results: If True, generates a visualization video.
        method: Detection backend to use, "pupil-labs" or "pylids".
        visualization_output_filepath: Output path for the visualization video.
        safe_execution: When True, asserts frame count consistency.
        keypoint_threshold: Minimum confidence for valid perimeter points.

    Returns:
        A tuple of (pupil_features, perimeter_info_dict) where pupil_features
        is a list of per-frame detection dictionaries and perimeter_info_dict
        contains MATLAB-compatible perimeter data.
    """

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
    if(safe_execution is True): 
        actual_video_framecount: int = Pi_util.inspect_video_frame_count(video) if isinstance(video, str) else len(video)
        assert actual_video_framecount == len(pupil_features), f"Video framecount ({actual_video_framecount}) not equal to analyzed frames ({len(pupil_features)}). Some frames may be corrupted"

    # Populate the perimeter info dict 
    # which we will use to massage data into 
    # the existing code's MATLAB readable formatb 
    perimeter_info_dict: dict = {"size": tuple([]), 
                                 "data": [], 
                                 "meta": {}
                                }
    if(method == "pylids"):
        # Retrieve the frame shape 
        perimeter_info_dict["size"] = Pi_util.inspect_video_framesize(video) if isinstance(video, str) else video.shape[1:3]

        # Reformat the data into a list of dicts with X/Y points of valid perimeter points 
        perimeter_data_arr: np.ndarray = np.empty( (len(pupil_features), 1) , dtype=object)
        for frame_num, frame_features in enumerate(pupil_features):
            # First, we will extract the X points and Y points of the pupil fit 
            Xp: np.ndarray = np.array([ x for point_num, x in enumerate(frame_features["dlc_kpts_x"]) if frame_features["dlc_confidence"][point_num] >= keypoint_threshold ]) 
            Yp: np.ndarray = np.array([ y for point_num, y in enumerate(frame_features["dlc_kpts_y"]) if frame_features["dlc_confidence"][point_num] >= keypoint_threshold ])
            confidence: np.ndarray = np.array([ confidence_value for confidence_value in frame_features["dlc_confidence"] if confidence_value >= keypoint_threshold ])

            perimeter_data_arr[frame_num, 0] = {"Xp": Xp.reshape(len(Xp), 1).astype(np.float64),
                                                "Yp": Yp.reshape(len(Yp), 1).astype(np.float64),
                                                "confidence": confidence.reshape(len(confidence), 1).astype(np.float64)
                                               }

        perimeter_info_dict["data"] = perimeter_data_arr


    # If we want to visualize the results 
    if(visualize_results is True):
        generate_playable_video(video, 
                                [ {"pupil": pupil_features[frame]} for frame in range(len(pupil_features))], 
                                visualization_output_filepath, 
                                is_grayscale=is_grayscale, 
                                safe_execution=safe_execution,
                                keypoint_threshold=keypoint_threshold
                               )

    return pupil_features, perimeter_info_dict

def visualize_eyelids(frame: np.ndarray,
                      frame_eyelid_features: dict[str, object],
                      plot_output: bool=False
                     ) -> None:
    """Plot detected upper and lower eyelid contours on a given frame.

    Draws polyline overlays for the upper and lower eyelid curves onto the
    frame image.

    Args:
        frame: Single video frame as a numpy array.
        frame_eyelid_features: Dictionary containing 'eyelid_x',
            'eyelid_lo_y', and 'eyelid_up_y' arrays.
        plot_output: If True, also generates and displays a matplotlib figure.

    Returns:
        The annotated frame as a uint8 numpy array if plot_output is False,
        otherwise a tuple of (annotated_frame, fig, ax).
    """
    

    # Convert to color if not already 
    frame_colored: np.ndarray = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR) if frame.ndim == 2 else frame.copy()

    # NOTE: For some reason, applying these back onto the original video as opposed to the cropped 
    #       video we were using with the BLNK pipeline, we had to apply offset -30 to y and -5 to x

    # Apply the visualization to the frame
    pts_lo: np.ndarray = np.column_stack((frame_eyelid_features['eyelid_x'], frame_eyelid_features['eyelid_lo_y'])).astype(np.int32)
    pts_up: np.ndarray = np.column_stack((frame_eyelid_features['eyelid_x'], frame_eyelid_features['eyelid_up_y'])).astype(np.int32)
    cv2.polylines(frame_colored, [pts_lo], isClosed=False, color=(0,0,255), thickness=2, lineType=cv2.LINE_AA)
    cv2.polylines(frame_colored, [pts_up], isClosed=False, color=(0,0,255), thickness=2, lineType=cv2.LINE_AA)

    # Ensure we are back to uint 8 
    frame_colored = frame_colored.astype(np.uint8)

    # Generate a figure if not given an axis 
    # (if we want to show output)
    if(plot_output is True):
        fig, ax = plt.subplots()

        # Show the frame on the axis
        ax.imshow(cv2.cvtColor(frame_colored, cv2.COLOR_BGR2RGB))

        # Pretty the plot 
        ax.set_title("Eyelids")

        # Show the plot if we generated the ax here 
        plt.show()

    # Return the figure if generated from here
    return frame_colored if plot_output is False else (frame_colored, fig, ax)

def extract_eyelid_features(video: str | np.ndarray,
                            is_grayscale: bool=False,
                            visualize_results: bool=False,
                            visualization_output_filepath: str="visualized_eyelidfeatures.avi",
			                safe_execution: bool=True,
                            keypoint_threshold: float=0.65
                           )-> list[dict]:
    """Extract eyelid features from a video on a per-frame basis.

    Uses the pylids backend to detect upper and lower eyelid contours.
    Optionally generates a visualization video of the detected eyelids.

    Args:
        video: Path to the video file or numpy array of frames.
        is_grayscale: Whether the input video is grayscale.
        visualize_results: If True, generates a visualization video.
        visualization_output_filepath: Output path for the visualization video.
        safe_execution: When True, asserts frame count consistency.
        keypoint_threshold: Minimum confidence for rendering keypoints in
            the visualization.

    Returns:
        List of per-frame eyelid feature dictionaries.
    """
    # Extract eyelid features with pylids
    eyelid_features: dict[str, dict] = pylids_analyze_video(video, "eyelid")

    # Check to ensure that the video is well formed
    if(safe_execution is True): 
        actual_video_framecount: int = Pi_util.inspect_video_frame_count(video) if isinstance(video, str) else len(video)
        assert actual_video_framecount == len(eyelid_features), f"Video framecount ({actual_video_framecount}) not equal to analyzed frames ({len(eyelid_features)}). Some frames may be corrupted"

    # Visualize results if desired 
    if(visualize_results is True):
        generate_playable_video(video, 
                                [ {"eyelids": eyelid_features[frame]} for frame in range(len(eyelid_features))], 
                                visualization_output_filepath, 
                                is_grayscale=is_grayscale, 
                                safe_execution=safe_execution,
                                keypoint_threshold=keypoint_threshold
                               )

    return eyelid_features

def extract_eye_features(video: str | np.ndarray,
                         is_grayscale: bool=True,
                         visualize_results: bool=False,
                         visualization_output_filepath: str="visualized_eyefeatures.avi",
			             safe_execution: bool=True,
                         pupil_feature_method: Literal["pylids", "pupil-labs"] = "pylids",
                         keypoint_threshold: float=0.65
                        ) -> list[dict]:
    """Extract combined pupil and eyelid features from a pupil camera video.

    Runs pupil detection (via pylids or Pupil Labs) and eyelid detection
    (via pylids) on each frame, then merges the results into a single
    per-frame feature dictionary. Optionally generates a visualization video.

    Args:
        video: Path to the video file or numpy array of frames.
        is_grayscale: Whether the input video is grayscale.
        visualize_results: If True, generates a visualization video.
        visualization_output_filepath: Output path for the visualization video.
        safe_execution: When True, asserts frame count consistency.
        pupil_feature_method: Backend for pupil detection, "pylids" or
            "pupil-labs".
        keypoint_threshold: Minimum confidence for valid keypoints.

    Returns:
        A tuple of (eye_features, perimeter_info_dict) where eye_features
        is a list of per-frame dictionaries with 'pupil' and 'eyelids' keys,
        and perimeter_info_dict contains MATLAB-compatible perimeter data.
    """
    
    # Retrieve the framecount of the video. We will use this for safeguarding 
    # silent errors when running with safe execution
    actual_video_framecount: int = Pi_util.inspect_video_frame_count(video) if isinstance(video, str) else len(video)

    # First, extract the pupil features of the video 
    pupil_features, perimeter_info_dict = extract_pupil_features(video, 
                                                                 is_grayscale, # Do not visualize single features if we want all features  
                                                                 not visualize_results if visualize_results is True else visualize_results, 
							                                     safe_execution=safe_execution,
                                                                 method=pupil_feature_method,
                                                                 keypoint_threshold=keypoint_threshold
                                                                )
    
    # Ensure the analysis was properly done (e.g. there were no silently corrupted frames not explictly caught with error by Pylids)
    if(safe_execution is True): 
        assert actual_video_framecount == len(pupil_features), f"Video framecount ({actual_video_framecount}) not equal to analyzed frames ({len(pupil_features)}). Some frames may be corrupted"

    # Then, extract the eyelid features 
    eyelid_features: list[dict] = extract_eyelid_features(video, 
                                                          is_grayscale,# Do not visualize single features if we want all features  
                                                          not visualize_results if visualize_results is True else visualize_results,
							                              safe_execution=safe_execution,
                                                          keypoint_threshold=keypoint_threshold
                                                         )
    
     # Ensure the analysis was properly done (e.g. there were no silently corrupted frames not explictly caught with error by Pylids)
    if(safe_execution is True):
        assert actual_video_framecount == len(eyelid_features), f"Video framecount ({actual_video_framecount}) not equal to analyzed frames ({len(eyelid_features)}). Some frames may be corrupted"


    # Assert we have the same number of features for both eyelid and pupil features 
    assert len(pupil_features) == len(eyelid_features), "Frame number difference between pupil features and eyelid features"

    # Combine pupil and eyelid feature dictionaries 
    eye_features: dict[str, list[dict]] = [ {'pupil': frame_pupil_features, 'eyelids': frame_eyelid_features} 
                                            for frame_pupil_features, frame_eyelid_features
                                            in zip(pupil_features, eyelid_features)
                                          ]

    # Visualize features if desired 
    if(visualize_results is True):
        generate_playable_video(video, eye_features, 
                                output_path=visualization_output_filepath,
                                is_grayscale=is_grayscale,
                                safe_execution=safe_execution,
                                keypoint_threshold=keypoint_threshold
                               )

    return eye_features, perimeter_info_dict

def main():
    """Entry point for command-line execution."""
    pass

if(__name__ == "__main__"):
    main()
