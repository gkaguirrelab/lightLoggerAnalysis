# Import libraries 
import extract_eye_features
import dill
import scipy
import os
import numpy as np
import sys 
import subprocess
import shlex
from typing import Literal

light_logger_dir_path: str = os.path.expanduser("~/Documents/MATLAB/projects/lightLogger")
Pi_util_dir_path: str = os.path.join(light_logger_dir_path, "raspberry_pi_firmware", "utility")
sys.path.append(Pi_util_dir_path)
import Pi_util
import matplotlib.pyplot as plt

"""Use this command to copy data off of the linux laptop to the analysis machine"""
"""
while true; do
  /opt/homebrew/bin/rsync -av --partial --append-verify --info=progress2 --timeout=60 --compress-level=1 \
    --bwlimit=60000 \
    -e 'ssh -p 22 -o ServerAliveInterval=30 -o ServerAliveCountMax=10 -o TCPKeepAlive=yes -o IPQoS=throughput -o RekeyLimit=4G -o Compression=no' \
    gka@10.30.8.122:/media/gka/EYEVIDEOS/lightLevel/HERO_gka/  ./HERO_gka/ && break
  echo "Retrying in 10s…"; sleep 10
done

"""

"""Predict eye features (pupil and eyelid) for a given BLNK video"""
def predict_eye_features(filepath: str,
                         output_folder_path: str, 
                         crop_box: tuple[int]=(140, 275, 190, 425), 
                         target_size: tuple[int]=(480, 640),   
                         temp_dir_path: str="./BLNK_temp", 
                         threshold_value: int=225,
                         visualize_results: bool=False,
                         is_grayscale: bool=False,
                         safe_execution: bool=True,
                         pupil_feature_method: Literal["pylids", "pupil-labs"] = "pylids"
                        ) -> list[dict]:

    # Define the portion of the video to crop out 
    t, b, l, r = crop_box

    # Extract the name of this video
    video_name: str = os.path.splitext(os.path.basename(filepath.rstrip('/')))[0]

    # Find the FPS of the video 
    video_fps: float = Pi_util.inspect_video_FPS(filepath)

    # Videos are small, so we can load them entirely in from memory 
    video_as_arr: np.ndarray = Pi_util.destruct_video(filepath, is_grayscale=True)

    # Crop out only the eye from the video
    # and set the rest of the frame to black 
    video_cropped: np.ndarray = video_as_arr[:, t:b, l:r].copy()
    white_pixels: np.ndarray = np.mean(video_cropped, axis=(0,))  > threshold_value
    video_cropped[:, white_pixels] = 0

    # Resize the video to not a small resolution 
    y_offset = (target_size[0] - video_cropped.shape[1]) // 2
    x_offset = (target_size[1] - video_cropped.shape[2]) // 2

    video_resized: np.ndarray = np.zeros((len(video_cropped), *target_size), dtype=np.uint8)
    video_resized[:, y_offset:y_offset + video_cropped.shape[1], x_offset:x_offset + video_cropped.shape[2]] = video_cropped

    # Generate a temp video from this cropped video 
    if(not os.path.exists(temp_dir_path)):
        os.mkdir(temp_dir_path)

    temp_video_path: str = os.path.join(temp_dir_path, f"temp_{video_name}.avi")
    Pi_util.frames_to_video(video_resized, temp_video_path, video_fps)

    # Extract the eye features for this video
    eye_features: list[dict] = extract_eye_features.extract_eye_features(temp_video_path, 
                                                                         is_grayscale=is_grayscale, 
                                                                         visualize_results=visualize_results, 
                                                                         pupil_feature_method=pupil_feature_method, 
                                                                         safe_execution=safe_execution
                                                                        )

    # Repackage the features along with their metadata
    eye_features_dict: dict = {}
    eye_features_dict["eye_features"] = eye_features
    eye_features_dict["metadata"] = {'threshold_value': {'v': threshold_value, 
                                                    'desc': "binary mask constructed from avg cropped frame. Pixels above this value=0. Done to remove light around the eye"
                                                   },
                                'crop_box': {'v': crop_box, 
                                             'desc': "box cropped out from original video to isolate the eye (t, b, l, r)"
                                            },
                                'target_size': {'v': target_size,
                                                'desc': "target size after crop + threshold. Eye video is centered via padding to reach this size"},
                                'model_names': {'v': ('pytorch-pupil-v1', 'pytorch-eyelid-v1'), 
                                                'desc': "models used for pupil/eyelid fitting"
                                               }
                               }
    
    
    # Output the features
    scipy.io.savemat(os.path.join(output_folder_path, f"{video_name}_eye_features.mat"), 
                    {"eye_features": eye_features_dict}
                    )
    
    # Remvove the temp avi video 
    os.remove(temp_video_path)

    return eye_features_dict


"""Predict eyelid features for a given BLNK video"""
def predict_eyelid_features(filepath: str,
                            output_folder_path: str, 
                            crop_box: tuple[int]=(140, 275, 190, 425), 
                            target_size: tuple[int]=(480, 640),   
                            temp_dir_path: str="./BLNK_temp", 
                            threshold_value: int=225,
                            visualize_results: bool=False,
                            is_grayscale: bool=False,
                            safe_execution: bool=True,
                          ) -> list[dict]:

    # Define the portion of the video to crop out 
    t, b, l, r = crop_box

    # Extract the name of this video
    video_name: str = os.path.splitext(os.path.basename(filepath.rstrip('/')))[0]

    # Find the FPS of the video 
    video_fps: float = Pi_util.inspect_video_FPS(filepath)

    # Videos are small, so we can load them entirely in from memory 
    video_as_arr: np.ndarray = Pi_util.destruct_video(filepath, is_grayscale=True)

    # Crop out only the eye from the video
    # and set the rest of the frame to black 
    video_cropped: np.ndarray = video_as_arr[:, t:b, l:r].copy()
    white_pixels: np.ndarray = np.mean(video_cropped, axis=(0,))  > threshold_value
    video_cropped[:, white_pixels] = 0

    # Resize the video to not a small resolution 
    y_offset = (target_size[0] - video_cropped.shape[1]) // 2
    x_offset = (target_size[1] - video_cropped.shape[2]) // 2

    video_resized: np.ndarray = np.zeros((len(video_cropped), *target_size), dtype=np.uint8)
    video_resized[:, y_offset:y_offset + video_cropped.shape[1], x_offset:x_offset + video_cropped.shape[2]] = video_cropped

    # Generate a temp video from this cropped video 
    if(not os.path.exists(temp_dir_path)):
        os.mkdir(temp_dir_path)

    temp_video_path: str = os.path.join(temp_dir_path, f"temp_{video_name}.avi")
    Pi_util.frames_to_video(video_resized, temp_video_path, video_fps)

    # Extract the eye features for this video
    eyelid_features: list[dict] = extract_eye_features.extract_eyelid_features(temp_video_path, 
                                                                            is_grayscale=is_grayscale, 
                                                                            visualize_results=visualize_results,
                                                                            safe_execution=safe_execution
                                                                          )

    # Repackage the features along with their metadata
    eyelid_features_dict: dict = {}
    eyelid_features_dict["eyelid_features"] = eyelid_features
    eyelid_features_dict["metadata"] = {'threshold_value': {'v': threshold_value, 
                                                    'desc': "binary mask constructed from avg cropped frame. Pixels above this value=0. Done to remove light around the eye"
                                                   },
                                      'crop_box': {'v': crop_box, 
                                             'desc': "box cropped out from original video to isolate the eye (t, b, l, r)"
                                            },
                                      'target_size': {'v': target_size,
                                                'desc': "target size after crop + threshold. Eye video is centered via padding to reach this size"},
                                      'model_names': {'v': ('pytorch-pupil-v1', 'pytorch-eyelid-v1'), 
                                                'desc': "models used for pupil/eyelid fitting"
                                               }
                               }
    
    
    # Output the features
    scipy.io.savemat(os.path.join(output_folder_path, f"{video_name}_eyelid_features.mat"), 
                    {"eyelid_features": eyelid_features_dict}
                    )
    
    # Remvove the temp avi video 
    os.remove(temp_video_path)

    return eyelid_features_dict

def main():
    pass 

if(__name__ == "__main__"):
    main()