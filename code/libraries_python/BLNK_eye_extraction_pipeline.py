# Import libraries 
import extract_eye_features
import dill
import scipy
import os
import numpy as np
import sys 
import subprocess
import shlex

# TODO: Make this dynamic
sys.path.append("/Users/zacharykelly/Documents/MATLAB/projects/lightLogger/raspberry_pi_firmware/utility")
import Pi_util
import matplotlib.pyplot as plt

"""Predict eye features (pupil and eyelid) for a given BLNK video"""
def predict_eye_features(filepath: str, 
                         crop_box: tuple[int], 
                         target_size: tuple[int],  
                         output_folder_path: str, 
                         temp_dir_path: str="./BLNK_temp", 
                         threshold_value: int=225
                        ) -> None:

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
                                                                         is_grayscale=True, 
                                                                         visualize_results=False, 
                                                                         pupil_feature_method='pylids', 
                                                                         safe_execution=True
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



"""
    Uses scp to copy all files from remote_dir to local_dir.
    Assumes no subdirectories and no spaces in paths. For paths with spaces,
    prefer the Paramiko version above.
"""
def scp_download_dir(
    remote_dir: str,
    local_dir: str,
    host: str,
    username: str,
    port: int = 22,
    password: str | None = None,
):
    """
    Copy *each file* from `remote_dir` on the remote host to `local_dir` locally,
    iterating file-by-file instead of using a wildcard.

    Assumptions:
      - No subdirectories
      - No spaces or odd characters in filenames (for spaces, use Paramiko/SFTP)
    Requires:
      - `sshpass` if you provide `password`
    Returns:
      (downloaded_paths, failed_filenames)
    """
    os.makedirs(local_dir, exist_ok=True)

    # 1) List plain files on the remote (filenames only, one per line)
    list_cmd = ["ssh", "-p", str(port), f"{username}@{host}",
                "find", remote_dir, "-maxdepth", "1", "-type", "f", "-printf", "%f\n"]
    if password is not None:
        list_cmd = ["sshpass", "-p", password] + list_cmd

    try:
        out = subprocess.check_output(list_cmd, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Could not list files in {remote_dir!r} on {host}: {e}") from e

    filenames = [line.strip() for line in out.splitlines() if line.strip()]

    downloaded, failed = [], []

    # 2) Copy each file individually
    for name in filenames:
        remote_path = f"{remote_dir.rstrip('/')}/{name}"
        remote_spec = f"{username}@{host}:{remote_path}"

        scp_cmd = ["scp", "-P", str(port), remote_spec, local_dir]
        if password is not None:
            scp_cmd = ["sshpass", "-p", password] + scp_cmd

        print(f"Downloading {name}...")
        try:
            subprocess.run(scp_cmd, check=True)
            downloaded.append(os.path.join(local_dir, name))
        except subprocess.CalledProcessError:
            failed.append(name)

    return downloaded, failed


def main():
    pass 

if(__name__ == "__main__"):
    main()