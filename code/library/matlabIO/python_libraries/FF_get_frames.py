import numpy as np
from scipy.io import savemat

from video_io import inspect_video_FPS, extract_frames_from_video

video_path = "/Users/sophiamirabal/Library/CloudStorage/Dropbox-Aguirre-BrainardLab/Sophia Mirabal/temp_for_sharing/sophia_fielding_function.avi"

fps = inspect_video_FPS(video_path)
print("FPS:", fps)

times_sec = [54, 70, 88, 150]
frames_idx = [round(t * fps) for t in times_sec]
print("Frame indices:", frames_idx)

frames = extract_frames_from_video(video_path, frames_idx, is_grayscale=True)

savemat("selected_twilight_frames.mat", {
    "frames": frames,
    "times_sec": np.array(times_sec),
    "frames_idx": np.array(frames_idx),
    "fps": fps
})