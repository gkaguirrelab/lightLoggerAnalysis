import os
import cv2
import numpy as np
from scipy.io import savemat

video_path = "/Users/sophiamirabal/Library/CloudStorage/Dropbox-Aguirre-BrainardLab/Sophia Mirabal/temp_for_sharing/sophia_fielding_function.avi"

assert os.path.exists(video_path), f"Video not found: {video_path}"

cap = cv2.VideoCapture(video_path)
assert cap.isOpened(), "Could not open video"

fps = cap.get(cv2.CAP_PROP_FPS)
print("FPS:", fps)

times_sec = [40, 41, 94, 165, 187, 64, 78, 79, 179, 182, 74, 76, 77, 135, 200, 53, 55, 171, 172, 195]

frames_idx = [round(t * fps) for t in times_sec]
print("Frame indices:", frames_idx)

frames = []

for idx in frames_idx:
    cap.set(cv2.CAP_PROP_POS_FRAMES, idx)
    ret, frame = cap.read()

    if not ret:
        raise RuntimeError(f"Could not read frame {idx}")

    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    frames.append(gray)

cap.release()

frames = np.array(frames, dtype=np.uint8)

savemat("selected_twilight_frames.mat", {
    "frames": frames,
    "times_sec": np.array(times_sec),
    "frames_idx": np.array(frames_idx),
    "fps": fps
})

print("Saved selected_twilight_frames.mat")
print("frames shape:", frames.shape)